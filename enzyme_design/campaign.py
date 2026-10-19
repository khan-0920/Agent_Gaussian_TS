"""Bounded, resumable model execution with an immutable input snapshot."""
import contextlib
import copy
import datetime
import fcntl
import os
import signal
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path
from .io import read_json, write_json, sha256, within
from .spec import STAGES, validate_spec, validate_profile


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def process_identity(pid):
    """Linux birth time prevents treating a recycled PID as the model process."""
    try:
        # comm can contain spaces and parentheses; fields after the final ')' start at field 3.
        return Path('/proc/%d/stat' % pid).read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, OSError, IndexError):
        return None


def process_alive(attempt):
    if attempt.get('host') != socket.gethostname():
        raise ValueError('running job belongs to another host; inspect it there before recovery')
    pid = attempt.get('pid')
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    recorded = attempt.get('process_identity')
    current = process_identity(pid)
    return recorded is None or current is None or recorded == current


class Campaign:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.state_path = self.root / 'campaign.json'

    @classmethod
    def create(cls, config_path, profile_path, run_dir):
        config_path, profile_path = Path(config_path).resolve(), Path(profile_path).resolve()
        spec, profile = validate_spec(read_json(config_path)), validate_profile(read_json(profile_path))
        root = Path(run_dir).resolve()
        # Validate readable sources before creating the run directory.
        sources = [('theozyme', spec['theozyme'])]
        for index, ligand in enumerate(spec['ligands']):
            for field in ('sdf', 'rosetta_params'):
                if field in ligand:
                    sources.append(('ligands.%d.%s' % (index, field), ligand[field]))
        for index, report in enumerate(spec['evidence'].get('reports', [])):
            sources.append(('evidence.reports.%d' % index, report))
        originals = []
        for key, raw in sources:
            path = (config_path.parent / raw).resolve()
            if not path.is_file():
                raise ValueError('missing input: ' + key + ': ' + str(path))
            originals.append((key, path))
        # Keep virtualenv interpreter symlinks: resolving them would select the base environment.
        for tool in profile['tools'].values():
            for key in ('python', 'executable', 'checkpoint', 'packer_checkpoint', 'repo'):
                if key in tool:
                    tool[key] = os.path.abspath(profile_path.parent / tool[key])
        root.mkdir(parents=True, exist_ok=False)
        (root / 'inputs').mkdir()
        files = {}
        for index, (key, source) in enumerate(originals):
            relative = 'inputs/%03d_%s' % (index, source.name)
            target = root / relative
            target.write_bytes(source.read_bytes())
            files[relative] = {'sha256': sha256(target), 'source': str(source), 'role': key}
            if key == 'theozyme':
                spec['theozyme'] = relative
            elif key.startswith('ligands.'):
                _, number, field = key.split('.')
                spec['ligands'][int(number)][field] = relative
            else:
                spec['evidence']['reports'][int(key.rsplit('.', 1)[1])] = relative
        write_json(root / 'spec.json', spec)
        write_json(root / 'profile.json', profile)
        for name in ('spec.json', 'profile.json'):
            files[name] = {'sha256': sha256(root / name), 'role': 'configuration'}
        state = {'schema_version': 1, 'name': spec['name'], 'created_at': utc(), 'updated_at': utc(),
                 'input_files': files, 'stages': {stage: {'status': 'pending', 'attempts': []} for stage in STAGES},
                 'usage': {'attempts': 0, 'seconds': 0.0}, 'history': [], 'paused': False}
        write_json(root / 'campaign.json', state)
        return cls(root)

    @contextlib.contextmanager
    def lock(self):
        with (self.root / '.lock').open('a+') as stream:
            try:
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise ValueError('campaign is already being modified by another process') from None
            try:
                yield
            finally:
                fcntl.flock(stream, fcntl.LOCK_UN)

    def load(self):
        state = read_json(self.state_path)
        for relative, record in state['input_files'].items():
            path = within(self.root, relative)
            if not path.is_file() or sha256(path) != record['sha256']:
                raise ValueError('campaign input changed: ' + relative + '; start a new campaign for changed inputs')
        return state, validate_spec(read_json(self.root / 'spec.json')), validate_profile(read_json(self.root / 'profile.json'))

    def save(self, state):
        state['updated_at'] = utc()
        write_json(self.state_path, state)

    def verified_result(self, stage, state):
        record = state['stages'][stage]
        if record['status'] != 'completed':
            raise ValueError('dependency not completed: ' + stage)
        attempt = record['attempts'][-1]
        for relative, digest in attempt['artifacts'].items():
            path = within(self.root, relative)
            if not path.is_file() or sha256(path) != digest:
                raise ValueError('stage artifact changed or missing: ' + relative)
        return read_json(within(self.root, attempt['result']))

    def context(self):
        state, spec, _ = self.load()
        stages = {}
        for name, record in state['stages'].items():
            stages[name] = {'status': record['status'], 'attempts': len(record['attempts'])}
            if record['attempts']:
                last = record['attempts'][-1]
                stages[name].update({key: last.get(key) for key in ('pid', 'log', 'error', 'seconds', 'result')})
            if record['status'] == 'completed':
                result = self.verified_result(name, state)
                stages[name].update(candidate_count=len(result.get('candidates', [])), summary=result.get('summary'))
        next_stage = next((name for name in STAGES if stages[name]['status'] != 'completed'), None)
        return {'name': state['name'], 'run_dir': str(self.root), 'paused': state['paused'],
                'stages': stages, 'next_stage': next_stage, 'usage': state['usage'], 'limits': spec['limits'],
                'evidence': spec['evidence'], 'history': state['history'][-20:],
                'status': 'COMPUTATION_COMPLETE' if next_stage is None else 'IN_PROGRESS',
                'interpretation': 'Design and predicted geometry are computational evidence, not experimental enzyme activity.'}

    def doctor(self):
        _, _, profile = self.load()
        checks = []
        for stage, tool in profile['tools'].items():
            for key in ('python', 'executable', 'checkpoint', 'packer_checkpoint', 'repo'):
                if key not in tool:
                    continue
                path = Path(tool[key])
                ok = path.is_dir() if key == 'repo' else path.is_file()
                if key in ('python', 'executable'):
                    ok = ok and os.access(path, os.X_OK)
                checks.append({'stage': stage, 'kind': key, 'path': str(path), 'available': ok})
        return {'ready': all(check['available'] for check in checks), 'checks': checks,
                'scope': 'Paths/executability only; model imports, GPU, weights and actual output are checked by execution.'}

    def run_stage(self, stage, reason, retry=False):
        if stage not in STAGES:
            raise ValueError('unknown stage')
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError('a research or engineering reason is required')
        with self.lock():
            state, spec, profile = self.load()
            if state['paused']:
                raise ValueError('campaign is paused')
            predecessors = STAGES[:STAGES.index(stage)]
            previous, dependencies = None, {}
            for parent in predecessors:
                previous = self.verified_result(parent, state)
                dependencies[parent] = previous
            record = state['stages'][stage]
            if record['status'] == 'completed':
                return self.verified_result(stage, state)
            if record['status'] == 'running':
                attempt = record['attempts'][-1]
                if process_alive(attempt):
                    raise ValueError('previous worker is still running; do not submit a duplicate')
                record['status'] = attempt['status'] = 'interrupted'
                # The reserved timeout remains charged after a lost controller.
                attempt['error'] = 'worker no longer live; prior reservation retained'
                self.save(state)
            if record['status'] in ('failed', 'interrupted') and not retry:
                raise ValueError('inspect failure and explicitly request retry before another attempt')
            limits = spec['limits']
            reservation = min(limits['stage_timeout_seconds'], limits['total_seconds'] - state['usage']['seconds'])
            if state['usage']['attempts'] >= limits['max_attempts'] or reservation < 1:
                raise ValueError('campaign execution budget exhausted')
            index = len(record['attempts']) + 1
            relative = 'steps/%s/%04d' % (stage, index)
            directory = self.root / relative
            directory.mkdir(parents=True, exist_ok=False)
            runtime = directory / 'runtime' / 'enzyme_design'
            runtime.mkdir(parents=True)
            for source in Path(__file__).parent.glob('*.py'):
                shutil.copy2(source, runtime / source.name)
            request = {'schema_version': 1, 'root': str(self.root), 'stage': stage, 'attempt_dir': relative,
                       'spec': spec, 'profile': profile, 'previous': previous, 'dependencies': dependencies}
            write_json(directory / 'request.json', request)
            command = [profile['tools'][stage]['python'], '-m', 'enzyme_design.worker', str(directory / 'request.json')]
            attempt = {'index': index, 'status': 'running', 'started_at': utc(), 'host': socket.gethostname(),
                       'reason': reason, 'command': command, 'request_sha256': sha256(directory / 'request.json'),
                       'reserved_seconds': reservation, 'log': relative + '/worker.log', 'pid': None}
            record['status'] = 'running'
            record['attempts'].append(attempt)
            state['usage']['attempts'] += 1
            state['usage']['seconds'] += reservation
            state['history'].append({'stage': stage, 'attempt': index, 'reason': reason, 'at': utc()})
            self.save(state)
            environment = os.environ.copy()
            environment.update(profile['environment'])
            environment.update(CUDA_VISIBLE_DEVICES=str(profile['gpu']), OMP_NUM_THREADS=str(profile['cpus']),
                               MKL_NUM_THREADS=str(profile['cpus']), OPENBLAS_NUM_THREADS=str(profile['cpus']),
                               PYTHONDONTWRITEBYTECODE='1', PYTHONPATH=str(runtime.parent))
            started = time.monotonic()
            process = None
            interrupted = False
            try:
                with (directory / 'worker.log').open('w') as log:
                    process = subprocess.Popen(command, cwd=directory, env=environment, stdout=log,
                                               stderr=subprocess.STDOUT, start_new_session=True)
                    attempt['pid'] = process.pid
                    attempt['process_identity'] = process_identity(process.pid)
                    self.save(state)
                    returncode = process.wait(timeout=reservation)
                attempt['returncode'] = returncode
                if returncode:
                    raise ValueError('worker exited with code %d; inspect %s' % (returncode, attempt['log']))
                result = read_json(directory / 'result.json')
                if result.get('stage') != stage or result.get('schema_version') != 1:
                    raise ValueError('worker result schema/stage mismatch')
                if not isinstance(result.get('candidates'), list) or not isinstance(result.get('artifacts'), list):
                    raise ValueError('worker result lacks candidates or artifacts')
                artifacts = {}
                # Include all raw artifacts, not just the worker's preferred subset.
                for path in directory.rglob('*'):
                    if path.is_file():
                        rel = str(path.relative_to(self.root))
                        within(self.root, rel)
                        artifacts[rel] = sha256(path)
                for rel in result['artifacts']:
                    if rel not in artifacts:
                        raise ValueError('worker declared missing or foreign artifact: ' + str(rel))
                for candidate in result['candidates']:
                    if not isinstance(candidate, dict) or not isinstance(candidate.get('id'), str):
                        raise ValueError('invalid candidate record')
                    for field in ('structure',):
                        if field in candidate:
                            path = within(self.root, candidate[field])
                            if not path.is_file():
                                raise ValueError('candidate structure missing')
                if len({c['id'] for c in result['candidates']}) != len(result['candidates']):
                    raise ValueError('duplicate candidate IDs')
                attempt['result'] = relative + '/result.json'
                attempt['artifacts'] = artifacts
                record['status'] = attempt['status'] = 'completed'
            except BaseException as error:
                # Cleanup covers timeout, cancellation and parser failures; never kill unrelated jobs.
                interrupted = isinstance(error, KeyboardInterrupt)
                if process is not None and process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                record['status'] = attempt['status'] = 'interrupted' if interrupted else 'failed'
                attempt['error'] = type(error).__name__ + ': ' + str(error)
                raise
            finally:
                elapsed = time.monotonic() - started
                attempt['seconds'] = elapsed
                attempt['finished_at'] = utc()
                state['usage']['seconds'] += elapsed - reservation
                self.save(state)
            return result

    def run(self, reason, retry=False):
        for stage in STAGES:
            self.run_stage(stage, reason, retry=retry)
        return self.context()

    def set_paused(self, value, reason):
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError('reason is required')
        with self.lock():
            state, _, _ = self.load()
            state['paused'] = bool(value)
            state['history'].append({'action': 'pause' if value else 'resume', 'reason': reason, 'at': utc()})
            self.save(state)
        return self.context()

    def report(self):
        state, _, _ = self.load()
        if state['stages']['select']['status'] != 'completed':
            return {'status': 'INCOMPLETE', 'context': self.context()}
        for stage in STAGES:
            result = self.verified_result(stage, state)
        return copy.deepcopy(result)
