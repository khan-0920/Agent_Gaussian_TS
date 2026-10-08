"""Auditable non-linear tool execution; only this engine constructs evidence."""
import copy
import hashlib
import json
import time
from dataclasses import asdict
from pathlib import Path
from .models import Structure, Reaction, Limits, finite
from .geometry import edit_geometry, inspect_clashes, distance
from .gaussian import JobSpec, prepare_input, TS_KINDS, OPT_KINDS
from .parser import parse_log
from .runner import LocalRunner
from .store import Store
from .validation import validate_ts, same_geometry

GEOMETRY_TOOLS = ('set_distance', 'set_angle', 'set_dihedral', 'translate_fragment', 'rotate_fragment')
TOOL_PARAMETERS = {
    'inspect_geometry': [], 'use_structure': ['structure_id'], 'undo_geometry_edit': [],
    'save_hypothesis': ['hypothesis'], 'request_human_review': ['question'],
    'verify_ts': ['ts_job_id'], 'pause': [], 'resume': [],
    'run_gaussian': list(JobSpec.__dataclass_fields__) + ['structure_id', 'parent_job', 'timeout_seconds'],
}


class Session:
    def __init__(self, store, state, runner=None):
        self.store, self.state = store, state
        self.runner = runner or LocalRunner(state['gaussian_command'])

    @classmethod
    def create(cls, root, reaction, limits=None, command=None, job_defaults=None, problem_definition=None):
        limits = limits or Limits()
        command = command or ['g16']
        job_defaults = job_defaults or {}
        if not isinstance(job_defaults, dict) or set(job_defaults) - (set(JobSpec.__dataclass_fields__) - {'kind', 'constraints'}):
            raise ValueError('invalid Gaussian defaults')
        JobSpec('ts', **job_defaults)
        LocalRunner(command)  # Validate argv before creating artifacts.
        root = Path(root)
        root.mkdir(parents=True, exist_ok=False)
        store = Store(root)
        for name in ('actions', 'structures', 'jobs'):
            store.path(name).mkdir()
        state = {'version': 1, 'reaction': reaction.to_dict(), 'limits': asdict(limits), 'gaussian_command': command, 'job_defaults': job_defaults,
                 'problem_definition': problem_definition,
                 'structures': {'reactant': {'data': reaction.reactant.to_dict(), 'source_job': None},
                                'product': {'data': reaction.product.to_dict(), 'source_job': None}},
                 'current_structure_id': 'reactant', 'undo': [], 'jobs': [], 'fingerprints': [], 'hypotheses': [],
                 'usage': {'actions': 0, 'jobs': 0, 'wall_seconds': 0., 'core_seconds': 0.},
                 'paused': False, 'pending_action': None, 'last_error': None, 'report': None,
                 'history': [], 'created_at': time.time()}
        for key, s in [('reactant', reaction.reactant), ('product', reaction.product)]:
            store.write_text('structures/%s.xyz' % key, s.to_xyz())
        store.write_json('state.json', state)
        return cls(store, state)

    @classmethod
    def open(cls, root, runner=None):
        store = Store(root)
        with store.lock():
            state = store.read_json('state.json')
            if state.get('version') != 1:
                raise ValueError('unsupported checkpoint version')
            Reaction.from_dict(state['reaction'])
            Limits(**state['limits'])
            session = cls(store, state, runner)
            session._recover()
            return session

    def _recover(self):
        # Recover artifacts left between directory creation and a durable checkpoint.
        changed = False
        for directory in sorted(self.store.path('actions').iterdir()):
            if not directory.name.isdigit() or not directory.is_dir():
                continue
            sequence = int(directory.name)
            if sequence <= self.state['usage']['actions']:
                continue
            audit_path = 'actions/%s/action.json' % directory.name
            audit = self.store.read_json(audit_path) if self.store.path(audit_path).exists() else {'id': directory.name}
            audit.update(status='interrupted', error='orphaned action journal recovered; no automatic execution')
            self.store.write_json(audit_path, audit)
            self.state['usage']['actions'] = sequence
            changed = True
        pending = self.state.get('pending_action')
        if not pending:
            if changed:
                self.store.write_json('state.json', self.state)
            return
        # Reserved resource allowance is retained: a lost process is not assumed free.
        job_id = pending.get('job_id')
        for j in self.state['jobs']:
            if j['id'] == job_id and j['status'] == 'running':
                j['status'] = 'interrupted'
        path = 'actions/%s/action.json' % pending['id']
        audit = self.store.read_json(path) if self.store.path(path).exists() else pending.get('audit', {'id': pending['id']})
        audit.update(status='interrupted', error='interrupted action; reservation retained; no blind resubmission')
        self.store.write_json(path, audit)
        self.state['last_error'] = 'interrupted action; inspect raw log and Gaussian processes before continuing'
        self.state['pending_action'] = None
        self.store.write_json('state.json', self.state)

    def observe(self):
        with self.store.lock():
            self.state = self.store.read_json('state.json')
            self._recover()
            return self._observation()

    def _observation(self):
        state = copy.deepcopy(self.state)
        sid = state['current_structure_id']
        state['current_structure'] = state['structures'][sid]['data']
        s = Structure.from_dict(state['current_structure'])
        state['clashes'] = inspect_clashes(s)
        state['reaction_distances'] = [{'atoms': c['atoms'], 'distance_angstrom': distance(s, *c['atoms']), 'change': c['change']}
                                       for c in state['reaction']['bond_changes']]
        # Raw trajectories live in artifacts; keep prompts bounded without losing evidence.
        for j in state['jobs']:
            if 'parsed' in j:
                parsed = j['parsed']
                parsed['trajectory_frames'] = len(parsed.pop('trajectory', []))
        state['tools'] = dict(TOOL_PARAMETERS)
        state['tools'].update({name: {'set_distance': ['atoms', 'value', 'moving'], 'set_angle': ['atoms', 'value', 'moving'],
                                      'set_dihedral': ['atoms', 'value', 'moving'], 'translate_fragment': ['moving', 'vector'],
                                      'rotate_fragment': ['moving', 'axis', 'origin', 'degrees']}[name] for name in GEOMETRY_TOOLS})
        return state

    def execute(self, action):
        try:
            json.dumps(action, allow_nan=False)
        except (ValueError, TypeError):
            return {'ok': False, 'error': 'action must be JSON-compatible and contain only finite numbers'}
        with self.store.lock():
            self.state = self.store.read_json('state.json')
            self._recover()
            limits = Limits(**self.state['limits'])
            if self.state['usage']['actions'] >= limits.max_actions:
                return {'ok': False, 'error': 'action budget exhausted'}
            if self.state['paused'] and isinstance(action, dict) and action.get('actor', 'agent') == 'agent' and action.get('tool') != 'resume':
                return {'ok': False, 'error': 'paused: human takeover or explicit resume required'}
            self.state['usage']['actions'] += 1
            aid = '%06d' % self.state['usage']['actions']
            audit = {'id': aid, 'action': copy.deepcopy(action), 'actor': action.get('actor', 'agent') if isinstance(action, dict) else 'unknown',
                     'input_structure_id': self.state['current_structure_id'], 'status': 'running', 'started_at': time.time()}
            self.state['pending_action'] = {'id': aid, 'job_id': None, 'audit': copy.deepcopy(audit)}
            self.store.write_json('state.json', self.state)
            self.store.path('actions/' + aid).mkdir(exist_ok=False)
            self.store.write_json('actions/%s/action.json' % aid, audit)
            try:
                result = self._dispatch(action, aid, limits)
                answer = {'ok': True, 'action_id': aid, 'result': result}
                audit.update(status='completed', result=result)
                self.state['last_error'] = None
            except (ValueError, TypeError, KeyError, OSError) as exc:
                message = '%s: %s' % (type(exc).__name__, exc)
                answer = {'ok': False, 'action_id': aid, 'error': message}
                audit.update(status='failed', error=message)
                self.state['last_error'] = message
            except BaseException:
                audit.update(status='interrupted', error='execution interrupted; reservation retained')
                self.store.write_json('actions/%s/action.json' % aid, audit)
                self.store.write_json('state.json', self.state)
                raise
            audit.update(output_structure_id=self.state['current_structure_id'], finished_at=time.time(), usage=copy.deepcopy(self.state['usage']))
            self.store.write_json('actions/%s/action.json' % aid, audit)
            self.state['pending_action'] = None
            self.state['history'].append({'action_id': aid, 'tool': action.get('tool') if isinstance(action, dict) else None,
                                          'reason': action.get('reason') if isinstance(action, dict) else None,
                                          'ok': answer['ok'], 'error': answer.get('error')})
            self.store.write_json('state.json', self.state)
            return answer

    def _save_structure(self, structure, sid, source_job=None, make_current=True):
        self.store.write_text('structures/%s.xyz' % sid, structure.to_xyz())
        self.state['structures'][sid] = {'data': structure.to_dict(), 'source_job': source_job}
        if make_current:
            self.state['undo'].append(self.state['current_structure_id'])
            self.state['current_structure_id'] = sid
        self.state['report'] = None

    def _dispatch(self, action, aid, limits):
        if not isinstance(action, dict) or set(action) - {'tool', 'params', 'reason', 'actor'}:
            raise ValueError('action requires tool, params, reason and optional actor')
        if not isinstance(action.get('reason'), str) or not action['reason'].strip():
            raise ValueError('every action needs an explicit reason')
        if not isinstance(action.get('actor', 'agent'), str) or not action.get('actor', 'agent').strip():
            raise ValueError('actor must be a non-empty string')
        tool, params = action.get('tool'), action.get('params', {})
        if not isinstance(params, dict):
            raise ValueError('tool params must be an object')
        sid = self.state['current_structure_id']
        s = Structure.from_dict(self.state['structures'][sid]['data'])
        if tool in GEOMETRY_TOOLS:
            edited = edit_geometry(s, tool, params)
            clashes = inspect_clashes(edited)
            if clashes:
                raise ValueError('geometry edit introduces gross overlap: %s' % clashes)
            if edited.fingerprint == s.fingerprint:
                raise ValueError('geometry edit has no effect')
            self._save_structure(edited, 'structure_' + aid)
            return {'structure_id': self.state['current_structure_id'], 'clashes': clashes}
        if tool not in TOOL_PARAMETERS or set(params) - set(TOOL_PARAMETERS[tool]):
            raise ValueError('unknown tool or unexpected parameters')
        if tool == 'inspect_geometry':
            return {'structure': s.to_dict(), 'clashes': inspect_clashes(s)}
        if tool == 'use_structure':
            target = params['structure_id']
            if target not in self.state['structures']:
                raise ValueError('unknown registered structure')
            self.state['undo'].append(sid)
            self.state['current_structure_id'] = target
            self.state['report'] = None
            return {'structure_id': target}
        if tool == 'undo_geometry_edit':
            if not self.state['undo']:
                raise ValueError('no previous structure to restore')
            self.state['current_structure_id'] = self.state['undo'].pop()
            self.state['report'] = None
            return {'structure_id': self.state['current_structure_id']}
        if tool == 'save_hypothesis':
            hypothesis = params['hypothesis']
            if not isinstance(hypothesis, str) or not hypothesis.strip():
                raise ValueError('hypothesis must be non-empty text')
            self.state['hypotheses'].append({'text': hypothesis, 'structure_id': sid, 'action_id': aid})
            return self.state['hypotheses'][-1]
        if tool in ('pause', 'resume', 'request_human_review'):
            if tool == 'request_human_review' and (not isinstance(params.get('question'), str) or not params['question'].strip()):
                raise ValueError('human review requires a question')
            self.state['paused'] = tool != 'resume'
            return {'paused': self.state['paused'], 'question': params.get('question')}
        if tool == 'run_gaussian':
            return self._run_gaussian(params, aid, limits)
        if tool == 'verify_ts':
            evidence = self._verified_evidence()
            report = validate_ts(Reaction.from_dict(self.state['reaction']), evidence, params.get('ts_job_id'))
            self.state['report'] = report
            self.store.write_json('report.json', report)
            return report
        raise ValueError('unsupported tool')

    def _run_gaussian(self, params, aid, limits):
        settings = dict(self.state.get('job_defaults', {}))
        settings.update({k: v for k, v in params.items() if k in JobSpec.__dataclass_fields__})
        spec = JobSpec(**settings)
        timeout = finite(params.get('timeout_seconds', limits.job_timeout_seconds), 'timeout')
        if timeout <= 0 or timeout > limits.job_timeout_seconds or spec.cpus > limits.max_cpus or spec.memory_mb > limits.max_memory_mb:
            raise ValueError('requested job exceeds per-job CPU/memory/time limits')
        usage = self.state['usage']
        if usage['jobs'] >= limits.max_jobs or usage['wall_seconds'] + timeout > limits.max_wall_seconds or usage['core_seconds'] + timeout*spec.cpus > limits.max_core_seconds:
            raise ValueError('Gaussian job budget exhausted; full timeout allowance must fit')
        sid = params.get('structure_id', self.state['current_structure_id'])
        if sid not in self.state['structures']:
            raise ValueError('unknown input structure')
        entry = self.state['structures'][sid]
        s = Structure.from_dict(entry['data'])
        if inspect_clashes(s):
            raise ValueError('cannot submit geometry with gross atom overlaps')
        parent_id = params.get('parent_job', entry['source_job'])
        if spec.kind == 'freq' or spec.kind.startswith('irc_') or spec.kind.startswith('endpoint_'):
            parent = next((j for j in self.state['jobs'] if j['id'] == parent_id), None)
            allowed_parent = TS_KINDS if not spec.kind.startswith('endpoint_') else ('irc_' + spec.kind.split('_', 1)[1],)
            if not parent or parent['kind'] not in allowed_parent or parent['status'] != 'completed' or not same_geometry(entry['data'], parent.get('output_structure')):
                raise ValueError('frequency/IRC/endpoint requires the corresponding completed parent and its output geometry')
            if parent['level'] != spec.level:
                raise ValueError('child calculation must use parent calculation level')
        reaction = Reaction.from_dict(self.state['reaction'])
        input_text = prepare_input(reaction, s, spec)
        fingerprint = hashlib.sha256(input_text.encode()).hexdigest()
        if fingerprint in self.state['fingerprints']:
            raise ValueError('duplicate calculation rejected; change geometry or scientific settings')
        usage['jobs'] += 1
        jid = 'job_%04d' % usage['jobs']
        relative = 'jobs/' + jid
        directory = self.store.path(relative)
        sections = [s]
        if spec.kind in ('qst2', 'qst3'):
            sections = [reaction.reactant, reaction.product] + ([s] if spec.kind == 'qst3' else [])
        actual_input = sections[0] if spec.kind == 'qst2' else s
        actual_input_id = 'reactant' if spec.kind == 'qst2' else sid
        job = {'id': jid, 'kind': spec.kind, 'level': spec.level, 'spec': spec.to_dict(), 'constraints': spec.constraints,
               'parent_job': parent_id, 'input_structure': actual_input.to_dict(), 'input_structures': [x.to_dict() for x in sections],
               'input_structure_id': actual_input_id, 'output_structure': None,
               'output_structure_id': None, 'synthetic': bool(self.runner.synthetic), 'status': 'running', 'directory': relative}
        self.state['jobs'].append(job)
        self.state['fingerprints'].append(fingerprint)
        self.state['pending_action']['job_id'] = jid
        usage['wall_seconds'] += timeout
        usage['core_seconds'] += timeout * spec.cpus
        self.store.write_json('state.json', self.state)  # Reserve BEFORE spawning.
        started = time.monotonic()
        try:
            result = self.runner.run(directory, input_text, timeout)
        except (OSError, ValueError) as exc:
            job.update(status='failed', error=str(exc))
            elapsed = time.monotonic() - started
            usage['wall_seconds'] += elapsed - timeout
            usage['core_seconds'] += (elapsed - timeout)*spec.cpus
            self.store.write_json(relative + '/job.json', job)
            raise
        elapsed = result['wall_seconds']
        usage['wall_seconds'] += elapsed - timeout
        usage['core_seconds'] += (elapsed - timeout) * spec.cpus
        job.update(status=result['status'], wall_seconds=elapsed, allocated_core_seconds=elapsed*spec.cpus,
                   returncode=result['returncode'], log_path=relative + '/output.log')
        raw = self.store.path(job['log_path']).read_text(encoding='utf-8', errors='replace')
        parsed = parse_log(raw, atom_ids=s.atom_ids)
        job['parsed'] = parsed
        if job['status'] == 'completed' and not parsed['normal_termination']:
            job['status'] = 'failed'
        if parsed['final_structure']:
            out = Structure.from_dict(parsed['final_structure'])
            if out.symbols != s.symbols:
                job['status'] = 'failed'
                job['error'] = 'Gaussian output atom sequence changed'
            else:
                output_id = jid + '_output'
                job.update(output_structure=out.to_dict(), output_structure_id=output_id)
                promote = spec.kind in OPT_KINDS and job['status'] == 'completed' and parsed['optimization_converged']
                self._save_structure(out, output_id, jid, make_current=bool(promote))
        self.store.write_json(relative + '/parsed.json', parsed)
        self.store.write_json(relative + '/job.json', job)
        return {'job': copy.deepcopy(job)}

    def _verified_evidence(self):
        evidence = copy.deepcopy(self.state['jobs'])
        for job in evidence:
            path = job.get('log_path')
            if not path or 'parsed' not in job:
                job['status'] = 'failed'
                job['parsed'] = {}
                continue
            try:
                raw = self.store.path(path).read_text(encoding='utf-8', errors='replace')
                parsed = parse_log(raw, atom_ids=Reaction.from_dict(self.state['reaction']).reactant.atom_ids)
            except (OSError, ValueError):
                job['status'] = 'failed'
                job['parsed'] = {'diagnostics': ['MISSING_OR_UNREADABLE_EVIDENCE']}
                continue
            if parsed['log_sha256'] != job['parsed'].get('log_sha256'):
                job['status'] = 'failed'
                parsed['diagnostics'].append('LOG_CHANGED_SINCE_EXECUTION')
            job['parsed'] = parsed
        return evidence

    def evaluate(self, ts_job_id=None):
        """Read-only scientific evaluation with a saved report; no calculation or LLM action."""
        with self.store.lock():
            self.state = self.store.read_json('state.json')
            self._recover()
            report = validate_ts(Reaction.from_dict(self.state['reaction']), self._verified_evidence(), ts_job_id)
            report['usage'] = copy.deepcopy(self.state['usage'])
            report['last_error'] = self.state['last_error']
            report['hypotheses'] = copy.deepcopy(self.state['hypotheses'])
            self.state['report'] = report
            self.store.write_json('report.json', report)
            self.store.write_json('state.json', self.state)
            return report
