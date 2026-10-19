"""Stage entry point, run in the selected model's isolated Python environment."""
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import traceback
from pathlib import Path
from .io import read_json, write_json, sha256, within


class Work:
    def __init__(self, request):
        self.request = request
        self.root = Path(request['root']).resolve()
        self.stage = request['stage']
        self.out = within(self.root, request['attempt_dir'])
        self.spec, self.profile = request['spec'], request['profile']
        self.tool = self.profile['tools'][self.stage]
        self.previous = request['previous']
        self.dependencies = request.get('dependencies', {})
        self.commands = []

    def path(self, relative):
        return within(self.root, relative)

    def relative(self, path):
        relative = str(Path(path).resolve().relative_to(self.root))
        within(self.root, relative)
        return relative

    def run(self, command, cwd=None):
        if not all(isinstance(part, str) for part in command):
            raise ValueError('command arguments must be strings')
        self.commands.append({'argv': command, 'cwd': str(cwd or self.out)})
        write_json(self.out / 'commands.json', self.commands)
        subprocess.run(command, cwd=cwd or self.out, check=True)

    def manifest(self, candidates, rejected=None, summary=None):
        artifacts = sorted(self.relative(p) for p in self.out.rglob('*') if p.is_file()
                           and p.name not in ('result.json', 'worker.log'))
        result = {'schema_version': 1, 'stage': self.stage, 'candidates': candidates,
                  'rejected': rejected or [], 'summary': summary or {}, 'artifacts': artifacts,
                  'interpretation': 'Computational design/screening only; no experimental activity claim.'}
        write_json(self.out / 'result.json', result)
        return result

    def provenance(self):
        versions = {}
        for name in ('rc-foundry', 'torch', 'atomworks', 'biotite', 'rdkit', 'pyrosetta', 'numpy'):
            try:
                versions[name] = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:
                pass
        weights = {}
        for key in ('checkpoint', 'packer_checkpoint'):
            if key in self.tool:
                path = Path(self.tool[key])
                weights[key] = {'path': str(path), 'sha256': sha256(path), 'bytes': path.stat().st_size}
        source = {'worker_sha256': sha256(__file__)}
        for name in ('adapters.py', 'structure.py', 'spec.py'):
            file = Path(__file__).with_name(name)
            if file.exists():
                source[name] = sha256(file)
        if 'repo' in self.tool:
            command = subprocess.run(['git', '-C', self.tool['repo'], 'rev-parse', 'HEAD'],
                                     text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if command.returncode == 0:
                source['upstream_commit'] = command.stdout.strip()
        record = {'python': sys.executable, 'python_version': platform.python_version(), 'versions': versions,
                  'weights': weights, 'sources': source, 'gpu': os.environ.get('CUDA_VISIBLE_DEVICES'),
                  'generation_mode': self.spec['generation']['inference_mode'],
                  'seed': self.spec['generation']['seed'], 'runtime': self.tool}
        write_json(self.out / 'provenance.json', record)
        if self.stage == 'generate' and versions.get('rc-foundry') != self.tool['expected_foundry']:
            raise ValueError('RFD3 runtime version mismatch: expected ' + self.tool['expected_foundry'])
        return record


def main():
    work = Work(read_json(sys.argv[1]))
    try:
        work.provenance()
        from . import adapters
        getattr(adapters, work.stage)(work)
        return 0
    except Exception as exc:
        traceback.print_exc()
        write_json(work.out / 'error.json', {'error': type(exc).__name__, 'message': str(exc)})
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
