"""Command line entry points for manual actions and provider-driven research."""
import argparse
import json
import shutil
import sys
from pathlib import Path
from .agent import Agent, HTTPPolicy, ReplayPolicy
from .engine import Session
from .models import Structure, Reaction, Limits
from .toolkit import define_problem, export_context, tool_catalog, function_tools


def read_json(path):
    with Path(path).open(encoding='utf-8') as source:
        return json.load(source)


def init_session(config_path, run_dir):
    path = Path(config_path).resolve()
    config = read_json(path)
    required = {'reactant', 'product', 'charge', 'multiplicity', 'bond_changes'}
    allowed = required | {'atom_ids', 'limits', 'gaussian_command', 'job_defaults', 'problem'}
    if not isinstance(config, dict) or not required.issubset(config) or set(config) - allowed:
        raise ValueError('reaction config has missing or unknown fields')
    mapping = config.get('atom_ids')
    r = Structure.from_xyz((path.parent / config['reactant']).read_text(encoding='utf-8'), mapping)
    p = Structure.from_xyz((path.parent / config['product']).read_text(encoding='utf-8'), mapping)
    reaction = Reaction(r, p, config['charge'], config['multiplicity'], tuple(config['bond_changes']))
    problem = define_problem(config['problem'], config) if 'problem' in config else None
    return Session.create(run_dir, reaction, Limits(**config.get('limits', {})),
                          config.get('gaussian_command', ['g16']), config.get('job_defaults', {}), problem)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Gaussian transition-state problem and tool kit for Codex/GPT')
    sub = parser.add_subparsers(dest='command', required=True)
    init = sub.add_parser('init', help='initialize a human-guided research session')
    init.add_argument('--config', required=True)
    init.add_argument('--run-dir', required=True)
    apply = sub.add_parser('apply', help='apply one audited human action')
    apply.add_argument('--run-dir', required=True)
    apply.add_argument('--action', required=True)
    apply.add_argument('--actor', default='human')
    for name in ('status', 'context', 'report', 'pause', 'resume'):
        cmd = sub.add_parser(name)
        cmd.add_argument('--run-dir', required=True)
    run = sub.add_parser('run', help='decision loop using an LLM or explicit action replay')
    run.add_argument('--run-dir', required=True)
    policy = run.add_mutually_exclusive_group(required=True)
    policy.add_argument('--llm-config')
    policy.add_argument('--actions')
    run.add_argument('--steps', type=int, default=20)
    doctor = sub.add_parser('doctor', help='check Python and Gaussian availability')
    doctor.add_argument('--gaussian-executable', default='g16')
    define = sub.add_parser('define', help='frame a research problem and expose missing calculation inputs')
    define.add_argument('--problem', required=True)
    define.add_argument('--config')
    catalog = sub.add_parser('tools', help='export typed tools for Codex/GPT')
    catalog.add_argument('--format', choices=['catalog', 'functions'], default='catalog')
    args = parser.parse_args(argv)
    try:
        if args.command == 'define':
            result = define_problem(read_json(args.problem), read_json(args.config) if args.config else None)
        elif args.command == 'tools':
            result = function_tools() if args.format == 'functions' else tool_catalog()
        elif args.command == 'doctor':
            executable = shutil.which(args.gaussian_executable)
            result = {'python': sys.version.split()[0], 'gaussian_executable': executable, 'gaussian_available': bool(executable),
                      'note': 'Availability only; license/environment and a real smoke calculation must be checked on the compute host.'}
        elif args.command == 'init':
            session = init_session(args.config, args.run_dir)
            result = {'run_dir': str(session.store.root), 'status': 'initialized', 'mode': 'human_guided'}
        else:
            session = Session.open(args.run_dir)
            if args.command == 'apply':
                action = read_json(args.action)
                if not isinstance(action, dict):
                    raise ValueError('action must be a JSON object')
                result = session.execute(dict(action, actor=args.actor))
            elif args.command in ('pause', 'resume'):
                result = session.execute({'tool': args.command, 'params': {}, 'reason': 'user checkpoint control', 'actor': 'human'})
            elif args.command == 'status':
                result = session.observe()
            elif args.command == 'context':
                result = export_context(session)
            elif args.command == 'report':
                result = session.evaluate()
            else:
                decision_policy = ReplayPolicy(read_json(args.actions)) if args.actions else HTTPPolicy(**read_json(args.llm_config))
                def progress(action, outcome):
                    print(json.dumps({'tool': action.get('tool') if isinstance(action, dict) else None,
                                      'ok': outcome['ok'], 'action_id': outcome.get('action_id'), 'error': outcome.get('error')}, ensure_ascii=False), file=sys.stderr, flush=True)
                result = Agent(session, decision_policy, progress).run(args.steps)
        print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
        return 1 if isinstance(result, dict) and result.get('ok') is False else 0
    except (ValueError, TypeError, KeyError, OSError) as exc:
        print(json.dumps({'ok': False, 'error': '%s: %s' % (type(exc).__name__, exc)}, ensure_ascii=False), file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print('Interrupted; checkpoint and reserved budget preserved.', file=sys.stderr)
        return 130
