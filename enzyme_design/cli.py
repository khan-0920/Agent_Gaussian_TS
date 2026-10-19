"""Agent-friendly entry points for complete campaigns and individually chosen stages."""
import argparse
import json
import subprocess
import sys
from .campaign import Campaign
from .io import read_json, object_fields
from .spec import STAGES, validate_spec, validate_profile


def tools():
    def schema(properties=None, required=None):
        return {'type': 'object', 'properties': properties or {}, 'required': required or [], 'additionalProperties': False}
    return {
        'run_stage': {'description': 'Execute one dependency-checked, budgeted stage; inspect failures before retrying.',
                      'parameters': schema({'stage': {'type': 'string', 'enum': list(STAGES)}, 'retry': {'type': 'boolean'}}, ['stage'])},
        'run_pipeline': {'description': 'Execute all remaining stages through candidate evaluation and delivery.',
                         'parameters': schema({'retry': {'type': 'boolean'}})},
        'context': {'description': 'Read current stage states, budget, failure logs and next action.', 'parameters': schema()},
        'report': {'description': 'Verify artifact hashes and return candidate evaluation, never an enzyme-activity claim.', 'parameters': schema()},
        'pause': {'description': 'Pause at a stage boundary; a currently running worker holds the campaign lock.', 'parameters': schema()},
        'resume': {'description': 'Resume a paused campaign at a stage boundary.', 'parameters': schema()},
        'doctor': {'description': 'Check configured runtime paths, without submitting model inference.', 'parameters': schema()},
    }


def execute(campaign, action):
    object_fields(action, ('tool', 'params', 'reason'), label='enzyme action')
    if action['tool'] not in tools():
        raise ValueError('unknown enzyme tool')
    parameters = tools()[action['tool']]['parameters']
    object_fields(action['params'], parameters['required'], set(parameters['properties'])-set(parameters['required']), 'tool params')
    if not isinstance(action['reason'], str) or not action['reason'].strip():
        raise ValueError('action reason must be non-empty')
    if 'retry' in action['params'] and not isinstance(action['params']['retry'], bool):
        raise ValueError('retry must be boolean')
    if action['tool'] == 'run_stage':
        return campaign.run_stage(reason=action['reason'], **action['params'])
    if action['tool'] == 'run_pipeline':
        return campaign.run(reason=action['reason'], **action['params'])
    if action['tool'] in ('pause', 'resume'):
        return campaign.set_paused(action['tool'] == 'pause', action['reason'])
    return getattr(campaign, action['tool'])()


def main(argv=None):
    parser = argparse.ArgumentParser(description='Auditable theozyme → de novo enzyme generation and evaluation')
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('validate', 'init'):
        command = sub.add_parser(name)
        command.add_argument('--config', required=True)
        command.add_argument('--profile', required=True)
        if name == 'init':
            command.add_argument('--run-dir', required=True)
    catalog = sub.add_parser('tools')
    catalog.add_argument('--format', choices=('catalog', 'functions'), default='catalog')
    for name in ('doctor', 'context', 'report', 'pause', 'resume', 'run', 'stage', 'apply'):
        command = sub.add_parser(name)
        command.add_argument('--run-dir', required=True)
        if name in ('pause', 'resume', 'run', 'stage'):
            command.add_argument('--reason', required=True)
        if name in ('run', 'stage'):
            command.add_argument('--retry', action='store_true')
        if name == 'stage':
            command.add_argument('stage', choices=STAGES)
        if name == 'apply':
            command.add_argument('--action', required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == 'tools':
            result = tools()
            if args.format == 'functions':
                result = [{'type': 'function', 'function': {'name': name, **definition}} for name, definition in result.items()]
        elif args.command == 'validate':
            spec, profile = validate_spec(read_json(args.config)), validate_profile(read_json(args.profile))
            result = {'valid': True, 'name': spec['name'], 'stages': list(STAGES),
                      'expected_rfd3_foundry': profile['tools']['generate']['expected_foundry'],
                      'scope': 'Schema validation only; init snapshots inputs and prepare validates structures/chemistry.'}
        elif args.command == 'init':
            result = Campaign.create(args.config, args.profile, args.run_dir).context()
        else:
            campaign = Campaign(args.run_dir)
            if args.command == 'apply':
                result = execute(campaign, read_json(args.action))
            elif args.command == 'stage':
                result = campaign.run_stage(args.stage, args.reason, retry=args.retry)
            elif args.command == 'run':
                result = campaign.run(args.reason, retry=args.retry)
            elif args.command in ('pause', 'resume'):
                result = campaign.set_paused(args.command == 'pause', args.reason)
            else:
                result = getattr(campaign, args.command)()
        print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
        return 0
    except (ValueError, TypeError, KeyError, OSError, subprocess.SubprocessError) as error:
        print(json.dumps({'ok': False, 'error': type(error).__name__ + ': ' + str(error)}, ensure_ascii=False), file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print('Interrupted; worker process group stopped and attempt retained.', file=sys.stderr)
        return 130
