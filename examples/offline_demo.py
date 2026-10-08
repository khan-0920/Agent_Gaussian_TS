"""Feedback-dependent engineering demonstration. No LLM or Gaussian results claimed."""
import argparse
import json
import sys
from pathlib import Path
from gaussian_ts.agent import Agent
from gaussian_ts.cli import init_session
from gaussian_ts.runner import LocalRunner


class SyntheticRunner(LocalRunner):
    synthetic = True


def action(tool, params, reason):
    return {'tool': tool, 'params': params, 'reason': reason}


class FeedbackDemoPolicy:
    """Example policy reads actual process feedback; it is intentionally not an LLM."""
    def decide(self, observation):
        jobs = observation['jobs']
        if not observation['hypotheses']:
            return action('save_hypothesis', {'hypothesis': 'Illustrative collinear hydrogen exchange; synthetic only'}, 'record the demonstration hypothesis')
        if not jobs:
            return action('run_gaussian', {'kind': 'ts', 'timeout_seconds': 5}, 'test the initial candidate to obtain numerical feedback')
        if len(jobs) == 1:
            if observation['current_structure']['positions'][1][0] < 1.:
                return action('set_distance', {'atoms': [1, 2], 'value': 1.5, 'moving': [2]},
                              'first TS attempt reported Hessian failure; try a symmetric exchange geometry')
            return action('run_gaussian', {'kind': 'ts', 'timeout_seconds': 5}, 'retry after the feedback-driven geometry change')
        ts = jobs[1]
        plan = [('freq', ts), ('irc_forward', ts), ('irc_reverse', ts)]
        if len(jobs) < 5:
            kind, parent = plan[len(jobs)-2]
            return action('run_gaussian', {'kind': kind, 'structure_id': parent['output_structure_id'], 'parent_job': parent['id'], 'timeout_seconds': 5},
                          'collect separate synthetic validation evidence at the same TS geometry')
        if len(jobs) < 7:
            parent = jobs[3] if len(jobs) == 5 else jobs[4]
            direction = 'forward' if len(jobs) == 5 else 'reverse'
            return action('run_gaussian', {'kind': 'endpoint_' + direction, 'structure_id': parent['output_structure_id'],
                                           'parent_job': parent['id'], 'timeout_seconds': 5}, 'optimize the selected synthetic IRC endpoint')
        if observation['report'] is None:
            return action('verify_ts', {}, 'the validator must reject synthetic evidence even with a complete-looking trace')
        return action('pause', {}, 'offline demonstration finished; real Gaussian/LLM scientific acceptance is outstanding')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', default='runs/offline-demo')
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent
    session = init_session(directory/'hydrogen_exchange'/'reaction.json', args.run_dir)
    session.runner = SyntheticRunner([sys.executable, str(directory/'synthetic_gaussian.py')])

    def progress(action, result):
        print('%s: %s' % (action['tool'], 'ok' if result['ok'] else result.get('error')), file=sys.stderr, flush=True)

    result = Agent(session, FeedbackDemoPolicy(), progress).run(20)
    result['demo'] = 'SYNTHETIC ENGINEERING DEMO; NO REAL GAUSSIAN OR LLM CALCULATIONS'
    result['run_dir'] = str(session.store.root)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
