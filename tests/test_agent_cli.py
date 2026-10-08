import json
import os
import subprocess
import sys
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

from gaussian_ts.agent import Agent, HTTPPolicy, ReplayPolicy
from gaussian_ts.engine import Session
from tests.test_engine import reaction, action


class FeedbackPolicy:
    def decide(self, observation):
        if observation['usage']['actions'] == 0:
            return action('set_distance', {'atoms': [1, 2], 'value': .1})
        if observation['last_error']:
            return action('set_distance', {'atoms': [1, 2], 'value': 1.2})
        return action('request_human_review', {'question': 'geometry is ready; confirm calculation settings'})


class AgentTests(unittest.TestCase):
    def test_feedback_drives_correction_and_takeover(self):
        with tempfile.TemporaryDirectory() as temp:
            session = Session.create(Path(temp)/'run', reaction())
            result = Agent(session, FeedbackPolicy()).run(max_steps=5)
            self.assertEqual(result['stop_reason'], 'human_review')
            self.assertEqual(session.observe()['usage']['actions'], 3)
            self.assertAlmostEqual(session.observe()['current_structure']['positions'][1][0], 1.2)
            self.assertEqual(session.observe()['history'][0]['ok'], False)

    def test_replay_pause_and_exhaustion_are_honest(self):
        with tempfile.TemporaryDirectory() as temp:
            s = Session.create(Path(temp)/'run', reaction())
            result = Agent(s, ReplayPolicy([action('verify_ts')])).run(5)
            self.assertEqual(result['stop_reason'], 'policy_exhausted')
            self.assertFalse(result['report']['validated'])

    def test_http_policy_real_server_json_and_prompt(self):
        captured = []

        class Handler(BaseHTTPRequestHandler):
            def do_POST(self):
                captured.append(json.loads(self.rfile.read(int(self.headers['Content-Length']))))
                self.send_response(200)
                self.send_header('Content-Type', 'application/json')
                self.end_headers()
                self.wfile.write(json.dumps({'choices': [{'message': {'content': json.dumps(action('inspect_geometry'))}}]}).encode())

            def log_message(self, *args):
                pass

        server = HTTPServer(('127.0.0.1', 0), Handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            policy = HTTPPolicy('http://127.0.0.1:%d/chat/completions' % server.server_port, 'test-model', api_key_env=None)
            result = policy.decide({'tools': {}, 'usage': {'actions': 0}})
            self.assertEqual(result['tool'], 'inspect_geometry')
            self.assertIn('VALIDATED_TS', captured[0]['messages'][0]['content'])
            self.assertEqual(captured[0]['model'], 'test-model')
        finally:
            server.shutdown()
            server.server_close()
            thread.join()

    def test_http_policy_rejects_insecure_remote_url_and_missing_key(self):
        with self.assertRaises(ValueError):
            HTTPPolicy('http://example.com/chat/completions', 'model')
        with self.assertRaises(ValueError):
            HTTPPolicy('https://user:secret@example.com/chat/completions', 'model')
        with self.assertRaises(ValueError):
            HTTPPolicy('https://example.com/chat/completions', 'model', api_key_env='GAUSSIAN_TS_TEST_MISSING_KEY')

    def test_cli_init_apply_status_and_report(self):
        project = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'reactant.xyz').write_text(reaction().reactant.to_xyz())
            (root/'product.xyz').write_text(reaction().product.to_xyz())
            config = {'reactant': 'reactant.xyz', 'product': 'product.xyz', 'charge': 0, 'multiplicity': 2,
                      'bond_changes': list(reaction().bond_changes), 'gaussian_command': ['g16']}
            (root/'reaction.json').write_text(json.dumps(config))
            (root/'action.json').write_text(json.dumps(action('inspect_geometry')))
            commands = [['init', '--config', str(root/'reaction.json'), '--run-dir', str(root/'run')],
                        ['apply', '--run-dir', str(root/'run'), '--action', str(root/'action.json')],
                        ['status', '--run-dir', str(root/'run')], ['report', '--run-dir', str(root/'run')]]
            for cmd in commands:
                result = subprocess.run([sys.executable, '-m', 'gaussian_ts'] + cmd, cwd=project, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIsInstance(json.loads(result.stdout), dict)


if __name__ == '__main__':
    unittest.main()
