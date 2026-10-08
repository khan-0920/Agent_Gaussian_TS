import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class EndToEndTests(unittest.TestCase):
    def test_offline_feedback_demo_preserves_trace_and_never_validates(self):
        project = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temp:
            result = subprocess.run([sys.executable, '-m', 'examples.offline_demo', '--run-dir', str(Path(temp)/'run')],
                                    cwd=project, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            output = json.loads(result.stdout)
            self.assertEqual(output['stop_reason'], 'paused')
            self.assertEqual(output['report']['label'], 'INCONCLUSIVE')
            state = json.loads((Path(temp)/'run'/'state.json').read_text())
            self.assertEqual(state['jobs'][0]['status'], 'failed')
            self.assertEqual(state['jobs'][1]['status'], 'completed')
            self.assertEqual(len(state['jobs']), 7)
            self.assertTrue(all(job['synthetic'] for job in state['jobs']))
            self.assertGreaterEqual(len(list((Path(temp)/'run'/'actions').glob('*/action.json'))), 10)
