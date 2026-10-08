import json
import sys
import tempfile
import unittest
from pathlib import Path
from gaussian_ts.models import Structure, Reaction, Limits
from gaussian_ts.runner import LocalRunner
from gaussian_ts.engine import Session
from gaussian_ts.store import Store


def reaction():
    r = Structure.from_xyz('3\nr\nH 0 0 0\nH 0.75 0 0\nH 3 0 0\n')
    p = Structure.from_xyz('3\np\nH 0 0 0\nH 2.25 0 0\nH 3 0 0\n')
    return Reaction(r, p, 0, 2, ({'atoms': [1, 2], 'change': -1}, {'atoms': [2, 3], 'change': 1}))


def action(tool, params=None, actor='human'):
    return {'tool': tool, 'params': params or {}, 'reason': 'test controlled operation', 'actor': actor}


class EngineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name) / 'run'
        self.session = Session.create(self.root, reaction(), Limits(), ['missing-gaussian-test'])

    def tearDown(self):
        self.temp.cleanup()

    def test_geometry_undo_checkpoint_and_audit(self):
        result = self.session.execute(action('set_distance', {'atoms': [1, 2], 'value': 1.2, 'moving': [2]}))
        self.assertTrue(result['ok'])
        s = Session.open(self.root)
        self.assertAlmostEqual(s.observe()['current_structure']['positions'][1][0], 1.2)
        self.assertTrue(s.execute(action('undo_geometry_edit'))['ok'])
        self.assertAlmostEqual(s.observe()['current_structure']['positions'][1][0], .75)
        audit = json.loads((self.root / 'actions' / '000001' / 'action.json').read_text())
        self.assertEqual(audit['actor'], 'human')
        self.assertNotEqual(audit['input_structure_id'], audit['output_structure_id'])

    def test_failed_edit_is_audited_and_does_not_change_structure(self):
        result = self.session.execute(action('set_distance', {'atoms': [1, 2], 'value': .1}))
        self.assertFalse(result['ok'])
        self.assertEqual(self.session.observe()['current_structure_id'], 'reactant')
        self.assertEqual(self.session.observe()['usage']['actions'], 1)

    def test_action_limits_and_duplicate_jobs(self):
        self.session.state['limits']['max_actions'] = 1
        self.session.store.write_json('state.json', self.session.state)
        self.assertTrue(self.session.execute(action('inspect_geometry'))['ok'])
        self.assertFalse(self.session.execute(action('inspect_geometry'))['ok'])

    def test_job_resource_limits_prevent_process_start(self):
        for params in [{'kind': 'ts', 'cpus': 999}, {'kind': 'ts', 'memory_mb': 999999}, {'kind': 'ts', 'timeout_seconds': 99999}]:
            result = self.session.execute(action('run_gaussian', params))
            self.assertFalse(result['ok'])
        self.assertEqual(self.session.observe()['usage']['jobs'], 0)

    def test_missing_executable_is_audited_with_input(self):
        result = self.session.execute(action('run_gaussian', {'kind': 'ts'}))
        self.assertFalse(result['ok'])
        self.assertTrue((self.root / 'jobs' / 'job_0001' / 'input.gjf').exists())
        self.assertEqual(len(self.session.observe()['jobs']), 1)
        self.assertFalse(self.session.execute(action('run_gaussian', {'kind': 'ts'}))['ok'])

    def test_budget_reservation_and_recovery_are_conservative(self):
        runner = LocalRunner([sys.executable, '-c', 'import time; time.sleep(5)'])
        self.session.runner = runner
        result = self.session.execute(action('run_gaussian', {'kind': 'ts', 'timeout_seconds': .03}))
        self.assertTrue(result['ok'])
        self.assertEqual(result['result']['job']['status'], 'timed_out')
        self.assertGreater(self.session.observe()['usage']['core_seconds'], 0)
        state = self.session.store.read_json('state.json')
        state['pending_action'] = {'id': '000099', 'job_id': None}
        self.session.store.write_json('state.json', state)
        (self.root / 'actions' / '000099').mkdir()
        self.session.store.write_json('actions/000099/action.json', {'status': 'running'})
        s = Session.open(self.root)
        self.assertIsNone(s.observe()['pending_action'])
        self.assertEqual(s.store.read_json('actions/000099/action.json')['status'], 'interrupted')

    def test_pause_resume_and_human_override_record(self):
        self.assertTrue(self.session.execute(action('pause'))['ok'])
        self.assertFalse(self.session.execute(action('inspect_geometry', actor='agent'))['ok'])
        self.assertTrue(self.session.execute(action('use_structure', {'structure_id': 'product'}, actor='chemist'))['ok'])
        self.assertTrue(self.session.execute(action('resume'))['ok'])
        self.assertEqual(self.session.observe()['current_structure_id'], 'product')

    def test_llm_cannot_set_validation_or_execute_shell(self):
        self.assertFalse(self.session.execute(action('shell', {'command': 'rm -rf x'}))['ok'])
        self.assertFalse(self.session.execute(action('verify_ts', {'label': 'VALIDATED_TS'}))['ok'])
        result = self.session.execute(action('verify_ts'))
        self.assertEqual(result['result']['label'], 'INCONCLUSIVE')

    def test_store_rejects_traversal_absolute_paths_and_symlinks(self):
        for path in ['../escape.json', str(Path(self.temp.name) / 'escape.json')]:
            with self.assertRaises(ValueError):
                self.session.store.write_json(path, {})
        outside = Path(self.temp.name) / 'outside'
        outside.mkdir()
        (self.root / 'escape').symlink_to(outside, target_is_directory=True)
        with self.assertRaises(ValueError):
            self.session.store.write_json('escape/file.json', {})

    def test_stale_session_instance_reloads_checkpoint(self):
        other = Session.open(self.root)
        self.session.execute(action('use_structure', {'structure_id': 'product'}))
        other.execute(action('save_hypothesis', {'hypothesis': 'alternate geometry'}))
        self.assertEqual(other.observe()['current_structure_id'], 'product')

    def test_non_json_action_does_not_corrupt_next_audit_sequence(self):
        result = self.session.execute(action('set_distance', {'atoms': [1, 2], 'value': float('nan')}))
        self.assertFalse(result['ok'])
        session = Session.open(self.root)
        self.assertTrue(session.execute(action('inspect_geometry'))['ok'])

    def test_orphan_action_directory_is_recovered_before_next_action(self):
        self.session.store.path('actions/000001').mkdir()
        self.session.store.write_json('actions/000001/action.json', {'status': 'running', 'id': '000001'})
        session = Session.open(self.root)
        self.assertTrue(session.execute(action('inspect_geometry'))['ok'])
        self.assertEqual(session.store.read_json('actions/000001/action.json')['status'], 'interrupted')

    def test_qst2_dedup_uses_actual_reactant_product_not_current_guess(self):
        self.session.execute(action('run_gaussian', {'kind': 'qst2'}))
        self.session.execute(action('set_distance', {'atoms': [1, 2], 'value': 1.2}))
        result = self.session.execute(action('run_gaussian', {'kind': 'qst2'}))
        self.assertFalse(result['ok'])
        self.assertEqual(self.session.observe()['usage']['jobs'], 1)
        job = self.session.state['jobs'][0]
        self.assertEqual(job['input_structure_id'], 'reactant')
        self.assertEqual(len(job['input_structures']), 2)

    def test_pending_without_directory_recovers_full_action_audit(self):
        state = self.session.store.read_json('state.json')
        original = action('save_hypothesis', {'hypothesis': 'interrupted intent'}, actor='chemist')
        state['usage']['actions'] = 1
        state['pending_action'] = {'id': '000001', 'job_id': None,
                                   'audit': {'id': '000001', 'action': original, 'actor': 'chemist',
                                             'input_structure_id': 'reactant', 'status': 'running'}}
        self.session.store.write_json('state.json', state)
        session = Session.open(self.root)
        recovered = session.store.read_json('actions/000001/action.json')
        self.assertEqual(recovered['status'], 'interrupted')
        self.assertEqual(recovered['action'], original)
        self.assertTrue(session.execute(action('inspect_geometry'))['ok'])

    def test_removed_log_is_reported_as_missing_evidence(self):
        from tests.logs import log
        self.session.runner = LocalRunner([sys.executable, '-c', 'print(%r)' % log()])
        result = self.session.execute(action('run_gaussian', {'kind': 'ts'}))
        self.assertTrue(result['ok'])
        self.session.store.path('jobs/job_0001/output.log').unlink()
        report = self.session.evaluate()
        self.assertFalse(report['validated'])


if __name__ == '__main__':
    unittest.main()
