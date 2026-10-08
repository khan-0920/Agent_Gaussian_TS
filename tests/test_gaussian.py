import sys
import tempfile
import unittest
from pathlib import Path

from gaussian_ts.models import Structure, Reaction
from gaussian_ts.gaussian import JobSpec, prepare_input
from gaussian_ts.parser import parse_log
from gaussian_ts.runner import LocalRunner
from tests.logs import log


class GaussianTests(unittest.TestCase):
    def setUp(self):
        self.s = Structure.from_xyz('3\nhydrogen exchange\nH 0 0 0\nH 1.5 0 0\nH 3 0 0\n', atom_ids=[10, 20, 30])
        self.r = Reaction(self.s, self.s, 0, 2, ({'atoms': [10, 20], 'change': -1}, {'atoms': [20, 30], 'change': 1}))

    def test_ts_has_calcfc_no_constraints_or_noeigentest(self):
        inp = prepare_input(self.r, self.s, JobSpec('ts'))
        self.assertIn('Opt=(TS,CalcFC,Tight,MaxCycles=100)', inp)
        self.assertNotIn('NoEigenTest', inp)
        self.assertIn('0 2\nH', inp)

    def test_modredundant_maps_stable_ids_to_gaussian_indices(self):
        inp = prepare_input(self.r, self.s, JobSpec('frozen_opt', constraints=({'type': 'B', 'atoms': [10, 30], 'operation': 'F'},)))
        self.assertIn('\nB 1 3 F\n', inp)

    def test_scan_and_qst_molecule_sections(self):
        inp = prepare_input(self.r, self.s, JobSpec('scan', constraints=({'type': 'B', 'atoms': [10, 20], 'operation': 'S', 'steps': 10, 'step_size': -.1},)))
        self.assertIn('B 1 2 S 10 -0.1', inp)
        for kind, count in [('qst2', 2), ('qst3', 3)]:
            inp = prepare_input(self.r, self.s, JobSpec(kind))
            self.assertEqual(inp.count('\n0 2\n'), count)
            self.assertNotIn('Opt=(TS', inp)

    def test_injection_and_constraints_on_ts_rejected(self):
        for kwargs in [{'method': 'HF\n--Link1--'}, {'memory_mb': 0}, {'cpus': True}, {'kind': 'shell'}, {'kind': 'ts', 'constraints': ({'type': 'B', 'atoms': [10, 20], 'operation': 'F'},)}]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                spec = JobSpec(**dict({'kind': 'ts'}, **kwargs))
                prepare_input(self.r, self.s, spec)

    def test_parser_energy_convergence_and_mapping(self):
        result = parse_log(log(), atom_ids=self.s.atom_ids)
        self.assertTrue(result['normal_termination'])
        self.assertTrue(result['optimization_converged'])
        self.assertAlmostEqual(result['energy_hartree'], -1.5)
        self.assertEqual(result['final_structure']['atom_ids'], (10, 20, 30))

    def test_parser_frequencies_modes_without_low_frequency_duplicates(self):
        result = parse_log(' Low frequencies --- -0.1 0.0 0.1\n' + log('freq'), atom_ids=self.s.atom_ids)
        self.assertEqual(result['frequencies_cm1'], [-500., 100., 200., 300.])
        self.assertEqual(result['modes'][0]['displacements'][1], [1., 0., 0.])

    def test_truncated_failed_and_composite_logs_never_converge(self):
        text = log().replace('RMS Displacement 0.000001 0.001200 YES', 'RMS Displacement 0.9 0.001200 NO')
        self.assertFalse(parse_log(text)['optimization_converged'])
        self.assertFalse(parse_log(log(normal=False))['normal_termination'])
        for text in [log() + '--Link1--\n' + log('freq'), log() + ' Link1: Proceeding to internal job step number 2.\n unfinished']:
            self.assertTrue(parse_log(text)['multi_job'])

    def test_scf_failure_diagnosis_and_irc_points(self):
        result = parse_log(' Convergence failure -- run terminated.\n Error termination via Lnk1e')
        self.assertIn('SCF_FAILURE', result['diagnostics'])
        result = parse_log(log('irc'), atom_ids=self.s.atom_ids)
        self.assertEqual(result['irc_points'], [1, 2])
        self.assertTrue(result['irc_complete'])

    def test_runner_real_process_stdin_and_timeout(self):
        with tempfile.TemporaryDirectory() as temp:
            runner = LocalRunner([sys.executable, '-c', 'import sys; print(sys.stdin.read())'])
            job = Path(temp) / 'job'
            handle = runner.start(job, 'input test\n')
            result = runner.wait(handle, 2)
            self.assertEqual(result['status'], 'completed')
            self.assertIn('input test', (job / 'output.log').read_text())
            runner = LocalRunner([sys.executable, '-c', 'import time; time.sleep(5)'])
            handle = runner.start(Path(temp) / 'timeout', 'x')
            self.assertEqual(runner.wait(handle, .03)['status'], 'timed_out')
            self.assertIsNotNone(handle.process.poll())

    def test_runner_cannot_overwrite_existing_directory(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(FileExistsError):
                LocalRunner([sys.executable]).start(Path(temp), 'x')

    def test_second_job_start_without_termination_invalidates_previous_success(self):
        text = log() + '\n Entering Gaussian System, Link 0=g16\n' + log('freq', normal=False).split('Error termination')[0]
        parsed = parse_log(text)
        self.assertTrue(parsed['multi_job'])
        self.assertFalse(parsed['normal_termination'])
        self.assertFalse(parsed['optimization_converged'])

    def test_nan_in_frequency_table_is_parse_warning(self):
        parsed = parse_log(log('freq').replace('100.0', 'nan'))
        self.assertTrue(parsed['parse_warnings'])


if __name__ == '__main__':
    unittest.main()
