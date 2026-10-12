"""Regression against a real local G16 C.02 ammonia HF/STO-3G run.
Fixture is unmodified output, not a validated transition-state benchmark.
"""
from pathlib import Path
import re
import unittest
from gaussian_ts.parser import parse_log
from gaussian_ts.validation import usable


class RealGaussianParserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = (Path(__file__).parent / 'fixtures/g16_c02_ammonia_opt_freq.log').read_text()

    def test_internal_steps_are_separate_and_complete(self):
        p = parse_log(self.text, [1, 2, 3, 4])
        self.assertTrue(p['normal_termination'])
        self.assertTrue(p['multi_job'])
        self.assertEqual(len(p['steps']), 2)
        opt, freq = p['steps']
        self.assertTrue(opt['optimization_converged'])
        self.assertEqual(len(opt['convergence']), 4)
        self.assertEqual(opt['frequencies_cm1'], [])
        self.assertEqual(len(freq['frequencies_cm1']), 6)
        self.assertTrue(all(v > 0 for v in freq['frequencies_cm1']))
        self.assertTrue(all(len(m['displacements']) == 4 for m in freq['modes']))
        self.assertFalse(opt['parse_warnings'])
        self.assertFalse(freq['parse_warnings'])
        self.assertFalse(p['optimization_converged'])
        self.assertIsNone(p['final_structure'])
        self.assertEqual(p['frequencies_cm1'], [])
        self.assertFalse(usable({'status': 'completed', 'synthetic': False, 'parsed': p}))

    def test_truncated_or_failed_second_step_is_not_success(self):
        for ending in ['', 'Error termination via Lnk1e']:
            text = self.text.rsplit('Normal termination of Gaussian', 1)[0] + ending
            p = parse_log(text)
            self.assertFalse(p['normal_termination'])
            self.assertTrue(p['steps'][0]['normal_termination'])
            self.assertFalse(p['steps'][1]['normal_termination'])

    def test_rms_whitespace_and_failed_last_table(self):
        first = re.split(r'Link1:\s*Proceeding to internal job step number\s+2\.', self.text)[0]
        self.assertTrue(parse_log(first)['optimization_converged'])
        failed = first + '\n Item Value Threshold Converged?\n RMS\t Force 0.9 0.000300 NO\n'
        p = parse_log(failed)
        self.assertFalse(p['optimization_converged'])
        self.assertEqual(len(p['convergence']), 1)
