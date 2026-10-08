import copy
import unittest
from gaussian_ts.models import Structure, Reaction
from gaussian_ts.parser import parse_log
from gaussian_ts.validation import validate_ts
from tests.logs import log


def evidence(job_id, kind, xs, parent=None):
    parsed = parse_log(log('irc' if kind.startswith('irc_') else ('freq' if kind == 'freq' else 'opt'), xs))
    s = parsed['final_structure']
    return {'id': job_id, 'kind': kind, 'level': 'hf/sto-3g|SMD:|Integral:UltraFine', 'synthetic': False,
            'status': 'completed', 'constraints': [], 'parent_job': parent, 'input_structure': s,
            'output_structure': s, 'parsed': parsed}


class ValidationTests(unittest.TestCase):
    def setUp(self):
        r = Structure.from_xyz('3\nr\nH 0 0 0\nH 0.75 0 0\nH 3 0 0\n')
        p = Structure.from_xyz('3\np\nH 0 0 0\nH 2.25 0 0\nH 3 0 0\n')
        self.reaction = Reaction(r, p, 0, 2, ({'atoms': [1, 2], 'change': -1}, {'atoms': [2, 3], 'change': 1}))
        self.ts = evidence('ts1', 'ts', (0, 1.5, 3))
        self.freq = evidence('freq1', 'freq', (0, 1.5, 3), 'ts1')
        self.forward = evidence('f1', 'irc_forward', (0, 2.25, 3), 'ts1')
        self.reverse = evidence('r1', 'irc_reverse', (0, .75, 3), 'ts1')
        self.forward['input_structure'] = self.ts['output_structure']
        self.reverse['input_structure'] = self.ts['output_structure']
        self.end_f = evidence('ef1', 'endpoint_forward', (0, 2.25, 3), 'f1')
        self.end_r = evidence('er1', 'endpoint_reverse', (0, .75, 3), 'r1')
        self.jobs = [self.ts, self.freq, self.forward, self.reverse, self.end_f, self.end_r]

    def test_normal_termination_alone_is_not_a_ts(self):
        self.ts['parsed']['optimization_converged'] = False
        self.assertNotEqual(validate_ts(self.reaction, self.jobs)['label'], 'VALIDATED_TS')

    def test_one_frequency_without_irc_is_only_candidate(self):
        self.assertEqual(validate_ts(self.reaction, [self.ts, self.freq])['label'], 'CANDIDATE_TS')

    def test_full_chain_matches_endpoints_in_either_direction(self):
        self.assertEqual(validate_ts(self.reaction, self.jobs)['label'], 'VALIDATED_TS')
        self.end_f['parsed']['final_structure'] = self.reaction.reactant.to_dict()
        self.end_r['parsed']['final_structure'] = self.reaction.product.to_dict()
        self.assertEqual(validate_ts(self.reaction, self.jobs)['label'], 'VALIDATED_TS')

    def test_wrong_endpoints_are_rejected(self):
        self.end_f['parsed']['final_structure'] = self.reaction.reactant.to_dict()
        self.end_f['output_structure'] = self.reaction.reactant.to_dict()
        self.assertEqual(validate_ts(self.reaction, self.jobs)['label'], 'WRONG_CHANNEL')

    def test_wrong_vibration_cannot_pass(self):
        self.freq['parsed']['modes'][0]['displacements'] = [[0, 1, 0], [0, -1, 0], [0, 1, 0]]
        self.assertEqual(validate_ts(self.reaction, self.jobs)['label'], 'WRONG_CHANNEL')

    def test_tiny_or_multiple_negative_frequencies_are_inconclusive(self):
        for values in [[-2, 100, 200, 300], [-500, -2, 200, 300], [100, 200, 300, 400]]:
            jobs = copy.deepcopy(self.jobs)
            jobs[1]['parsed']['frequencies_cm1'] = values
            self.assertNotEqual(validate_ts(self.reaction, jobs)['label'], 'VALIDATED_TS')

    def test_constraints_level_geometry_and_synthetic_provenance(self):
        for key, value in [('constraints', [{'x': 1}]), ('synthetic', True), ('status', 'timed_out')]:
            jobs = copy.deepcopy(self.jobs)
            jobs[0][key] = value
            self.assertNotEqual(validate_ts(self.reaction, jobs)['label'], 'VALIDATED_TS')
        self.freq['level'] = 'other'
        self.assertNotEqual(validate_ts(self.reaction, self.jobs)['label'], 'VALIDATED_TS')
        self.freq['level'] = self.ts['level']
        self.freq['input_structure'] = self.reaction.product.to_dict()
        self.assertNotEqual(validate_ts(self.reaction, self.jobs)['label'], 'VALIDATED_TS')

    def test_unfinished_irc_or_missing_endpoint_convergence_cannot_pass(self):
        self.forward['parsed']['irc_complete'] = False
        self.assertEqual(validate_ts(self.reaction, self.jobs)['label'], 'CANDIDATE_TS')
        self.forward['parsed']['irc_complete'] = True
        self.end_f['parsed']['optimization_converged'] = False
        self.assertNotEqual(validate_ts(self.reaction, self.jobs)['label'], 'VALIDATED_TS')

    def test_linear_structure_requires_3n_minus_5_modes(self):
        self.freq['parsed']['frequencies_cm1'] = [-500., 100., 200.]
        self.freq['parsed']['modes'] = self.freq['parsed']['modes'][:3]
        self.assertNotEqual(validate_ts(self.reaction, self.jobs)['label'], 'VALIDATED_TS')

    def test_mode_projection_uses_frequency_coordinate_frame(self):
        # Gaussian can rotate the structure before printing normal modes.
        data = self.freq['parsed']['final_structure']
        data['positions'] = tuple((0., p[0], 0.) for p in data['positions'])
        self.freq['parsed']['modes'][0]['displacements'] = [[0., 0., 0.], [0., 1., 0.], [0., 0., 0.]]
        self.assertEqual(validate_ts(self.reaction, self.jobs)['label'], 'VALIDATED_TS')

    def test_rotated_frequency_geometry_with_transverse_mode_is_wrong_channel(self):
        data = self.freq['parsed']['final_structure']
        data['positions'] = tuple((0., p[0], 0.) for p in data['positions'])
        self.assertEqual(validate_ts(self.reaction, self.jobs)['label'], 'WRONG_CHANNEL')


if __name__ == '__main__':
    unittest.main()
