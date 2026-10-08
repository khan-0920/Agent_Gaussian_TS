import math
import unittest

from gaussian_ts.models import Structure, Reaction
from gaussian_ts.geometry import edit_geometry, distance, angle, dihedral, inspect_clashes


class GeometryTests(unittest.TestCase):
    def setUp(self):
        self.s = Structure.from_xyz('4\ntest\nC 0 0 0\nC 1 0 0\nH 1 1 0\nH 1 1 1\n', atom_ids=[10, 20, 30, 40])

    def test_xyz_roundtrip_preserves_mapping(self):
        restored = Structure.from_dict(self.s.to_dict())
        self.assertEqual(restored.atom_ids, (10, 20, 30, 40))
        self.assertEqual(restored.positions, self.s.positions)

    def test_xyz_rejects_bad_counts_nan_and_duplicate_mapping(self):
        for xyz, ids in [('2\nx\nH 0 0 0\n', None), ('1\nx\nH nan 0 0\n', None), ('2\nx\nH 0 0 0\nH 1 0 0\n', [1, 1])]:
            with self.subTest(xyz=xyz), self.assertRaises(ValueError):
                Structure.from_xyz(xyz, atom_ids=ids)

    def test_fragment_distance_edit_preserves_internal_geometry(self):
        s = edit_geometry(self.s, 'set_distance', {'atoms': [10, 20], 'value': 2.0, 'moving': [20, 30, 40]})
        self.assertAlmostEqual(distance(s, 10, 20), 2.0)
        self.assertAlmostEqual(distance(s, 30, 40), 1.0)
        self.assertEqual(s.positions[0], (0., 0., 0.))
        self.assertAlmostEqual(distance(self.s, 10, 20), 1.0)

    def test_angle_and_dihedral_targets(self):
        a = edit_geometry(self.s, 'set_angle', {'atoms': [10, 20, 30], 'value': 120., 'moving': [30, 40]})
        self.assertAlmostEqual(angle(a, 10, 20, 30), 120., places=6)
        d = edit_geometry(self.s, 'set_dihedral', {'atoms': [10, 20, 30, 40], 'value': -60., 'moving': [40]})
        self.assertAlmostEqual(dihedral(d, 10, 20, 30, 40), -60., places=6)

    def test_translate_and_rotate_are_rigid(self):
        s = edit_geometry(self.s, 'translate_fragment', {'moving': [30, 40], 'vector': [2, 0, 0]})
        s = edit_geometry(s, 'rotate_fragment', {'moving': [30, 40], 'origin': [0, 0, 0], 'axis': [0, 0, 1], 'degrees': 90})
        self.assertAlmostEqual(distance(s, 30, 40), 1.)
        self.assertAlmostEqual(s.positions[2][0], -1.)
        self.assertAlmostEqual(s.positions[2][1], 3.)

    def test_invalid_atoms_and_moving_anchor_rejected(self):
        for p in [{'atoms': [10, 999], 'value': 2}, {'atoms': [10, 20], 'value': 2, 'moving': [10, 20]}, {'atoms': [10, 20], 'value': -1}]:
            with self.subTest(params=p), self.assertRaises(ValueError):
                edit_geometry(self.s, 'set_distance', p)

    def test_clash_is_detected(self):
        s = Structure.from_xyz('2\nx\nH 0 0 0\nH 0.1 0 0\n')
        self.assertTrue(inspect_clashes(s))

    def test_reaction_rejects_atom_order_and_invalid_spin(self):
        r = Structure.from_xyz('2\nx\nH 0 0 0\nH 1 0 0\n')
        with self.assertRaises(ValueError):
            Reaction(r, r, 0, 2, ({'atoms': [1, 2], 'change': 1},))
        with self.assertRaises(ValueError):
            Reaction(r, Structure.from_xyz('2\nx\nH 0 0 0\nHe 1 0 0\n'), 0, 1, ())


if __name__ == '__main__':
    unittest.main()
