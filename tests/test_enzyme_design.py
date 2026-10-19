import copy
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from dataclasses import replace
from enzyme_design.adapters import preserve_fixed_residues, select
from enzyme_design.worker import Work
from enzyme_design.campaign import Campaign
from enzyme_design.cli import execute, tools
from enzyme_design.io import write_json, read_json, within
from enzyme_design.spec import STAGES, validate_spec, validate_profile
from enzyme_design.structure import Atom, write_pdb, load_atoms, geometry, check_identity, aligned_rmsd


def specification():
    return {'schema_version': 1, 'name': 'test_contract', 'theozyme': 'theozyme.pdb',
            'evidence': {'status': 'hypothesis', 'description': 'Engineering fixture, not a physical enzyme model.'},
            'motif': {'A1': ['NE2']}, 'ligands': [{'residue': 'Z1', 'name': 'ZN', 'ccd': 'ZN', 'predict': True}],
            'constraints': [{'name': 'coordination', 'type': 'distance', 'atoms': ['A1:NE2', 'Z1:ZN'],
                             'min': 1.8, 'max': 2.3, 'states': ['design', 'holo']}],
            'states': [{'name': 'holo', 'description': 'Metal complex', 'ligands': ['Z1']}],
            'generation': {'length': [80, 100], 'backbones': 1, 'seed': 1, 'timesteps': 200,
                           'inference_mode': 'fast', 'ligands': ['Z1']},
            'sequence': {'per_backbone': 2, 'temperature': 0.1},
            'relax': {'repeats': 1, 'coordinate_sd': 0.2}, 'refold': {'recycles': 3, 'steps': 10},
            'assess': {'samples': 2},
            'selection': {'max_candidates': 2, 'max_sequence_identity': 0.8, 'max_backbone_rmsd': 2.0,
                          'min_plddt': 0.7, 'min_geometry_fraction': 0.5},
            'limits': {'max_attempts': 14, 'stage_timeout_seconds': 5, 'total_seconds': 60}}


def profile():
    tools = {stage: {'python': sys.executable} for stage in STAGES}
    tools['generate'].update(executable=sys.executable, checkpoint='weights', expected_foundry='1.0.0')
    tools['sequence'].update(repo='.', checkpoint='weights', packer_checkpoint='weights')
    tools['refold'].update(executable=sys.executable, checkpoint='weights')
    tools['assess'].update(repo='.', checkpoint='weights')
    return {'schema_version': 1, 'tools': tools, 'gpu': 0, 'cpus': 1}


class EnzymeGeometryTests(unittest.TestCase):
    def setUp(self):
        self.atoms = [Atom('A', 1, 'HIS', 'NE2', 'N', (0., 0., 0.)),
                      Atom('Z', 1, 'ZN', 'ZN', 'ZN', (2.1, 0., 0.), True)]

    def test_pdb_roundtrip_does_not_confuse_atom_name_with_metal_element(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'input.pdb'
            write_pdb(path, self.atoms)
            self.assertEqual(load_atoms(path), self.atoms)

    def test_fixed_cartesian_coordinates_survive_torsion_reconstruction(self):
        displaced = [replace(a, xyz=(1., 2., 3.)) for a in self.atoms]
        atoms, audit = preserve_fixed_residues(displaced, self.atoms, ['A1'])
        self.assertEqual(atoms[0], self.atoms[0])
        self.assertEqual(atoms[1], displaced[1])
        self.assertGreater(audit[0]['displacement_before_restoration'], 0)
        with self.assertRaisesRegex(ValueError, 'inventory'):
            preserve_fixed_residues(displaced[1:], self.atoms, ['A1'])
        with self.assertRaisesRegex(ValueError, 'identity'):
            preserve_fixed_residues([replace(displaced[0], residue='GLN'), displaced[1]], self.atoms, ['A1'])

    def test_missing_or_wrong_mapping_never_passes_catalytic_geometry(self):
        constraints = specification()['constraints']
        self.assertTrue(geometry(self.atoms, constraints, {'A1': 'A1', 'Z1': 'Z1'}, 'holo')['passed'])
        result = geometry(self.atoms, constraints, {'A1': 'A9', 'Z1': 'Z1'}, 'holo')
        self.assertFalse(result['passed'])
        self.assertTrue(result['missing_geometry'])
        self.assertFalse(geometry(self.atoms, constraints, {}, 'absent_state')['passed'])

    def test_fixed_residue_mutation_and_ligand_atom_loss_are_rejected(self):
        spec = specification()
        check_identity(self.atoms, self.atoms, {'A1': 'A1', 'Z1': 'Z1'}, spec)
        mutant = [Atom('A', 1, 'GLN', 'NE2', 'N', (0, 0, 0)), self.atoms[1]]
        with self.assertRaisesRegex(ValueError, 'identity changed'):
            check_identity(mutant, self.atoms, {'A1': 'A1', 'Z1': 'Z1'}, spec)
        with self.assertRaisesRegex(ValueError, 'missing'):
            check_identity(self.atoms[:1], self.atoms, {'A1': 'A1', 'Z1': 'Z1'}, spec)

    def test_rigid_alignment_removes_rotation_translation_but_not_reflection(self):
        reference = [(0, 0, 0), (1, 0, 0), (0, 2, 0), (0, 0, 3)]
        rotated = [(10-y, 20+x, 30+z) for x, y, z in reference]
        reflected = [(-x, y, z) for x, y, z in reference]
        self.assertAlmostEqual(aligned_rmsd(reference, rotated), 0, places=6)
        self.assertGreater(aligned_rmsd(reference, reflected), 0.1)

    def test_duplicate_atoms_and_nonfinite_coordinates_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'input.pdb'
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                write_pdb(path, self.atoms + self.atoms[:1])
            with self.assertRaisesRegex(ValueError, 'finite'):
                write_pdb(path, [Atom('A', 1, 'HIS', 'NE2', 'N', (float('nan'), 0, 0))])


class EnzymeContractTests(unittest.TestCase):
    def test_valid_spec_and_new_version_profile(self):
        self.assertEqual(validate_spec(specification())['ligands'][0]['pdb_name'], 'ZN')
        self.assertEqual(validate_profile(profile())['tools']['generate']['expected_foundry'], '1.0.0')

    def test_old_rfd3_and_unbounded_resource_inputs_rejected(self):
        old = profile(); old['tools']['generate']['expected_foundry'] = '0.2.0'
        with self.assertRaises(ValueError): validate_profile(old)
        for field, value in [('backbones', True), ('timesteps', 0), ('seed', -1)]:
            spec = specification(); spec['generation'][field] = value
            with self.assertRaises(ValueError): validate_spec(spec)
        spec = specification(); spec['limits']['stage_timeout_seconds'] = math.inf
        with self.assertRaises(ValueError): validate_spec(spec)

    def test_unmapped_custom_ligand_and_absent_state_components_rejected(self):
        spec = specification(); ligand = spec['ligands'][0]
        del ligand['ccd']; ligand['sdf'] = 'ligand.sdf'
        with self.assertRaisesRegex(ValueError, 'mapping'): validate_spec(spec)
        spec = specification(); spec['states'][0]['ligands'] = []; spec['states'][0]['assess'] = False
        with self.assertRaisesRegex(ValueError, 'absent'): validate_spec(spec)

    def test_apo_assessment_requires_explicit_opt_out_and_valid_geometry(self):
        spec = specification(); spec['states'][0]['ligands'] = []
        with self.assertRaisesRegex(ValueError, 'assess=false'): validate_spec(spec)
        spec['states'][0]['assess'] = False
        spec['constraints'][0]['states'] = ['design']
        spec['constraints'].append({'name': 'apo_internal', 'type': 'distance',
                                   'atoms': ['A1:NE2', 'A1:CG'], 'min': 1, 'max': 4, 'states': ['holo']})
        spec['motif']['A1'].append('CG')
        self.assertFalse(validate_spec(spec)['states'][0]['assess'])

    def test_external_success_labels_cannot_bypass_scientific_validation(self):
        spec = specification(); spec['evidence']['status'] = 'VALIDATED_TS'
        with self.assertRaises(ValueError): validate_spec(spec)
        config = profile(); config['environment'] = {'LD_PRELOAD': '/arbitrary/library'}
        with self.assertRaises(ValueError): validate_profile(config)


class EnzymeSelectionTests(unittest.TestCase):
    def test_confidence_cannot_override_failed_geometry_or_missing_ensemble(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root/'selection').mkdir()
            write_json(root/'campaign.json', {'input_files': {}})
            baseline = {'id': 'candidate', 'scaffold_id': 'backbone', 'sequence': 'AAA',
                        'structure': 'not_used.pdb', 'mapping': {},
                        'predictions': {'holo': {'plddt': .99, 'backbone_rmsd': .1,
                                                'geometry': {'passed': False}, 'mapping': {}}},
                        'ensembles': {'holo': {'geometry_fraction': 1., 'samples': [{'passed': False}, {'passed': False}]}}}
            missing = copy.deepcopy(baseline); missing['id'] = 'missing'
            missing['predictions']['holo']['geometry']['passed'] = True
            missing['ensembles'] = {}
            work = Work({'root': str(root), 'stage': 'select', 'attempt_dir': 'selection',
                         'spec': validate_spec(specification()), 'profile': validate_profile(profile()),
                         'previous': {'candidates': [baseline, missing]}, 'dependencies': {}})
            select(work)
            result = read_json(root/'selection/result.json')
            self.assertEqual(result['summary']['status'], 'NO_CANDIDATES_PASSED')
            self.assertEqual(result['candidates'], [])
            reasons = {c['id']: c['reason'] for c in result['rejected']}
            self.assertIn('predicted catalytic geometry failed', reasons['candidate'])
            self.assertIn('ensemble geometry fraction below', reasons['candidate'])
            self.assertIn('incomplete local ensemble', reasons['missing'])
            self.assertEqual((root/'selection/candidates.fasta').read_text(), '')


SUCCESS_WORKER = '''
import json,sys,pathlib
req=json.load(open(sys.argv[1]))
out=pathlib.Path(req['root'])/req['attempt_dir']
(out/'payload.txt').write_text('engineering test fixture')
(out/'result.json').write_text(json.dumps({'schema_version':1,'stage':req['stage'],'candidates':[],
 'artifacts':[req['attempt_dir']+'/payload.txt'],'summary':{'engineering_fixture':True}}))
'''


class EnzymeCampaignTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        write_pdb(self.root/'theozyme.pdb', [Atom('A', 1, 'HIS', 'NE2', 'N', (0, 0, 0)),
                                           Atom('Z', 1, 'ZN', 'ZN', 'ZN', (2.1, 0, 0), True)])
        (self.root/'weights').write_bytes(b'fixture-not-model')
        write_json(self.root/'spec.json', specification())
        write_json(self.root/'profile.json', profile())
        self.campaign = Campaign.create(self.root/'spec.json', self.root/'profile.json', self.root/'run')

    def tearDown(self):
        self.temp.cleanup()

    def run_fixture(self, script=SUCCESS_WORKER, retry=False):
        real_popen = subprocess.Popen
        def launch(command, **kwargs):
            return real_popen([sys.executable, '-c', script, command[-1]], **kwargs)
        with patch('enzyme_design.campaign.subprocess.Popen', side_effect=launch):
            return self.campaign.run_stage('prepare', 'test controller behavior with explicitly synthetic fixture', retry=retry)

    def test_dependency_order_immutability_and_duplicate_execution(self):
        with self.assertRaisesRegex(ValueError, 'dependency'):
            self.campaign.run_stage('generate', 'must not run before preparation')
        self.run_fixture()
        self.campaign.run_stage('prepare', 'reuse unchanged completed evidence')
        self.assertEqual(self.campaign.context()['usage']['attempts'], 1)
        (self.campaign.root/'steps/prepare/0001/payload.txt').write_text('tampered')
        with self.assertRaisesRegex(ValueError, 'artifact changed'):
            self.campaign.context()

    def test_failure_requires_explicit_retry_and_keeps_original_logs(self):
        with self.assertRaisesRegex(ValueError, 'code 3'):
            self.run_fixture('import sys; print("known failure"); sys.exit(3)')
        with self.assertRaisesRegex(ValueError, 'explicitly request retry'):
            self.campaign.run_stage('prepare', 'not implicitly retrying')
        self.run_fixture(retry=True)
        state = read_json(self.campaign.state_path)
        self.assertEqual([a['status'] for a in state['stages']['prepare']['attempts']], ['failed', 'completed'])
        self.assertIn('known failure', (self.campaign.root/'steps/prepare/0001/worker.log').read_text())

    def test_timeout_stops_worker_records_failure_and_retains_runtime(self):
        config = specification(); config['limits']['stage_timeout_seconds'] = 1
        write_json(self.root/'short.json', config)
        self.campaign = Campaign.create(self.root/'short.json', self.root/'profile.json', self.root/'short_run')
        with self.assertRaises(subprocess.TimeoutExpired):
            self.run_fixture('import time; time.sleep(60)')
        attempt = read_json(self.campaign.state_path)['stages']['prepare']['attempts'][0]
        self.assertEqual(attempt['status'], 'failed')
        with self.assertRaises(ProcessLookupError): os.kill(attempt['pid'], 0)
        self.assertTrue((self.campaign.root/'steps/prepare/0001/runtime/enzyme_design/worker.py').is_file())

    def test_invalid_worker_artifact_cannot_mark_stage_completed(self):
        script = SUCCESS_WORKER.replace("req['attempt_dir']+'/payload.txt'", "'../foreign.txt'")
        with self.assertRaisesRegex(ValueError, 'foreign artifact'):
            self.run_fixture(script)
        self.assertEqual(self.campaign.context()['stages']['prepare']['status'], 'failed')

    def test_input_edits_and_output_path_escape_are_not_accepted(self):
        data = read_json(self.campaign.root/'spec.json'); data['generation']['seed'] = 100
        write_json(self.campaign.root/'spec.json', data)
        with self.assertRaisesRegex(ValueError, 'input changed'): self.campaign.context()
        with self.assertRaisesRegex(ValueError, 'escapes'): within(self.root, '../outside')

    def test_lock_prevents_concurrent_submission_and_pause_stops_new_jobs(self):
        with self.campaign.lock():
            with self.assertRaisesRegex(ValueError, 'another process'):
                self.campaign.run_stage('prepare', 'duplicate worker must not launch')
        self.campaign.set_paused(True, 'user checkpoint')
        with self.assertRaisesRegex(ValueError, 'paused'):
            self.campaign.run_stage('prepare', 'must honor pause')

    def test_typed_agent_actions_cannot_inject_commands(self):
        self.assertIn('run_stage', tools())
        with self.assertRaises(ValueError):
            execute(self.campaign, {'tool': 'run_stage', 'params': {'stage': 'prepare', 'shell': 'anything'}, 'reason': 'invalid'})
        with self.assertRaises(ValueError):
            execute(self.campaign, {'tool': 'run_pipeline', 'params': {'retry': 'yes'}, 'reason': 'invalid'})

    def test_symlinked_virtualenv_interpreter_is_not_resolved_to_base_python(self):
        link = self.root/'venv_python'; link.symlink_to(sys.executable)
        config = profile()
        config['tools']['prepare']['python'] = str(link)
        write_json(self.root/'profile2.json', config)
        campaign = Campaign.create(self.root/'spec.json', self.root/'profile2.json', self.root/'run2')
        self.assertEqual(read_json(campaign.root/'profile.json')['tools']['prepare']['python'], str(link))


if __name__ == '__main__':
    unittest.main()
