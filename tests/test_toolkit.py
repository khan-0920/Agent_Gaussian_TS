import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from gaussian_ts.toolkit import define_problem, export_context, tool_catalog
from gaussian_ts.engine import Session
from tests.test_engine import reaction, action


class ToolkitTests(unittest.TestCase):
    def test_problem_definition_keeps_assumptions_separate_and_lists_missing_inputs(self):
        result = define_problem({'question': '寻找指定质子转移的目标过渡态', 'hypotheses': ['单步协同转移'],
                                 'fixed_conditions': ['保持指定自旋态'], 'uncertainties': ['是否需显式溶剂']})
        self.assertFalse(result['ready_for_calculation'])
        self.assertIn('charge', result['missing_inputs'])
        self.assertEqual(result['hypotheses'], ['单步协同转移'])
        self.assertIn('VALIDATED_TS', result['success_criteria']['label'])

    def test_full_problem_config_has_no_silent_defaults(self):
        config = {'reactant': 'r.xyz', 'product': 'p.xyz', 'charge': 0, 'multiplicity': 2, 'atom_ids': [1, 2, 3],
                  'bond_changes': list(reaction().bond_changes), 'gaussian_command': ['g16'],
                  'job_defaults': {'method': 'UHF', 'basis': 'STO-3G'}, 'limits': {'max_jobs': 5}}
        result = define_problem({'question': '验证氢交换路径'}, config)
        self.assertTrue(result['ready_for_calculation'])
        self.assertEqual(result['specified_inputs']['job_defaults']['method'], 'UHF')
        self.assertEqual(result['question'], '验证氢交换路径')

    def test_context_exports_feedback_and_tools_without_raw_trajectory(self):
        with tempfile.TemporaryDirectory() as temp:
            session = Session.create(Path(temp)/'run', reaction())
            session.execute(action('set_distance', {'atoms': [1, 2], 'value': .1}))
            context = export_context(session)
            self.assertIn('collision', context['last_error'].lower().replace('overlap', 'collision'))
            self.assertIn('set_distance', context['tools'])
            self.assertEqual(context['atom_mapping'][0], {'atom_id': 1, 'gaussian_index': 1, 'element': 'H'})
            self.assertEqual(context['validation']['label'], 'INCONCLUSIVE')
            self.assertNotIn('gaussian_command', context)

    def test_catalog_exposes_validatable_action_parameters(self):
        schema = tool_catalog()['set_distance']['parameters']
        self.assertEqual(schema['properties']['atoms']['minItems'], 2)
        self.assertFalse(schema['additionalProperties'])
        self.assertIn('run_gaussian', tool_catalog())

    def test_cli_problem_and_tool_catalog(self):
        project = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temp:
            problem = Path(temp)/'problem.json'
            problem.write_text(json.dumps({'question': '明确研究问题'}))
            for args in [['define', '--problem', str(problem)], ['tools']]:
                result = subprocess.run([sys.executable, '-m', 'gaussian_ts']+args, cwd=project, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIsInstance(json.loads(result.stdout), dict)
