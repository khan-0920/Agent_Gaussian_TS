"""Problem framing and compact tool context for Codex/GPT as the researcher."""
import copy
from .models import integer

PROTOCOL = {
    'label': 'VALIDATED_TS',
    'requirements': ['Unconstrained stationary point at the specified calculation level, all optimization criteria passed.',
                     'Complete spectrum with exactly one significant negative curvature, numerical artifacts reviewed.',
                     'Normal-mode displacement matches the specified reaction coordinates.',
                     'Forward/reverse IRC plus optimized endpoints match mapped reactant/product identities.'],
    'unverified_labels': ['CANDIDATE_TS', 'WRONG_CHANNEL', 'FAILED_SEARCH', 'INCONCLUSIVE'],
    'barriers': 'No activation energy/free energy without explicit energy definition, corrections and reference state.'
}


def define_problem(problem, reaction_config=None):
    """Expose known inputs, hypotheses and gaps; do not invent chemical conditions."""
    allowed = {'question', 'fixed_conditions', 'hypotheses', 'uncertainties', 'notes'}
    if not isinstance(problem, dict) or set(problem)-allowed:
        raise ValueError('problem fields: question, fixed_conditions, hypotheses, uncertainties, notes')
    question = problem.get('question', '')
    if not isinstance(question, str):
        raise ValueError('question must be text')
    for key in ('fixed_conditions', 'hypotheses', 'uncertainties'):
        values = problem.get(key, [])
        if not isinstance(values, list) or not all(isinstance(x, str) and x.strip() for x in values):
            raise ValueError('%s must be a list of non-empty strings' % key)
    if not isinstance(problem.get('notes', ''), str):
        raise ValueError('notes must be text')
    config = reaction_config or {}
    required = ('reactant', 'product', 'charge', 'multiplicity', 'atom_ids', 'bond_changes', 'gaussian_command', 'job_defaults', 'limits')
    missing = [key for key in required if key not in config or config[key] is None or (key != 'charge' and not config[key])]
    if not question.strip():
        missing.insert(0, 'question')
    for setting in ('method', 'basis'):
        if not isinstance(config.get('job_defaults'), dict) or not config['job_defaults'].get(setting):
            missing.append('job_defaults.' + setting)
    if 'multiplicity' in config:
        integer(config['multiplicity'], 'multiplicity', 1)
    if 'charge' in config and (isinstance(config['charge'], bool) or not isinstance(config['charge'], int)):
        raise ValueError('charge must be explicitly specified as an integer')
    known = {key: copy.deepcopy(config[key]) for key in required if key in config}
    return {'question': question, 'fixed_conditions': list(problem.get('fixed_conditions', [])),
            'hypotheses': list(problem.get('hypotheses', [])), 'uncertainties': list(problem.get('uncertainties', [])),
            'notes': problem.get('notes', ''), 'specified_inputs': known, 'missing_inputs': missing,
            'ready_for_calculation': not missing, 'success_criteria': copy.deepcopy(PROTOCOL),
            'readiness_scope': 'Input definition completeness only; init validates molecular data, doctor checks executable, real host checks license/environment.'}


def schema(properties=None, required=None):
    return {'type': 'object', 'properties': properties or {}, 'required': required or [], 'additionalProperties': False}


def tool_catalog():
    atom = {'type': 'integer', 'minimum': 1, 'description': 'Stable atom ID, not a mutable row index'}
    moving = {'type': 'array', 'items': atom, 'minItems': 1, 'uniqueItems': True}
    vec = {'type': 'array', 'items': {'type': 'number'}, 'minItems': 3, 'maxItems': 3}
    text = {'type': 'string', 'minLength': 1}
    tools = {}

    def register(name, description, parameters):
        tools[name] = {'description': description, 'parameters': parameters}

    for name, count, value_schema in [('set_distance', 2, {'type': 'number', 'exclusiveMinimum': 0}),
                                      ('set_angle', 3, {'type': 'number', 'exclusiveMinimum': 0, 'exclusiveMaximum': 180}),
                                      ('set_dihedral', 4, {'type': 'number'})]:
        atoms = {'type': 'array', 'items': atom, 'minItems': count, 'maxItems': count, 'uniqueItems': True}
        register(name, 'Edit a coordinate by moving the explicit target fragment; reference atoms stay fixed.',
                 schema({'atoms': atoms, 'value': value_schema, 'moving': moving}, ['atoms', 'value']))
    register('translate_fragment', 'Translate an explicit fragment in Angstrom.', schema({'moving': moving, 'vector': vec}, ['moving', 'vector']))
    register('rotate_fragment', 'Rigidly rotate an explicit fragment using an axis and degrees.',
             schema({'moving': moving, 'origin': vec, 'axis': vec, 'degrees': {'type': 'number'}}, ['moving', 'origin', 'axis', 'degrees']))
    register('inspect_geometry', 'Read mapped coordinates and gross-overlap checks.', schema())
    register('use_structure', 'Switch to a registered structure; preserve provenance and the previous version.', schema({'structure_id': text}, ['structure_id']))
    register('undo_geometry_edit', 'Restore the previous registered current structure.', schema())
    register('save_hypothesis', 'Record an interpretation separately from computed evidence.', schema({'hypothesis': text}, ['hypothesis']))
    register('request_human_review', 'Record a specific question and pause at this checkpoint.', schema({'question': text}, ['question']))
    register('pause', 'Pause at a tool boundary.', schema())
    register('resume', 'Resume a paused session.', schema())
    register('verify_ts', 'Run independent scientific validation from linked raw job evidence; labels cannot be supplied.', schema({'ts_job_id': text}))
    positive_int = {'type': 'integer', 'minimum': 1}
    constraint = schema({'type': {'enum': ['B', 'A', 'D']}, 'atoms': moving, 'operation': {'enum': ['F', 'S']},
                         'steps': positive_int, 'step_size': {'type': 'number'}}, ['type', 'atoms', 'operation'])
    from .gaussian import KINDS
    properties = {'kind': {'type': 'string', 'enum': list(KINDS)}, 'structure_id': text, 'parent_job': text,
                  'method': text, 'basis': text, 'cpus': positive_int, 'memory_mb': positive_int,
                  'max_cycles': positive_int, 'irc_max_points': positive_int, 'irc_step_size': positive_int,
                  'scf_xqc': {'type': 'boolean'}, 'solvent': {'enum': ['', 'Water', 'Acetonitrile', 'Methanol', 'Ethanol', 'Toluene', 'Dichloromethane']},
                  'constraints': {'type': 'array', 'items': constraint}, 'timeout_seconds': {'type': 'number', 'exclusiveMinimum': 0}}
    register('run_gaussian', 'Submit one budgeted local Gaussian job. Freq/IRC/endpoint require matching parent geometry and level. No raw shell or route.',
             schema(properties, ['kind']))
    return tools


def function_tools():
    return [{'type': 'function', 'function': {'name': name, **definition}} for name, definition in tool_catalog().items()]


def export_context(session):
    """Compact context for Codex: scientific state, open questions and executable tools."""
    report = session.evaluate()
    obs = session.observe()
    mapped = [{'atom_id': atom, 'gaussian_index': i+1, 'element': symbol}
              for i, (atom, symbol) in enumerate(zip(obs['current_structure']['atom_ids'], obs['current_structure']['symbols']))]
    jobs = []
    for j in obs['jobs']:
        p = j.get('parsed', {})
        jobs.append({key: j.get(key) for key in ('id', 'kind', 'level', 'status', 'synthetic', 'parent_job', 'input_structure_id',
                                               'output_structure_id', 'directory', 'log_path', 'error')})
        jobs[-1]['feedback'] = {key: p.get(key) for key in ('normal_termination', 'energy_hartree', 'optimization_converged',
                                                          'convergence', 'frequencies_cm1', 'irc_points', 'irc_complete', 'diagnostics', 'parse_warnings')}
        jobs[-1]['feedback']['negative_modes'] = [m for m in p.get('modes', []) if m['frequency_cm1'] < 0]
    return {'role': 'Codex/GPT is the research decision maker; this package supplies problem framing, tools and evidence.',
            'problem_definition': obs.get('problem_definition'), 'reaction': obs['reaction'], 'atom_mapping': mapped,
            'current_structure_id': obs['current_structure_id'], 'current_structure': obs['current_structure'],
            'registered_structures': {key: {'source_job': value['source_job']} for key, value in obs['structures'].items()},
            'reaction_distances': obs['reaction_distances'], 'clashes': obs['clashes'], 'jobs': jobs,
            'history': obs['history'], 'hypotheses': obs['hypotheses'], 'limits': obs['limits'], 'usage': obs['usage'],
            'paused': obs['paused'], 'last_error': obs['last_error'], 'validation': report,
            'success_criteria': copy.deepcopy(PROTOCOL), 'tools': tool_catalog(),
            'action_contract': {'tool': 'catalog key', 'params': 'matching tool parameters', 'reason': 'chemical/numerical rationale'},
            'artifact_root': str(session.store.root)}
