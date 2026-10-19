"""Explicit scientific inputs and bounded execution settings for enzyme campaigns."""
import copy
import re
from .io import object_fields, safe_name, finite, integer

STAGES = ('prepare', 'generate', 'sequence', 'relax', 'refold', 'assess', 'select')
RESIDUE = re.compile(r'[A-Za-z0-9]-?\d+')
ATOM = re.compile(r'([A-Za-z0-9]-?\d+):([A-Za-z0-9][A-Za-z0-9_\x27*+-]{0,7})')


def residue_id(value):
    if not isinstance(value, str) or not RESIDUE.fullmatch(value):
        raise ValueError('residue must be a single-character chain followed by integer residue number')
    return value[0], int(value[1:])


def atom_id(value):
    if not isinstance(value, str) or not ATOM.fullmatch(value):
        raise ValueError('atom must be an explicit chain/residue:atom identifier, e.g. A12:NE2')
    residue, atom = value.split(':')
    residue_id(residue)
    return residue, atom


def validate_spec(raw):
    spec = copy.deepcopy(raw)
    object_fields(spec, ('schema_version', 'name', 'theozyme', 'evidence', 'motif', 'ligands',
                         'constraints', 'states', 'generation', 'sequence', 'relax', 'refold',
                         'assess', 'selection', 'limits'), label='campaign')
    if spec['schema_version'] != 1:
        raise ValueError('unsupported campaign schema_version')
    safe_name(spec['name'])
    if not isinstance(spec['theozyme'], str) or not spec['theozyme']:
        raise ValueError('theozyme structure path is required')
    evidence = object_fields(spec['evidence'], ('status', 'description'), ('reports',), 'evidence')
    # A claim supplied in JSON cannot promote a geometry to independently validated TS.
    if evidence['status'] not in ('hypothesis', 'external_evidence', 'gaussian_session'):
        raise ValueError('evidence status must be hypothesis, external_evidence or gaussian_session')
    if not isinstance(evidence['description'], str) or not evidence['description'].strip():
        raise ValueError('evidence description is required')
    if 'reports' in evidence and (not isinstance(evidence['reports'], list) or
                                  not all(isinstance(p, str) and p for p in evidence['reports'])):
        raise ValueError('evidence reports must be paths')
    motif = spec['motif']
    if not isinstance(motif, dict) or not motif:
        raise ValueError('explicit catalytic motif residues and fixed atom names are required')
    for residue, atoms in motif.items():
        residue_id(residue)
        if not isinstance(atoms, list) or not atoms or len(set(atoms)) != len(atoms):
            raise ValueError('motif atom lists must be nonempty and unique')
        for atom in atoms:
            atom_id(residue + ':' + atom)
    if not isinstance(spec['ligands'], list) or not spec['ligands']:
        raise ValueError('at least one explicit ligand/cofactor is required')
    ligand_ids = set()
    for ligand in spec['ligands']:
        object_fields(ligand, ('residue', 'name', 'predict'),
                      ('ccd', 'sdf', 'atom_map', 'pdb_name', 'rosetta_params', 'rosetta_atom_map'), 'ligand')
        residue_id(ligand['residue'])
        if ligand['residue'] in ligand_ids or ligand['residue'] in motif:
            raise ValueError('duplicate ligand/motif residue identity')
        ligand_ids.add(ligand['residue'])
        safe_name(ligand['name'], 'ligand name')
        if not isinstance(ligand['predict'], bool):
            raise ValueError('ligand predict must be boolean')
        if ('ccd' in ligand) == ('sdf' in ligand):
            raise ValueError('each ligand needs exactly one CCD code or SDF file')
        if 'ccd' in ligand:
            safe_name(ligand['ccd'], 'CCD code')
        else:
            if not isinstance(ligand['sdf'], str) or not ligand['sdf']:
                raise ValueError('ligand SDF must be a path')
            mapping = ligand.get('atom_map')
            if not isinstance(mapping, dict) or not mapping:
                raise ValueError('custom ligands require atom-name to zero-based SDF index mapping')
            for atom, index in mapping.items():
                atom_id(ligand['residue'] + ':' + atom)
                integer(index, 'SDF atom index', 0)
            if len(set(mapping.values())) != len(mapping):
                raise ValueError('SDF atom mapping must be one-to-one')
        pdb_name = ligand.get('pdb_name', ligand['name'])
        if not isinstance(pdb_name, str) or not re.fullmatch('[A-Za-z0-9]{1,3}', pdb_name):
            raise ValueError('ligands need a 1-3 character pdb_name for downstream tools')
        ligand['pdb_name'] = pdb_name
        if 'rosetta_atom_map' in ligand:
            mapping = ligand['rosetta_atom_map']
            if not isinstance(mapping, dict) or len(set(mapping.values())) != len(mapping):
                raise ValueError('Rosetta atom mapping must be one-to-one')
            for source, target in mapping.items():
                atom_id(ligand['residue'] + ':' + source)
                atom_id(ligand['residue'] + ':' + target)
    if not isinstance(spec['states'], list) or not spec['states']:
        raise ValueError('at least one explicit reaction state is required')
    state_names = set()
    for state in spec['states']:
        object_fields(state, ('name', 'description', 'ligands'), ('assess',), label='state')
        safe_name(state['name'], 'state name')
        if state['name'] in state_names or state['name'] == 'design':
            raise ValueError('state names must be unique; design is reserved')
        state_names.add(state['name'])
        if not isinstance(state['description'], str) or not state['description'].strip():
            raise ValueError('each state needs a chemical interpretation')
        # State-specific chemistry is represented by separately declared ligand identities.
        if not isinstance(state['ligands'], list) or len(set(state['ligands'])) != len(state['ligands']):
            raise ValueError('state ligands must be a unique list of declared residue IDs')
        if set(state['ligands']) - ligand_ids:
            raise ValueError('state refers to undeclared ligand')
        state.setdefault('assess', True)
        if not isinstance(state['assess'], bool):
            raise ValueError('state assess must be boolean')
        if state['assess'] and not any(l['predict'] and l['residue'] in state['ligands'] for l in spec['ligands']):
            raise ValueError('PLACER requires a predicted ligand; set assess=false explicitly for apo-only states')
    if not isinstance(spec['constraints'], list) or not spec['constraints']:
        raise ValueError('explicit catalytic geometry constraints are required')
    for constraint in spec['constraints']:
        object_fields(constraint, ('name', 'type', 'atoms', 'min', 'max', 'states'), label='constraint')
        safe_name(constraint['name'], 'constraint name')
        count = {'distance': 2, 'angle': 3}.get(constraint['type'])
        if count is None or not isinstance(constraint['atoms'], list) or len(constraint['atoms']) != count:
            raise ValueError('constraint atom count does not match distance/angle')
        if len(set(constraint['atoms'])) != count:
            raise ValueError('constraint atoms must be distinct')
        for atom in constraint['atoms']:
            residue, _ = atom_id(atom)
            if residue not in set(motif) | ligand_ids:
                raise ValueError('constraint atoms must refer to declared motif/ligand residues')
        finite(constraint['min'], 'constraint min', 0)
        finite(constraint['max'], 'constraint max', constraint['min'], 180 if count == 3 else None)
        if not isinstance(constraint['states'], list) or not constraint['states'] or set(constraint['states']) - state_names - {'design'}:
            raise ValueError('constraint states must be explicitly declared')
        for state in spec['states']:
            if state['name'] in constraint['states']:
                needed = {atom_id(a)[0] for a in constraint['atoms']} & ligand_ids
                if needed - set(state['ligands']):
                    raise ValueError('constraint needs a ligand absent from its state')
    if len({c['name'] for c in spec['constraints']}) != len(spec['constraints']):
        raise ValueError('constraint names must be unique')
    if not any('design' in c['states'] for c in spec['constraints']):
        raise ValueError('design-stage geometry constraints are required')
    for state in state_names:
        if not any(state in c['states'] for c in spec['constraints']):
            raise ValueError('every predicted state needs explicit geometry constraints')
    gen = object_fields(spec['generation'], ('length', 'backbones', 'seed', 'inference_mode', 'timesteps', 'ligands'), label='generation')
    if not isinstance(gen['ligands'], list) or not gen['ligands'] or len(set(gen['ligands'])) != len(gen['ligands']) or set(gen['ligands']) - ligand_ids:
        raise ValueError('generation ligands must explicitly select declared residue IDs')
    for constraint in spec['constraints']:
        if 'design' in constraint['states']:
            needed = {atom_id(a)[0] for a in constraint['atoms']} & ligand_ids
            if needed - set(gen['ligands']):
                raise ValueError('design constraint requires a ligand not selected for generation')
    if not isinstance(gen['length'], list) or len(gen['length']) != 2:
        raise ValueError('generation length must be [minimum, maximum]')
    integer(gen['length'][0], 'minimum length', 20, 1000)
    integer(gen['length'][1], 'maximum length', gen['length'][0], 1000)
    integer(gen['backbones'], 'backbones', 1, 100000)
    integer(gen['seed'], 'seed', 0, 2**31 - 1)
    integer(gen['timesteps'], 'timesteps', 1, 1000)
    if gen['inference_mode'] not in ('fast', 'legacy'):
        raise ValueError('explicit RFD3 inference_mode must be fast or legacy')
    seq = object_fields(spec['sequence'], ('per_backbone', 'temperature'), label='sequence')
    integer(seq['per_backbone'], 'sequences per backbone', 1, 10000)
    finite(seq['temperature'], 'temperature', 0.001, 10)
    relax = object_fields(spec['relax'], ('repeats', 'coordinate_sd'), ('residue_types',), label='relax')
    integer(relax['repeats'], 'relax repeats', 1, 20)
    finite(relax['coordinate_sd'], 'coordinate constraint SD', 0.001, 10)
    residue_types = relax.get('residue_types', {})
    if not isinstance(residue_types, dict) or set(residue_types) - set(motif):
        raise ValueError('Rosetta residue types must refer to catalytic motif residues')
    for value in residue_types.values():
        if not isinstance(value, str) or not re.fullmatch('[A-Z0-9_]+', value):
            raise ValueError('invalid explicit Rosetta residue type')
    refold = object_fields(spec['refold'], ('recycles', 'steps'), label='refold')
    integer(refold['recycles'], 'RF3 recycles', 1, 100)
    integer(refold['steps'], 'RF3 diffusion steps', 1, 1000)
    assess = object_fields(spec['assess'], ('samples',), label='assess')
    integer(assess['samples'], 'PLACER samples', 1, 10000)
    selection = object_fields(spec['selection'], ('max_candidates', 'max_sequence_identity', 'max_backbone_rmsd',
                                                 'min_plddt', 'min_geometry_fraction'),
                              ('max_per_scaffold', 'min_structure_rmsd'), label='selection')
    integer(selection['max_candidates'], 'max candidates')
    finite(selection['max_sequence_identity'], 'max sequence identity', 0, 1)
    finite(selection['max_backbone_rmsd'], 'max backbone RMSD', 0)
    finite(selection['min_plddt'], 'min pLDDT (0-1)', 0, 1)
    finite(selection['min_geometry_fraction'], 'min geometry fraction', 0, 1)
    selection.setdefault('max_per_scaffold', 1)
    selection.setdefault('min_structure_rmsd', 0.0)
    integer(selection['max_per_scaffold'], 'max per scaffold')
    finite(selection['min_structure_rmsd'], 'min structural diversity RMSD', 0)
    limits = object_fields(spec['limits'], ('max_attempts', 'stage_timeout_seconds', 'total_seconds'), label='limits')
    integer(limits['max_attempts'], 'max attempts', 7)
    finite(limits['stage_timeout_seconds'], 'stage timeout', 1)
    finite(limits['total_seconds'], 'total wall budget', limits['stage_timeout_seconds'])
    return spec


def validate_profile(raw):
    profile = copy.deepcopy(raw)
    object_fields(profile, ('schema_version', 'tools', 'gpu', 'cpus'), ('environment',), 'runtime profile')
    if profile['schema_version'] != 1:
        raise ValueError('unsupported runtime profile schema_version')
    integer(profile['gpu'], 'GPU index', 0)
    integer(profile['cpus'], 'CPUs', 1, 256)
    tools = profile['tools']
    if not isinstance(tools, dict) or set(tools) != set(STAGES):
        raise ValueError('runtime tools must cover exactly the seven campaign stages')
    required = {'prepare': (), 'generate': ('executable', 'checkpoint', 'expected_foundry'),
                'sequence': ('repo', 'checkpoint', 'packer_checkpoint'), 'relax': (),
                'refold': ('executable', 'checkpoint'), 'assess': ('repo', 'checkpoint'), 'select': ()}
    for stage, tool in tools.items():
        object_fields(tool, ('python',) + required[stage], label=stage + ' runtime')
        for key, value in tool.items():
            if not isinstance(value, str) or not value:
                raise ValueError('runtime values must be non-empty strings')
        if stage == 'generate' and tool['expected_foundry'] != '1.0.0':
            raise ValueError('this adapter targets the audited Foundry 1.0.0 API')
    environment = profile.get('environment', {})
    allowed = {'HF_HOME', 'TORCH_HOME', 'FOUNDRY_CHECKPOINT_DIRS', 'TMPDIR', 'RFD3_COMPILE_CACHE_DIR',
               'TORCHINDUCTOR_COMPILE_THREADS', 'TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD'}
    if not isinstance(environment, dict) or set(environment) - allowed:
        raise ValueError('profile environment contains unsupported keys')
    if not all(isinstance(v, str) for v in environment.values()):
        raise ValueError('environment values must be strings')
    profile['environment'] = environment
    return profile
