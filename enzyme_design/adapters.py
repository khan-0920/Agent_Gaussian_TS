"""Model-specific adapters. Heavy libraries are imported only inside their runtime."""
import copy
import csv
import json
import math
import re
import shutil
import sys
from dataclasses import replace
from pathlib import Path
from .io import read_json, write_json, sha256
from .structure import (load_atoms, write_pdb, atom_lookup, sequence as protein_sequence,
                        check_identity, geometry, backbone_checks, backbone_rmsd, compare_designs, AMINO_ACIDS)


def selected_spec(spec, ligand_ids):
    result = copy.deepcopy(spec)
    result['ligands'] = [l for l in result['ligands'] if l['residue'] in ligand_ids]
    return result


def reference_atoms(work):
    if 'prepare' in work.dependencies:
        return load_atoms(work.path(work.dependencies['prepare']['summary']['reference_pdb']))
    return load_atoms(work.path(work.spec['theozyme']))


def normalize_resnames(atoms, spec, mapping):
    names = {mapping[l['residue']]: l['pdb_name'] for l in spec['ligands']}
    return [replace(a, residue=names.get(a.residue_id, a.residue)) for a in atoms]


def prepare(work):
    from rdkit import Chem
    atoms = reference_atoms(work)
    lookup = atom_lookup(atoms)
    identity = {a.residue_id: a.residue_id for a in atoms}
    check_identity(atoms, atoms, identity, work.spec)
    molecules = []
    for ligand in work.spec['ligands']:
        subset = [a for a in atoms if a.residue_id == ligand['residue']]
        if not subset or {a.residue for a in subset} != {ligand['name']}:
            raise ValueError('theozyme ligand residue identity mismatch: ' + ligand['residue'])
        if not all(a.hetero for a in subset):
            raise ValueError('ligand must be explicitly non-protein: ' + ligand['residue'])
        record = {'residue': ligand['residue'], 'name': ligand['name'], 'atoms': len(subset)}
        if 'sdf' in ligand:
            supplier = Chem.SDMolSupplier(str(work.path(ligand['sdf'])), removeHs=False)
            molecules_in_file = list(supplier)
            if len(molecules_in_file) != 1 or molecules_in_file[0] is None:
                raise ValueError('ligand SDF must contain exactly one valid molecule')
            molecule = molecules_in_file[0]
            heavy = {a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum() != 1}
            if set(ligand['atom_map'].values()) != heavy:
                raise ValueError('SDF atom_map must cover exactly every heavy atom')
            current_heavy = {a.name for a in subset if a.element != 'H'}
            if current_heavy != set(ligand['atom_map']):
                raise ValueError('ligand SDF mapping and theozyme atom names differ')
            for name, index in ligand['atom_map'].items():
                if lookup[ligand['residue'] + ':' + name].element != molecule.GetAtomWithIdx(index).GetSymbol().upper():
                    raise ValueError('SDF mapped element mismatch')
            record.update(formal_charge=Chem.GetFormalCharge(molecule),
                          smiles=Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True))
        else:
            record['ccd'] = ligand['ccd']
        molecules.append(record)
    for residue in work.spec['motif']:
        subset = [a for a in atoms if a.residue_id == residue]
        if any(a.hetero for a in subset) or not subset or subset[0].residue not in AMINO_ACIDS:
            raise ValueError('catalytic motif must identify canonical protein residues')
    checks = {state: geometry(atoms, work.spec['constraints'], identity, state)
              for state in ['design'] + [s['name'] for s in work.spec['states']]}
    if any(g['missing_geometry'] for g in checks.values()):
        raise ValueError('theozyme lacks atoms required by explicit geometry constraints')
    if not checks['design']['passed']:
        raise ValueError('theozyme itself violates the declared design geometry; inspect input constraints')
    generation = work.spec['generation']
    ligands = selected_spec(work.spec, generation['ligands'])['ligands']
    fixed = {residue: ','.join(names) for residue, names in work.spec['motif'].items()}
    fixed.update({l['residue']: 'ALL' for l in ligands})
    config = {work.spec['name']: {'input': str(work.path(work.spec['theozyme'])),
                                'ligand': ','.join(l['name'] for l in ligands),
                                'unindex': ','.join(work.spec['motif']),
                                'length': '%d-%d' % tuple(generation['length']), 'select_fixed_atoms': fixed}}
    write_json(work.out / 'rfd3_input.json', config)
    write_json(work.out / 'chemical_inventory.json', molecules)
    write_json(work.out / 'input_geometry.json', checks)
    write_pdb(work.out / 'reference.pdb', normalize_resnames(atoms, work.spec, identity))
    work.manifest([], summary={'rfd3_input': work.relative(work.out / 'rfd3_input.json'),
                              'reference_pdb': work.relative(work.out / 'reference.pdb'),
                              'chemical_inventory': molecules, 'geometry': checks,
                              'evidence': work.spec['evidence'],
                              'validated_ts': False,
                              'evidence_scope': 'External reports are retained as provenance; no TS validation is inferred from supplied status.'})


def generate(work):
    generation = work.spec['generation']
    out = work.out / 'raw'
    command = [work.tool['executable'], 'design',
               'inputs=' + str(work.path(work.previous['summary']['rfd3_input'])), 'out_dir=' + str(out),
               'ckpt_path=' + work.tool['checkpoint'], 'inference_mode=' + generation['inference_mode'],
               'diffusion_batch_size=1', 'n_batches=' + str(generation['backbones']),
               'seed=' + str(generation['seed']), 'prevalidate_inputs=True', 'skip_existing=False',
               'inference_sampler.num_timesteps=' + str(generation['timesteps'])]
    work.run(command)
    paths = sorted(list(out.glob('*.cif')) + list(out.glob('*.cif.gz')))
    if len(paths) != generation['backbones']:
        raise ValueError('RFD3 output count differs from requested backbone count')
    reference = reference_atoms(work)
    spec = selected_spec(work.spec, generation['ligands'])
    candidates, rejected = [], []
    for index, path in enumerate(paths, 1):
        identifier = 'backbone_%05d' % index
        meta_path = path.with_suffix('') if path.suffix == '.gz' else path
        meta = read_json(meta_path.with_suffix('.json'))
        atoms = load_atoms(path)
        mapping = meta.get('diffused_index_map', {})
        if not isinstance(mapping, dict):
            raise ValueError('RFD3 output has no usable residue mapping')
        mapping = dict(mapping)
        identities = {a.residue_id for a in atoms}
        # Unchanged ligand IDs are valid only if explicitly present; no nearest-neighbor guessing.
        for ligand in spec['ligands']:
            if ligand['residue'] not in mapping and ligand['residue'] in identities:
                mapping[ligand['residue']] = ligand['residue']
        try:
            check_identity(atoms, reference, mapping, spec)
            checks = geometry(atoms, spec['constraints'], mapping, 'design')
            topology = backbone_checks(atoms)
            metrics = meta.get('metrics', {})
            upstream_topology = {name: metrics.get(name) for name in ('n_chainbreaks', 'n_conjoined_residues')}
            if any(v is None for v in upstream_topology.values()):
                raise ValueError('RFD3 topology metrics missing')
            if not topology['passed'] or any(v != 0 for v in upstream_topology.values()) or not checks['passed']:
                rejected.append({'id': identifier, 'reason': 'generated topology or catalytic geometry failed',
                                 'geometry': checks, 'topology': topology, 'upstream': upstream_topology})
                continue
            ligand_clashes = metrics.get('n_clashing.ligand_clashes')
            if ligand_clashes is None or ligand_clashes != 0:
                rejected.append({'id': identifier, 'reason': 'RFD3 ligand clash check failed or missing',
                                 'rfd3_metrics': metrics, 'source_structure': work.relative(path)})
                continue
            target = work.out / (identifier + '.pdb')
            write_pdb(target, normalize_resnames(atoms, spec, mapping))
            candidates.append({'id': identifier, 'scaffold_id': identifier, 'structure': work.relative(target),
                               'mapping': mapping, 'fixed_residues': [mapping[r] for r in spec['motif']],
                               'geometry': checks, 'topology': topology, 'rfd3_metrics': metrics,
                               'source_structure': work.relative(path)})
        except ValueError as exc:
            rejected.append({'id': identifier, 'reason': str(exc), 'source_structure': work.relative(path)})
    work.manifest(candidates, rejected, {'requested_backbones': generation['backbones'],
                                       'generated_backbones': len(paths), 'passing_backbones': len(candidates)})


def fasta_records(path):
    records = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            records.append([line[1:], ''])
        elif line.strip():
            if not records:
                raise ValueError('FASTA sequence appears before header')
            records[-1][1] += line.strip()
    return records


def preserve_fixed_residues(packed, scaffold, fixed_residues):
    """Keep input coordinates for fixed residues after torsion-based MPNN reconstruction.

    LigandMPNN's repack_everything=0 preserves torsions, not Cartesian positions.
    The raw packed output remains available; every restoration is explicitly logged.
    """
    fixed = set(fixed_residues)
    original = {a.key: a for a in scaffold if a.residue_id in fixed and a.element != 'H'}
    current = {a.key: a for a in packed if a.residue_id in fixed and a.element != 'H'}
    if set(original) != set(current):
        raise ValueError('fixed residue heavy-atom inventory changed during packing')
    changes = []
    for key, atom in original.items():
        other = current[key]
        if atom.residue != other.residue or atom.element != other.element:
            raise ValueError('fixed residue identity changed during packing')
        changes.append({'atom': key, 'displacement_before_restoration':
                        math.sqrt(sum((x-y)**2 for x, y in zip(atom.xyz, other.xyz)))})
    return [original.get(a.key, a) for a in packed if a.element != 'H'], changes


def sequence(work):
    candidates, rejected = [], []
    per_backbone = work.spec['sequence']['per_backbone']
    reference = reference_atoms(work)
    spec = selected_spec(work.spec, work.spec['generation']['ligands'])
    for scaffold in work.previous['candidates']:
        out = work.out / scaffold['id']
        command = [work.tool['python'], str(Path(work.tool['repo']) / 'run.py'), '--model_type', 'ligand_mpnn',
                   '--checkpoint_ligand_mpnn', work.tool['checkpoint'], '--checkpoint_path_sc', work.tool['packer_checkpoint'],
                   '--pdb_path', str(work.path(scaffold['structure'])), '--out_folder', str(out),
                   '--fixed_residues', ' '.join(scaffold['fixed_residues']), '--batch_size', '1',
                   '--number_of_batches', str(per_backbone), '--temperature', str(work.spec['sequence']['temperature']),
                   '--seed', str(work.spec['generation']['seed']), '--ligand_mpnn_use_atom_context', '1',
                   '--pack_side_chains', '1', '--number_of_packs_per_design', '1', '--repack_everything', '0',
                   '--pack_with_ligand_context', '1', '--force_hetatm', '1']
        work.run(command, cwd=work.tool['repo'])
        fasta = out / 'seqs' / (Path(scaffold['structure']).stem + '.fa')
        records = fasta_records(fasta)
        if len(records) != per_backbone + 1:
            raise ValueError('LigandMPNN FASTA count differs from requested designs plus native record')
        for index, (header, sequence_string) in enumerate(records[1:], 1):
            identifier = scaffold['id'] + '_seq%04d' % index
            packed = out / 'packed' / (Path(scaffold['structure']).stem + '_packed_%d_1.pdb' % index)
            if not packed.is_file():
                raise ValueError('LigandMPNN packed structure missing')
            try:
                atoms = load_atoms(packed)
                observed, _ = protein_sequence(atoms)
                if observed != sequence_string:
                    raise ValueError('packed structure sequence differs from FASTA')
                atoms, restored = preserve_fixed_residues(atoms, load_atoms(work.path(scaffold['structure'])),
                                                         scaffold['fixed_residues'])
                normalized = out / (identifier + '_fixed.pdb')
                write_pdb(normalized, atoms)
                write_json(out / (identifier + '_restored_coordinates.json'), restored)
                check_identity(atoms, reference, scaffold['mapping'], spec)
                checks = geometry(atoms, spec['constraints'], scaffold['mapping'], 'design')
                if not checks['passed']:
                    raise ValueError('fixed catalytic geometry changed during sequence packing')
                candidate = dict(scaffold, id=identifier, sequence=sequence_string, structure=work.relative(normalized),
                                 raw_packed_structure=work.relative(packed),
                                 mpnn_header=header, geometry=checks)
                candidates.append(candidate)
            except ValueError as exc:
                rejected.append({'id': identifier, 'scaffold_id': scaffold['id'], 'reason': str(exc)})
    work.manifest(candidates, rejected, {'input_backbones': len(work.previous['candidates']),
                                       'requested_sequences': len(work.previous['candidates']) * per_backbone,
                                       'passing_sequences': len(candidates)})


def relax(work):
    if not work.previous['candidates']:
        return work.manifest([], summary={'optimized': 0})
    import pyrosetta
    R = pyrosetta.rosetta
    spec = selected_spec(work.spec, work.spec['generation']['ligands'])
    reference = reference_atoms(work)
    params, renaming = [], {}
    for ligand in spec['ligands']:
        names = {a.name for a in reference if a.residue_id == ligand['residue'] and a.element != 'H'}
        if 'sdf' in ligand and 'rosetta_params' not in ligand:
            raise ValueError('custom ligand requires explicit Rosetta params and atom mapping: ' + ligand['residue'])
        if 'rosetta_params' in ligand:
            path = work.path(ligand['rosetta_params'])
            if any(c.isspace() for c in str(path)):
                raise ValueError('Rosetta params paths must not contain whitespace')
            lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
            declared = {line[1] for line in lines if line[0] == 'ATOM'}
            io_names = [line[1] for line in lines if line[0] == 'IO_STRING']
            if len(io_names) != 1:
                raise ValueError('Rosetta params lacks unique IO_STRING')
            mapping = ligand.get('rosetta_atom_map', {name: name for name in names})
            if not names.issubset(mapping) or not {mapping[name] for name in names}.issubset(declared):
                raise ValueError('Rosetta params atom names do not cover mapped ligand heavy atoms')
            renaming[ligand['residue']] = (io_names[0], mapping)
            params.append(str(path))
    options = '-ignore_unrecognized_res false -in:file:load_PDB_components false -no_optH true -packing:flip_HNQ false -constant_seed -jran %d -mute all' % spec['generation']['seed']
    if params:
        options += ' -extra_res_fa ' + ' '.join(params)
    pyrosetta.init(options)
    write_json(work.out / 'rosetta_runtime.json', {'version': pyrosetta.version(), 'options': options,
                                                'score_function': 'ref2015 with explicit geometry constraints',
                                                'scope': 'Structural optimization; scores are REU, not activation energies or affinity.'})
    candidates, rejected = [], []
    for candidate in work.previous['candidates']:
        out = work.out / candidate['id']
        out.mkdir()
        mapping = candidate['mapping']
        before = load_atoms(work.path(candidate['structure']))
        changed = []
        inverse = {}
        for ligand in spec['ligands']:
            residue = mapping[ligand['residue']]
            if ligand['residue'] in renaming:
                name, atoms = renaming[ligand['residue']]
                inverse[residue] = (ligand['pdb_name'], {v: k for k, v in atoms.items()})
        for atom in before:
            source = next((key for key in renaming if mapping[key] == atom.residue_id), None)
            if source is None:
                changed.append(atom)
            else:
                name, atommap = renaming[source]
                if atom.element == 'H':
                    continue  # Hydrogens are rebuilt from the explicitly supplied Rosetta type.
                changed.append(replace(atom, residue=name, name=atommap[atom.name]))
        input_path = out / 'input.pdb'
        write_pdb(input_path, changed)
        pose = pyrosetta.pose_from_pdb(str(input_path))
        for source, restype in spec['relax'].get('residue_types', {}).items():
            residue = mapping[source]
            index = pose.pdb_info().pdb2pose(residue[0], int(residue[1:]))
            if not index:
                raise ValueError('catalytic residue missing from Rosetta pose')
            R.core.pose.replace_pose_residue_copying_existing_coordinates(pose, index,
                                                                          pose.residue_type_set_for_pose().name_map(restype))
        def ros_atom(key):
            source, atom = key.split(':')
            residue = mapping[source]
            atom = renaming.get(source, ('', {}))[1].get(atom, atom)
            index = pose.pdb_info().pdb2pose(residue[0], int(residue[1:]))
            if not index or not pose.residue(index).has(atom):
                raise ValueError('required mapped atom absent in Rosetta pose: ' + key)
            return R.core.id.AtomID(pose.residue(index).atom_index(atom), index)
        score = pyrosetta.get_fa_scorefxn()
        for score_type in (R.core.scoring.coordinate_constraint, R.core.scoring.atom_pair_constraint, R.core.scoring.angle_constraint):
            score.set_weight(score_type, 1.0)
        anchor = R.core.id.AtomID(1, 1)
        constrained = set()
        for i in range(1, pose.total_residue()+1):
            if pose.residue(i).is_protein() and pose.residue(i).has('CA'):
                constrained.add((pose.residue(i).atom_index('CA'), i))
        for source, names in spec['motif'].items():
            for name in names:
                atom = ros_atom(source + ':' + name)
                constrained.add((atom.atomno(), atom.rsd()))
        for ligand in spec['ligands']:
            for atom in reference:
                if atom.residue_id == ligand['residue'] and atom.element != 'H':
                    target = ros_atom(atom.key)
                    constrained.add((target.atomno(), target.rsd()))
        for number, residue in constrained:
            atom = R.core.id.AtomID(number, residue)
            pose.add_constraint(R.core.scoring.constraints.CoordinateConstraint(
                atom, anchor, pose.xyz(atom), R.core.scoring.func.HarmonicFunc(0.0, spec['relax']['coordinate_sd'])))
        for constraint in spec['constraints']:
            if 'design' not in constraint['states']:
                continue
            ids = [ros_atom(key) for key in constraint['atoms']]
            center = (constraint['min'] + constraint['max']) / 2
            sd = max((constraint['max'] - constraint['min']) / 4, 0.01)
            if constraint['type'] == 'distance':
                cst = R.core.scoring.constraints.AtomPairConstraint(*ids, R.core.scoring.func.HarmonicFunc(center, sd))
            else:
                cst = R.core.scoring.constraints.AngleConstraint(*ids, R.core.scoring.func.HarmonicFunc(math.radians(center), math.radians(sd)))
            pose.add_constraint(cst)
        movemap = R.core.kinematics.MoveMap()
        movemap.set_bb(True); movemap.set_chi(True); movemap.set_jump(False)
        for residue in candidate['fixed_residues']:
            index = pose.pdb_info().pdb2pose(residue[0], int(residue[1:]))
            movemap.set_chi(index, False)
        for i in range(1, pose.total_residue()+1):
            if not pose.residue(i).is_protein():
                movemap.set_bb(i, False); movemap.set_chi(i, False)
        initial = float(score(pose))
        relaxer = R.protocols.relax.FastRelax(score, spec['relax']['repeats'])
        relaxer.set_movemap(movemap)
        relaxer.apply(pose)
        final = float(score(pose))
        if not math.isfinite(initial) or not math.isfinite(final):
            raise ValueError('Rosetta returned non-finite scores')
        pose.dump_pdb(str(out / 'raw_relaxed.pdb'))
        atoms = load_atoms(out / 'raw_relaxed.pdb')
        restored = []
        for atom in atoms:
            if atom.residue_id in inverse:
                resname, atommap = inverse[atom.residue_id]
                if atom.element == 'H':
                    continue
                if atom.name not in atommap:
                    raise ValueError('Rosetta introduced an unmapped ligand heavy atom')
                atom = replace(atom, residue=resname, name=atommap[atom.name])
            # Protein/ligand hydrogen naming differs among tools; heavy-atom contracts are explicit.
            if atom.element != 'H':
                restored.append(atom)
        target = out / 'relaxed.pdb'
        write_pdb(target, restored)
        try:
            check_identity(restored, reference, mapping, spec)
            if protein_sequence(restored)[0] != candidate['sequence']:
                raise ValueError('Rosetta changed the designed sequence')
            checks = geometry(restored, spec['constraints'], mapping, 'design')
            if not checks['passed'] or not backbone_checks(restored)['passed']:
                raise ValueError('relaxed geometry/topology failed')
            candidates.append(dict(candidate, structure=work.relative(target), geometry=checks,
                                   rosetta={'initial_REU': initial, 'final_REU': final,
                                            'residue_types': spec['relax'].get('residue_types', {})}))
        except ValueError as exc:
            rejected.append({'id': candidate['id'], 'reason': str(exc), 'structure': work.relative(target)})
    work.manifest(candidates, rejected, {'input_sequences': len(work.previous['candidates']), 'optimized': len(candidates)})


def refold(work):
    if not work.previous['candidates']:
        return work.manifest([], summary={'predicted_sequences': 0, 'states': []})
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from atomworks.io.tools.inference import sdf_to_annotated_atom_array
    from atomworks.io.tools.rdkit import sdf_to_rdkit
    reference = reference_atoms(work)
    candidates, rejected = [], []
    ligands = {l['residue']: l for l in work.spec['ligands']}
    for candidate in work.previous['candidates']:
        designed = load_atoms(work.path(candidate['structure']))
        states = {}
        for state in work.spec['states']:
            out = work.out / candidate['id'] / state['name']
            out.mkdir(parents=True)
            name = candidate['id'] + '_' + state['name']
            components = [{'seq': candidate['sequence'], 'chain_id': 'A'}]
            ligand_chains, atom_names = {}, {}
            if len(state['ligands']) > 25:
                raise ValueError('PDB interchange supports at most 25 non-protein component chains')
            for index, source in enumerate(state['ligands']):
                ligand = ligands[source]
                chain = chr(ord('B') + index)
                ligand_chains[source] = chain
                if 'ccd' in ligand:
                    components.append({'ccd_code': ligand['ccd'], 'chain_id': chain})
                else:
                    molecule = Chem.SDMolSupplier(str(work.path(ligand['sdf'])), removeHs=False)[0]
                    molecule.RemoveAllConformers()
                    molecule = Chem.AddHs(molecule)
                    embed = AllChem.ETKDGv3()
                    embed.randomSeed = work.spec['generation']['seed']
                    if AllChem.EmbedMolecule(molecule, embed) != 0:
                        raise ValueError('could not build independent ligand conformer for RF3')
                    for atom in molecule.GetAtoms():
                        atom.SetAtomMapNum(atom.GetIdx() + 1)
                    path = out / (source + '.sdf')
                    path.write_text(Chem.MolToMolBlock(molecule) + '\n$$$$\n')
                    # Use AtomWorks' unknown-ligand sentinel. A convenient PDB alias
                    # such as LIG can be a real CCD entry and cause chemical replacement.
                    components.append({'path': str(path), 'chain_id': chain, 'res_name': 'UNL'})
                    # Ask the installed RF3/AtomWorks converter for its naming convention.
                    # Pair through explicit SDF indices, never by spatial proximity.
                    arr = sdf_to_annotated_atom_array(str(path), chain_id=chain, res_name='UNL')
                    mol = sdf_to_rdkit(str(path))
                    heavy_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
                    if len(arr) != len(heavy_atoms):
                        raise ValueError('RF3 SDF converter changed atom count; explicit mapping needed')
                    converted = {a.GetAtomMapNum()-1: i for i, a in enumerate(heavy_atoms)}
                    if len(converted) != len(heavy_atoms) or set(converted) != set(ligand['atom_map'].values()):
                        raise ValueError('RF3 SDF converter lost stable atom-map labels')
                    mapping = {}
                    for original, sdf_index in ligand['atom_map'].items():
                        array_index = converted[sdf_index]
                        if str(arr.element[array_index]).upper() != heavy_atoms[array_index].GetSymbol().upper():
                            raise ValueError('RF3 SDF conversion reordered atoms without explicit mapping')
                        mapping[str(arr.atom_name[array_index])] = original
                    if len(mapping) != len(ligand['atom_map']):
                        raise ValueError('RF3 generated non-unique ligand atom names')
                    atom_names[source] = mapping
            write_json(out / 'inputs.json', [{'name': name, 'components': components}])
            write_json(out / 'input_atom_names.json', atom_names)
            command = [work.tool['executable'], 'fold', 'inputs=' + str(out / 'inputs.json'),
                       'out_dir=' + str(out / 'raw'), 'ckpt_path=' + work.tool['checkpoint'],
                       'diffusion_batch_size=1', 'n_recycles=' + str(work.spec['refold']['recycles']),
                       'num_steps=' + str(work.spec['refold']['steps']), 'seed=' + str(work.spec['generation']['seed']),
                       'skip_existing=False']
            work.run(command)
            raw = out / 'raw' / name
            paths = list(raw.glob('*_model.cif'))
            if len(paths) != 1:
                raise ValueError('RF3 must emit exactly one independent predicted model for this state')
            atoms = load_atoms(paths[0])
            normalized = []
            mapped_ligands = {}
            for source, chain in ligand_chains.items():
                subset = [a for a in atoms if a.chain == chain]
                residue_ids = {a.residue_id for a in subset}
                if len(residue_ids) != 1:
                    raise ValueError('RF3 ligand component chain did not retain one residue identity')
                mapped_ligands[source] = next(iter(residue_ids))
            bychain = {chain: source for source, chain in ligand_chains.items()}
            for atom in atoms:
                if atom.chain in bychain:
                    source = bychain[atom.chain]
                    if atom.element == 'H':
                        continue
                    if source in atom_names and atom.name not in atom_names[source]:
                        raise ValueError('RF3 output contains an unmapped custom ligand atom')
                    atom = replace(atom, name=atom_names.get(source, {}).get(atom.name, atom.name),
                                   residue=ligands[source]['pdb_name'], hetero=True)
                normalized.append(atom)
            target = out / 'predicted.pdb'
            write_pdb(target, normalized)
            rmsd, residue_mapping = backbone_rmsd(designed, normalized)
            mapping = {source: residue_mapping[candidate['mapping'][source]] for source in work.spec['motif']}
            mapping.update(mapped_ligands)
            check_identity(normalized, reference, mapping, selected_spec(work.spec, state['ligands']))
            checks = geometry(normalized, work.spec['constraints'], mapping, state['name'])
            confidences = read_json(raw / (name + '_summary_confidences.json'))
            plddt = confidences.get('overall_plddt')
            if isinstance(plddt, bool) or not isinstance(plddt, (int, float)) or not math.isfinite(plddt) or not 0 <= plddt <= 1:
                raise ValueError('RF3 expected finite overall_plddt on its 0-1 scale')
            states[state['name']] = {'structure': work.relative(target), 'mapping': mapping, 'geometry': checks,
                                     'backbone_rmsd': rmsd, 'plddt': plddt, 'confidences': confidences,
                                     'source_structure': work.relative(paths[0]),
                                     'interpretation': state['description']}
        candidates.append(dict(candidate, predictions=states))
    work.manifest(candidates, rejected, {'predicted_sequences': len(candidates),
                                       'states': [state['name'] for state in work.spec['states']],
                                       'prediction_input': 'sequence and chemical components only; no design template coordinates'})


def assess(work):
    if not work.previous['candidates']:
        return work.manifest([], summary={'assessed_sequences': 0, 'samples': 0})
    import random
    import numpy as np
    import torch
    sys.path.insert(0, work.tool['repo'])
    import PLACER
    random.seed(work.spec['generation']['seed'])
    np.random.seed(work.spec['generation']['seed'])
    torch.manual_seed(work.spec['generation']['seed'])
    model = PLACER.PLACER(weights_path=work.tool['checkpoint'])
    reference = reference_atoms(work)
    ligands = {l['residue']: l for l in work.spec['ligands']}
    custom_codes, custom_data, chemical_request = {}, {}, {}
    reserved = set(model.mols())
    for ligand in work.spec['ligands']:
        if 'sdf' not in ligand:
            continue
        code = next(('Q%02d' % i for i in range(100) if 'Q%02d' % i not in reserved), None)
        if code is None:
            raise ValueError('no unambiguous three-character PLACER custom residue code available')
        reserved.add(code)
        custom_codes[ligand['residue']] = code
        chemical_request[code] = {'sdf': str(work.path(ligand['sdf'])), 'atom_map': ligand['atom_map']}
    if chemical_request:
        write_json(work.out / 'chemical_request.json', chemical_request)
        work.run([work.profile['tools']['prepare']['python'], '-m', 'enzyme_design.chemistry',
                  str(work.out / 'chemical_request.json'), str(work.out / 'chemical_result.json')])
        custom_data = read_json(work.out / 'chemical_result.json')['residues']
        for code, data in custom_data.items():
            (work.out / (code + '.sdf')).write_text(data['sdf'])
    write_json(work.out / 'custom_residues.json', custom_data)
    candidates = []
    # Each declared ligand is an independent chemical component. Prevent OpenBabel's
    # distance-based bond perception from turning coordination into undeclared covalent bonds.
    original_parse = PLACER.pdbparser.PDBParser.parse_ligand_from_pdb_to_obmol
    from openbabel import openbabel as ob
    graph_audits = []
    current_label = [None]
    def parse_separate_components(parser, *args, **kwargs):
        molecule = original_parse(parser, *args, **kwargs)
        if molecule is None:
            return molecule
        removed = []
        def identity(atom):
            residue = atom.GetResidue()
            if residue is None:
                raise ValueError('PLACER ligand graph lost residue identity')
            return (residue.GetChain(), residue.GetNum(), residue.GetName())
        for bond in list(ob.OBMolBondIter(molecule)):
            first, second = bond.GetBeginAtom(), bond.GetEndAtom()
            if identity(first) != identity(second):
                removed.append({'atoms': [first.GetIdx(), second.GetIdx()], 'residues': [identity(first), identity(second)],
                                'order': bond.GetBondOrder()})
                molecule.DeleteBond(bond)
        molecule.SetChainsPerceived(True)
        graph_audits.append({'sample_context': current_label[0], 'removed_undeclared_intercomponent_bonds': removed,
                             'scope': 'Declared separate ligand components only; protein-ligand covalent states are not supported.'})
        write_json(work.out / 'ligand_graph_audit.json', graph_audits)
        return molecule
    PLACER.pdbparser.PDBParser.parse_ligand_from_pdb_to_obmol = parse_separate_components
    try:
        for candidate in work.previous['candidates']:
            evaluations = {}
            for state in work.spec['states']:
                if not state.get('assess', True):
                    evaluations[state['name']] = {'status': 'not_requested', 'samples': [], 'geometry_fraction': None}
                    continue
                prediction = candidate['predictions'][state['name']]
                out = work.out / candidate['id'] / state['name']
                out.mkdir(parents=True)
                mapping = prediction['mapping']
                atoms = load_atoms(work.path(prediction['structure']))
                renamed = []
                for atom in atoms:
                    source = next((r for r in state['ligands'] if mapping[r] == atom.residue_id), None)
                    if source in custom_codes:
                        atom = replace(atom, residue=custom_codes[source])
                    renamed.append(atom)
                write_pdb(out / 'input.pdb', renamed)
                inputs = PLACER.PLACERinput()
                inputs.pdb((out / 'input.pdb').read_text())
                inputs.name(candidate['id'] + '_' + state['name'])
                inputs.ignore_ligand_hydrogens(True)
                references, predicted = {}, []
                for source in state['ligands']:
                    ligand = ligands[source]
                    code = custom_codes.get(source, ligand['pdb_name'])
                    if source in custom_codes:
                        if code not in model.mols():
                            inputs.add_custom_residues({code: custom_data[code]})
                        references[code] = str(work.out / (code + '.sdf'))
                    else:
                        if code != ligand['ccd']:
                            raise ValueError('PLACER CCD residue requires identical pdb_name and CCD code')
                        references[code] = 'CCD'
                    residue = mapping[source]
                    if ligand['predict']:
                        predicted.append((residue[0], code, int(residue[1:])))
                if not predicted:
                    raise ValueError('PLACER state must select at least one predicted ligand; apo-only ensemble assessment is not implemented')
                inputs.ligand_reference(references)
                inputs.predict_ligand(predicted)
                inputs.fixed_ligand_noise(0.0)
                if len(predicted) > 1:
                    inputs.predict_multi(True)
                current_label[0] = candidate['id'] + '/' + state['name']
                outputs = model.run(inputs, work.spec['assess']['samples'])
                if len(outputs) != work.spec['assess']['samples']:
                    raise ValueError('PLACER sample count differs from requested count')
                PLACER.protocol.dump_output(outputs, str(out / 'ensemble'))
                samples = []
                for index, output in outputs.items():
                    raw = out / ('sample_%04d_raw.pdb' % index)
                    raw.write_text(output['model'])
                    predicted_atoms = load_atoms(raw)
                    restored = []
                    for atom in predicted_atoms:
                        source = next((r for r in state['ligands'] if mapping[r] == atom.residue_id), None)
                        if source:
                            atom = replace(atom, residue=ligands[source]['pdb_name'])
                        if atom.element != 'H':
                            restored.append(atom)
                    path = out / ('sample_%04d.pdb' % index)
                    write_pdb(path, restored)
                    errors = []
                    try:
                        check_identity(restored, reference, mapping, selected_spec(work.spec, state['ligands']))
                    except ValueError as exc:
                        errors.append(str(exc))
                    checks = geometry(restored, work.spec['constraints'], mapping, state['name'])
                    scores = {}
                    for field in ('prmsd', 'plddt', 'plddt_pde', 'pde'):
                        if field in output:
                            value = float(output[field])
                            if not math.isfinite(value):
                                errors.append('non-finite PLACER score: ' + field)
                            else:
                                scores[field] = value
                    samples.append({'structure': work.relative(path), 'scores': scores, 'geometry': checks,
                                    'identity_errors': errors, 'passed': not errors and checks['passed']})
                fraction = sum(s['passed'] for s in samples) / len(samples)
                evaluations[state['name']] = {'samples': samples, 'geometry_fraction': fraction,
                                              'scope': 'Local predicted conformational geometry, not a free energy or rate estimate.'}
            candidates.append(dict(candidate, ensembles=evaluations))
    finally:
        PLACER.pdbparser.PDBParser.parse_ligand_from_pdb_to_obmol = original_parse
    work.manifest(candidates, summary={'assessed_sequences': len(candidates), 'samples_per_state': work.spec['assess']['samples'],
                                      'states': [s['name'] for s in work.spec['states']]})


def select(work):
    thresholds = work.spec['selection']
    evaluated, rejected = [], []
    for stage, result in work.dependencies.items():
        rejected.extend(dict(row, stage=stage) for row in result.get('rejected', []))
    for candidate in work.previous['candidates']:
        failures, plddts, fractions, rmsds = [], [], [], []
        for state in work.spec['states']:
            name = state['name']
            prediction = candidate.get('predictions', {}).get(name)
            if prediction is None:
                failures.append(name + ': prediction missing')
                continue
            plddt, rmsd = prediction.get('plddt'), prediction.get('backbone_rmsd')
            if not isinstance(plddt, (int, float)) or not math.isfinite(plddt) or not 0 <= plddt <= 1:
                failures.append(name + ': invalid pLDDT')
            else:
                plddts.append(plddt)
                if plddt < thresholds['min_plddt']:
                    failures.append(name + ': pLDDT below threshold')
            if not isinstance(rmsd, (int, float)) or not math.isfinite(rmsd) or rmsd < 0:
                failures.append(name + ': invalid backbone RMSD')
            else:
                rmsds.append(rmsd)
                if rmsd > thresholds['max_backbone_rmsd']:
                    failures.append(name + ': backbone recovery failed')
            if not prediction.get('geometry', {}).get('passed', False):
                failures.append(name + ': predicted catalytic geometry failed')
            if state.get('assess', True):
                ensemble = candidate.get('ensembles', {}).get(name)
                if not ensemble or len(ensemble.get('samples', [])) != work.spec['assess']['samples']:
                    failures.append(name + ': incomplete local ensemble')
                else:
                    # Recount from sample results rather than trusting a supplied summary fraction.
                    fraction = sum(s.get('passed') is True for s in ensemble['samples']) / len(ensemble['samples'])
                    fractions.append(fraction)
                    if fraction < thresholds['min_geometry_fraction']:
                        failures.append(name + ': ensemble geometry fraction below threshold')
        summary = {'min_plddt': min(plddts) if plddts else None,
                   'max_backbone_rmsd': max(rmsds) if rmsds else None,
                   'min_geometry_fraction': min(fractions) if fractions else None,
                   'failures': failures, 'passes_quality_gates': not failures}
        evaluated.append(dict(candidate, evaluation=summary))
    ranked = sorted(evaluated, key=lambda c: (bool(c['evaluation']['failures']),
                                             -(c['evaluation']['min_geometry_fraction'] or 0),
                                             -(c['evaluation']['min_plddt'] or 0),
                                             c['evaluation']['max_backbone_rmsd'] if c['evaluation']['max_backbone_rmsd'] is not None else math.inf,
                                             c['id']))
    selected, comparisons, scaffold_counts = [], [], {}
    for candidate in ranked:
        if candidate['evaluation']['failures']:
            rejected.append({'stage': 'select', 'id': candidate['id'], 'reason': '; '.join(candidate['evaluation']['failures'])})
            continue
        reason = None
        if scaffold_counts.get(candidate['scaffold_id'], 0) >= thresholds['max_per_scaffold']:
            reason = 'per-scaffold diversity limit'
        elif len(selected) >= thresholds['max_candidates']:
            reason = 'candidate delivery limit'
        else:
            atoms = load_atoms(work.path(candidate['structure']))
            for accepted in selected:
                compared = compare_designs(atoms, load_atoms(work.path(accepted['structure'])))
                comparisons.append({'candidate': candidate['id'], 'selected': accepted['id'], **compared})
                if compared['identity'] > thresholds['max_sequence_identity']:
                    reason = 'sequence diversity threshold'
                    break
                if (thresholds['min_structure_rmsd'] > 0 and compared['coverage'] >= 0.7 and
                        compared['ca_rmsd'] is not None and compared['ca_rmsd'] < thresholds['min_structure_rmsd']):
                    reason = 'structural diversity threshold'
                    break
        if reason:
            rejected.append({'stage': 'select', 'id': candidate['id'], 'reason': reason})
        else:
            selected.append(candidate)
            scaffold_counts[candidate['scaffold_id']] = scaffold_counts.get(candidate['scaffold_id'], 0)+1
    (work.out / 'structures').mkdir()
    with (work.out / 'candidates.fasta').open('w') as fasta:
        for candidate in selected:
            fasta.write('>%s scaffold=%s computational_candidate\n%s\n' %
                        (candidate['id'], candidate['scaffold_id'], candidate['sequence']))
            destination = work.out / 'structures' / (candidate['id'] + '.pdb')
            shutil.copy2(work.path(candidate['structure']), destination)
    fields = ['candidate_id', 'scaffold_id', 'selected', 'min_plddt', 'max_backbone_rmsd', 'min_geometry_fraction', 'failure_reasons']
    with (work.out / 'candidate_scores.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields); writer.writeheader()
        selected_ids = {c['id'] for c in selected}
        for candidate in ranked:
            row = candidate['evaluation']
            writer.writerow({'candidate_id': candidate['id'], 'scaffold_id': candidate['scaffold_id'],
                             'selected': candidate['id'] in selected_ids, 'min_plddt': row['min_plddt'],
                             'max_backbone_rmsd': row['max_backbone_rmsd'], 'min_geometry_fraction': row['min_geometry_fraction'],
                             'failure_reasons': '; '.join(row['failures'])})
    with (work.out / 'rejected.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, ['stage', 'id', 'reason']); writer.writeheader()
        for row in rejected:
            writer.writerow({key: row.get(key, '') for key in ('stage', 'id', 'reason')})
    with (work.out / 'experiment_manifest.csv').open('w', newline='') as stream:
        fields = ['candidate_id', 'scaffold_id', 'sequence', 'expression', 'solubility', 'product_identity',
                  'turnover', 'kinetics', 'assay_conditions', 'controls', 'replicates', 'experiment_status']
        writer = csv.DictWriter(stream, fields); writer.writeheader()
        for candidate in selected:
            writer.writerow({'candidate_id': candidate['id'], 'scaffold_id': candidate['scaffold_id'],
                             'sequence': candidate['sequence'], 'experiment_status': 'not_tested'})
    write_json(work.out / 'all_evaluations.json', ranked)
    write_json(work.out / 'diversity_comparisons.json', comparisons)
    write_json(work.out / 'residue_atom_mapping.json', {c['id']: {'design': c['mapping'],
               'states': {name: value['mapping'] for name, value in c['predictions'].items()}} for c in selected})
    funnel = {name: result.get('summary', {}) for name, result in work.dependencies.items()}
    write_json(work.out / 'delivery_manifest.json', {'name': work.spec['name'], 'thresholds': thresholds,
               'evidence': work.spec['evidence'], 'stage_funnel': funnel,
               'selected': [c['id'] for c in selected], 'experimentally_validated': False,
               'input_inventory': read_json(work.root / 'campaign.json')['input_files'],
               'state_descriptions': work.spec['states']})
    work.manifest(selected, rejected, {'status': 'COMPUTATIONAL_CANDIDATES' if selected else 'NO_CANDIDATES_PASSED',
                                      'evaluated_sequences': len(evaluated), 'selected_sequences': len(selected),
                                      'selected_scaffolds': len(scaffold_counts), 'thresholds': thresholds,
                                      'delivery': work.relative(work.out / 'delivery_manifest.json'),
                                      'experimental_activity': 'not measured'})
