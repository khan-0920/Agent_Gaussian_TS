"""Structure identity, geometry and rigid alignment; no distance-based residue guessing."""
import gzip
import math
from dataclasses import dataclass, replace
from pathlib import Path
from .io import finite
from .spec import atom_id, residue_id

AMINO_ACIDS = dict(zip(
    ('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL').split(),
    'ARNDCQEGHILKMFPSTWYV'))


@dataclass(frozen=True)
class Atom:
    chain: str
    number: int
    residue: str
    name: str
    element: str
    xyz: tuple
    hetero: bool = False

    @property
    def residue_id(self):
        return self.chain + str(self.number)

    @property
    def key(self):
        return self.residue_id + ':' + self.name


def validate_atoms(atoms):
    if not atoms:
        raise ValueError('empty structure')
    if len({a.key for a in atoms}) != len(atoms):
        raise ValueError('duplicate atom identities or unresolved alternate locations')
    residues = {}
    for atom in atoms:
        atom_id(atom.key)
        if len(atom.xyz) != 3:
            raise ValueError('invalid coordinate')
        for value in atom.xyz:
            finite(value, 'coordinate')
        if atom.residue_id in residues and residues[atom.residue_id] != atom.residue:
            raise ValueError('multiple residue identities at one chain/index')
        residues[atom.residue_id] = atom.residue
        if not atom.element:
            raise ValueError('missing element: ' + atom.key)
    return atoms


def load_atoms(path):
    path = Path(path)
    stream = gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()
    with stream:
        if '.cif' in path.suffixes:
            from biotite.structure.io import pdbx
            array = pdbx.get_structure(pdbx.CIFFile.read(stream), model=1, altloc='all')
            if any(str(i) for i in array.ins_code):
                raise ValueError('insertion codes need explicit normalization before this pipeline')
            atoms = [Atom(str(a.chain_id), int(a.res_id), str(a.res_name), str(a.atom_name),
                          str(a.element).upper(), tuple(float(x) for x in a.coord), bool(a.hetero)) for a in array]
        else:
            atoms, model_seen = [], False
            for line in stream:
                if line.startswith('MODEL'):
                    if model_seen:
                        raise ValueError('multi-model PDB requires explicit model selection')
                    model_seen = True
                if line[:6].strip() not in ('ATOM', 'HETATM'):
                    continue
                if line[16:17].strip() or line[26:27].strip():
                    raise ValueError('alternate locations/insertion codes must be resolved explicitly')
                element = line[76:78].strip().upper()
                if not element:
                    raise ValueError('PDB element column is required; atom names cannot identify metals reliably')
                atoms.append(Atom(line[21], int(line[22:26]), line[17:20].strip(), line[12:16].strip(),
                                  element, tuple(float(line[i:i+8]) for i in (30, 38, 46)), line.startswith('HETATM')))
    return validate_atoms(atoms)


def write_pdb(path, atoms):
    validate_atoms(atoms)
    lines = []
    for index, atom in enumerate(atoms, 1):
        if index > 99999 or not -999 <= atom.number <= 9999 or len(atom.name) > 4 or len(atom.residue) > 3:
            raise ValueError('structure cannot be represented losslessly in PDB fields')
        if any(len('%8.3f' % x) > 8 for x in atom.xyz):
            raise ValueError('coordinate exceeds PDB field width')
        record = 'HETATM' if atom.hetero else 'ATOM'
        # Use the conventional one-character element atom-name alignment.
        name = (' ' + atom.name.ljust(3)) if len(atom.element) == 1 and len(atom.name) < 4 else atom.name.rjust(4)
        lines.append('%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s' %
                     (record, index, name, atom.residue, atom.chain, atom.number,
                      *atom.xyz, 1.0, 0.0, atom.element))
    Path(path).write_text('\n'.join(lines) + '\nEND\n')


def atom_lookup(atoms):
    validate_atoms(atoms)
    return {a.key: a for a in atoms}


def sequence(atoms):
    residues = {}
    for atom in atoms:
        if not atom.hetero:
            residues.setdefault(atom.residue_id, atom.residue)
    if not residues:
        raise ValueError('no protein residues')
    if len({r[0] for r in residues}) != 1:
        raise ValueError('this campaign adapter designs one protein chain; normalize multichain targets separately')
    if any(name not in AMINO_ACIDS for name in residues.values()):
        raise ValueError('noncanonical protein residues require an explicit supported model')
    return ''.join(AMINO_ACIDS[r] for r in residues.values()), list(residues)


def distance(a, b):
    return math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))


def angle(a, b, c):
    u, v = [x-y for x, y in zip(a, b)], [x-y for x, y in zip(c, b)]
    norm = math.sqrt(sum(x*x for x in u) * sum(x*x for x in v))
    if norm < 1e-12:
        raise ValueError('angle has coincident atoms')
    return math.degrees(math.acos(max(-1, min(1, sum(x*y for x, y in zip(u, v)) / norm))))


def geometry(atoms, constraints, mapping, state):
    lookup, rows = atom_lookup(atoms), []
    for constraint in constraints:
        if state not in constraint['states']:
            continue
        resolved = []
        missing = []
        for key in constraint['atoms']:
            residue, atom = atom_id(key)
            mapped = mapping.get(residue)
            target = (mapped + ':' + atom) if mapped else None
            if target not in lookup:
                missing.append(key)
            else:
                resolved.append(lookup[target].xyz)
        value, error = None, None
        if missing:
            error = 'missing mapped atoms: ' + ', '.join(missing)
        else:
            try:
                value = distance(*resolved) if constraint['type'] == 'distance' else angle(*resolved)
            except ValueError as exc:
                error = str(exc)
        rows.append({'name': constraint['name'], 'type': constraint['type'], 'value': value,
                     'min': constraint['min'], 'max': constraint['max'], 'error': error,
                     'passed': value is not None and constraint['min'] <= value <= constraint['max']})
    return {'state': state, 'passed': bool(rows) and all(row['passed'] for row in rows),
            'constraints': rows, 'missing_geometry': not rows or any(row['error'] for row in rows)}


def check_identity(atoms, reference, mapping, spec, include_ligands=True):
    current, source = atom_lookup(atoms), atom_lookup(reference)
    required = {r: list(names) for r, names in spec['motif'].items()}
    if include_ligands:
        required.update({ligand['residue']: [a.name for a in reference if a.residue_id == ligand['residue'] and a.element != 'H']
                         for ligand in spec['ligands']})
    targets = [mapping.get(r) for r in required]
    if None in targets or len(set(targets)) != len(targets):
        raise ValueError('residue mapping is missing or not one-to-one')
    for residue, names in required.items():
        if not names:
            raise ValueError('reference residue is empty: ' + residue)
        for name in names:
            original, target = source.get(residue + ':' + name), current.get(mapping[residue] + ':' + name)
            if original is None or target is None:
                raise ValueError('required catalytic/ligand atom missing: ' + residue + ':' + name)
            if original.element != target.element:
                raise ValueError('mapped element mismatch: ' + residue + ':' + name)
            if residue in spec['motif'] and original.residue != target.residue:
                raise ValueError('catalytic residue identity changed: ' + residue)
        if residue not in spec['motif']:
            actual_names = {a.name for a in atoms if a.residue_id == mapping[residue] and a.element != 'H'}
            if actual_names != set(names):
                raise ValueError('ligand atom inventory changed: ' + residue)


def backbone_checks(atoms):
    seq, residues = sequence(atoms)
    lookup = atom_lookup(atoms)
    missing = [r + ':' + a for r in residues for a in ('N', 'CA', 'C', 'O') if r + ':' + a not in lookup]
    breaks = []
    for before, after in zip(residues, residues[1:]):
        c, n = lookup.get(before + ':C'), lookup.get(after + ':N')
        if c and n and not 1.1 <= distance(c.xyz, n.xyz) <= 1.8:
            breaks.append([before, after, distance(c.xyz, n.xyz)])
    return {'passed': not missing and not breaks, 'length': len(seq), 'missing_atoms': missing, 'chain_breaks': breaks}


def aligned_rmsd(reference, mobile):
    """Horn's optimal rigid-fit RMSD using a symmetric Jacobi eigensolver (stdlib)."""
    if len(reference) != len(mobile) or len(reference) < 3:
        raise ValueError('rigid fit requires at least three explicitly matched atoms')
    n = len(reference)
    xmean = [sum(row[i] for row in reference)/n for i in range(3)]
    ymean = [sum(row[i] for row in mobile)/n for i in range(3)]
    x = [[v-m for v, m in zip(row, xmean)] for row in reference]
    y = [[v-m for v, m in zip(row, ymean)] for row in mobile]
    s = [[sum(a[i]*b[j] for a, b in zip(x, y)) for j in range(3)] for i in range(3)]
    xx, xy, xz = s[0]; yx, yy, yz = s[1]; zx, zy, zz = s[2]
    a = [[xx+yy+zz, yz-zy, zx-xz, xy-yx],
         [yz-zy, xx-yy-zz, xy+yx, zx+xz],
         [zx-xz, xy+yx, -xx+yy-zz, yz+zy],
         [xy-yx, zx+xz, yz+zy, -xx-yy+zz]]
    for _ in range(100):
        p, q = max(((i, j) for i in range(4) for j in range(i+1, 4)), key=lambda ij: abs(a[ij[0]][ij[1]]))
        if abs(a[p][q]) < 1e-12:
            break
        theta = 0.5 * math.atan2(2*a[p][q], a[q][q]-a[p][p])
        c, sn = math.cos(theta), math.sin(theta)
        pp, qq, pq = a[p][p], a[q][q], a[p][q]
        for k in range(4):
            if k not in (p, q):
                kp, kq = a[k][p], a[k][q]
                a[k][p] = a[p][k] = c*kp-sn*kq
                a[k][q] = a[q][k] = sn*kp+c*kq
        a[p][p] = c*c*pp-2*c*sn*pq+sn*sn*qq
        a[q][q] = sn*sn*pp+2*c*sn*pq+c*c*qq
        a[p][q] = a[q][p] = 0
    eigenvalue = max(a[i][i] for i in range(4))
    residual = sum(v*v for row in x+y for v in row)-2*eigenvalue
    return math.sqrt(max(0, residual/n))


def backbone_rmsd(designed, predicted):
    seq1, res1 = sequence(designed)
    seq2, res2 = sequence(predicted)
    if seq1 != seq2:
        raise ValueError('prediction sequence does not match designed sequence')
    first, second = atom_lookup(designed), atom_lookup(predicted)
    try:
        x = [first[r + ':CA'].xyz for r in res1]
        y = [second[r + ':CA'].xyz for r in res2]
    except KeyError as exc:
        raise ValueError('missing CA for sequence-aligned structure comparison') from exc
    return aligned_rmsd(x, y), dict(zip(res1, res2))


def align_sequences(first, second):
    """Needleman-Wunsch; identity denominator includes gaps, not only matches."""
    rows, cols = len(first), len(second)
    score = [[0] * (cols+1) for _ in range(rows+1)]
    for i in range(rows+1): score[i][0] = -2*i
    for j in range(cols+1): score[0][j] = -2*j
    for i in range(1, rows+1):
        for j in range(1, cols+1):
            score[i][j] = max(score[i-1][j-1] + (2 if first[i-1] == second[j-1] else -1),
                              score[i-1][j]-2, score[i][j-1]-2)
    i, j, pairs, columns, matches = rows, cols, [], 0, 0
    while i or j:
        columns += 1
        if i and j and score[i][j] == score[i-1][j-1] + (2 if first[i-1] == second[j-1] else -1):
            i -= 1; j -= 1
            pairs.append((i, j))
            matches += first[i] == second[j]
        elif i and score[i][j] == score[i-1][j]-2:
            i -= 1
        else:
            j -= 1
    return {'identity': matches/columns if columns else 0.0, 'pairs': list(reversed(pairs)),
            'coverage': len(pairs)/max(rows, cols) if rows or cols else 0.0}


def compare_designs(first, second):
    seq1, res1 = sequence(first); seq2, res2 = sequence(second)
    alignment = align_sequences(seq1, seq2)
    a, b = atom_lookup(first), atom_lookup(second)
    pairs = [(a[res1[i]+':CA'].xyz, b[res2[j]+':CA'].xyz) for i, j in alignment['pairs']
             if res1[i]+':CA' in a and res2[j]+':CA' in b]
    alignment['ca_rmsd'] = aligned_rmsd([x for x, _ in pairs], [y for _, y in pairs]) if len(pairs) >= 3 else None
    alignment.pop('pairs')
    return alignment
