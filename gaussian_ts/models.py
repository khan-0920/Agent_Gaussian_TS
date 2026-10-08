"""Stable atom mapping and validated research inputs (Angstrom, degrees)."""
import hashlib
import json
import math
from dataclasses import dataclass
from typing import Tuple

ELEMENTS = ('X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn '
            'Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce '
            'Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn '
            'Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()


def finite(value, name='value'):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
        raise ValueError('%s must be a finite number' % name)
    return float(value)


def integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError('%s must be an integer >= %s' % (name, minimum))
    return value


@dataclass(frozen=True)
class Structure:
    symbols: Tuple[str, ...]
    positions: Tuple[Tuple[float, float, float], ...]
    atom_ids: Tuple[int, ...]
    title: str = 'structure'

    def __post_init__(self):
        object.__setattr__(self, 'symbols', tuple(self.symbols))
        object.__setattr__(self, 'atom_ids', tuple(self.atom_ids))
        object.__setattr__(self, 'positions', tuple(tuple(finite(x, 'coordinate') for x in p) for p in self.positions))
        n = len(self.symbols)
        if not n or n != len(self.positions) or n != len(self.atom_ids):
            raise ValueError('atom, coordinate and mapping counts must agree')
        if any(s not in ELEMENTS[1:] for s in self.symbols) or any(len(p) != 3 for p in self.positions):
            raise ValueError('invalid element or coordinate dimension')
        if len(set(self.atom_ids)) != n:
            raise ValueError('duplicate atom mapping')
        for atom in self.atom_ids:
            integer(atom, 'atom ID', 1)
        if not isinstance(self.title, str) or '\n' in self.title or '\r' in self.title:
            raise ValueError('title must be a single line')

    def index(self, atom):
        integer(atom, 'atom ID', 1)
        try:
            return self.atom_ids.index(atom)
        except ValueError:
            raise ValueError('unknown atom ID: %s' % atom) from None

    def point(self, atom):
        return self.positions[self.index(atom)]

    @property
    def fingerprint(self):
        payload = {'symbols': self.symbols, 'positions': [[round(x, 7) for x in p] for p in self.positions], 'atom_ids': self.atom_ids}
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    def to_dict(self):
        return {'symbols': self.symbols, 'positions': self.positions, 'atom_ids': self.atom_ids, 'title': self.title}

    @classmethod
    def from_dict(cls, data):
        return cls(**data)

    @classmethod
    def from_xyz(cls, text, atom_ids=None):
        lines = text.splitlines()
        if len(lines) < 3:
            raise ValueError('XYZ requires count, title and coordinates')
        try:
            count = int(lines[0])
        except ValueError:
            raise ValueError('invalid XYZ count') from None
        rows = [line.split() for line in lines[2:] if line.strip()]
        if count < 1 or count != len(rows) or any(len(row) != 4 for row in rows):
            raise ValueError('XYZ count or coordinate rows invalid')
        return cls(tuple(r[0] for r in rows), tuple(tuple(float(x) for x in r[1:]) for r in rows),
                   tuple(atom_ids if atom_ids is not None else range(1, count + 1)), lines[1])

    def to_xyz(self):
        rows = ['%s %.10f %.10f %.10f' % (s, *p) for s, p in zip(self.symbols, self.positions)]
        return '\n'.join([str(len(rows)), self.title] + rows) + '\n'


@dataclass(frozen=True)
class Reaction:
    reactant: Structure
    product: Structure
    charge: int
    multiplicity: int
    bond_changes: tuple

    def __post_init__(self):
        if self.reactant.symbols != self.product.symbols or self.reactant.atom_ids != self.product.atom_ids:
            raise ValueError('reactant/product must have identical ordered elements and stable atom mapping')
        if isinstance(self.charge, bool) or not isinstance(self.charge, int):
            raise ValueError('charge must be an integer')
        integer(self.multiplicity, 'multiplicity', 1)
        electrons = sum(ELEMENTS.index(s) for s in self.reactant.symbols) - self.charge
        if electrons < 1 or self.multiplicity > electrons + 1 or (electrons - self.multiplicity + 1) % 2:
            raise ValueError('electron count and multiplicity are inconsistent')
        seen = set()
        for change in self.bond_changes:
            if set(change) != {'atoms', 'change'} or change['change'] not in (-1, 1) or isinstance(change['change'], bool):
                raise ValueError('bond change requires atoms and change (-1 break, +1 form)')
            atoms = change['atoms']
            if len(atoms) != 2 or atoms[0] == atoms[1]:
                raise ValueError('bond requires two distinct atom IDs')
            for a in atoms:
                self.reactant.index(a)
            key = tuple(sorted(atoms))
            if key in seen:
                raise ValueError('duplicate reaction bond')
            seen.add(key)

    def to_dict(self):
        return {'reactant': self.reactant.to_dict(), 'product': self.product.to_dict(), 'charge': self.charge,
                'multiplicity': self.multiplicity, 'bond_changes': self.bond_changes}

    @classmethod
    def from_dict(cls, data):
        return cls(Structure.from_dict(data['reactant']), Structure.from_dict(data['product']),
                   data['charge'], data['multiplicity'], tuple(data['bond_changes']))


@dataclass(frozen=True)
class Limits:
    max_actions: int = 40
    max_jobs: int = 12
    max_wall_seconds: float = 3600.
    max_core_seconds: float = 7200.
    max_cpus: int = 2
    max_memory_mb: int = 2048
    job_timeout_seconds: float = 300.

    def __post_init__(self):
        for name in ('max_actions', 'max_jobs', 'max_cpus', 'max_memory_mb'):
            integer(getattr(self, name), name, 1)
        for name in ('max_wall_seconds', 'max_core_seconds', 'job_timeout_seconds'):
            if finite(getattr(self, name), name) <= 0:
                raise ValueError('%s must be positive' % name)
