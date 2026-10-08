"""Explicit rigid-fragment edits; no implicit bond graph or guessed movable atoms."""
import math
from .models import Structure, finite


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(a, factor):
    return tuple(x * factor for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    length = norm(a)
    if length < 1e-12:
        raise ValueError('undefined axis or overlapping reference atoms')
    return scale(a, 1 / length)


def vector(v):
    if len(v) != 3:
        raise ValueError('vector must have three components')
    return tuple(finite(x, 'vector') for x in v)


def distance(s, a, b):
    return norm(sub(s.point(a), s.point(b)))


def angle(s, a, b, c):
    cosine = dot(unit(sub(s.point(a), s.point(b))), unit(sub(s.point(c), s.point(b))))
    return math.degrees(math.acos(max(-1., min(1., cosine))))


def dihedral(s, a, b, c, d):
    axis = unit(sub(s.point(c), s.point(b)))
    v = sub(s.point(a), s.point(b))
    w = sub(s.point(d), s.point(c))
    v = unit(sub(v, scale(axis, dot(v, axis))))
    w = unit(sub(w, scale(axis, dot(w, axis))))
    return math.degrees(math.atan2(dot(cross(axis, v), w), dot(v, w)))


def rotate(point, origin, axis, degrees):
    axis = unit(axis)
    p = sub(point, origin)
    rad = math.radians(degrees)
    q = add(add(scale(p, math.cos(rad)), scale(cross(axis, p), math.sin(rad))), scale(axis, dot(axis, p)*(1-math.cos(rad))))
    return add(q, origin)


def edit_geometry(s, tool, params):
    allowed = {'set_distance': {'atoms', 'value', 'moving'}, 'set_angle': {'atoms', 'value', 'moving'},
               'set_dihedral': {'atoms', 'value', 'moving'}, 'translate_fragment': {'moving', 'vector'},
               'rotate_fragment': {'moving', 'origin', 'axis', 'degrees'}}
    if tool not in allowed or set(params) - allowed[tool]:
        raise ValueError('unsupported geometry tool or parameters')
    atoms = params.get('atoms', [])
    expected = {'set_distance': 2, 'set_angle': 3, 'set_dihedral': 4}.get(tool, 0)
    if len(atoms) != expected or len(set(atoms)) != len(atoms):
        raise ValueError('invalid geometry reference atoms')
    for a in atoms:
        s.index(a)
    moving = params.get('moving', [atoms[-1]] if atoms else [])
    if not moving or len(set(moving)) != len(moving):
        raise ValueError('moving fragment must contain distinct atom IDs')
    for a in moving:
        s.index(a)
    if atoms and (atoms[-1] not in moving or set(moving).intersection(atoms[:-1])):
        raise ValueError('moving fragment must include target and exclude reference anchors')
    transform = None
    if tool == 'set_distance':
        value = finite(params['value'])
        if value <= 0:
            raise ValueError('distance must be positive')
        delta = scale(unit(sub(s.point(atoms[1]), s.point(atoms[0]))), value - distance(s, *atoms))
        transform = lambda p: add(p, delta)
    elif tool == 'translate_fragment':
        delta = vector(params['vector'])
        transform = lambda p: add(p, delta)
    else:
        if tool == 'rotate_fragment':
            origin, axis = vector(params['origin']), unit(vector(params['axis']))
            degrees = finite(params['degrees'])
        elif tool == 'set_angle':
            value = finite(params['value'])
            if not 0 < value < 180:
                raise ValueError('angle must be between 0 and 180 degrees')
            origin = s.point(atoms[1])
            u, v = sub(s.point(atoms[0]), origin), sub(s.point(atoms[2]), origin)
            axis = cross(u, v)
            if norm(axis) < 1e-10:
                trial = (1., 0., 0.) if abs(unit(u)[0]) < .9 else (0., 1., 0.)
                axis = cross(u, trial)
            degrees = value - angle(s, *atoms)
        else:
            origin = s.point(atoms[1])
            axis = sub(s.point(atoms[2]), origin)
            degrees = finite(params['value']) - dihedral(s, *atoms)
        transform = lambda p: rotate(p, origin, axis, degrees)
    positions = tuple(transform(p) if a in moving else p for a, p in zip(s.atom_ids, s.positions))
    return Structure(s.symbols, positions, s.atom_ids, s.title)


def inspect_clashes(s, minimum_distance=.5):
    # Conservative gross-overlap check, not a valence or chemistry validator.
    minimum_distance = finite(minimum_distance)
    if minimum_distance <= 0:
        raise ValueError('clash threshold must be positive')
    return [{'atoms': [a, b], 'distance_angstrom': distance(s, a, b)}
            for i, a in enumerate(s.atom_ids) for b in s.atom_ids[i+1:] if distance(s, a, b) < minimum_distance]


def distance_rms(a, b):
    """Translation/rotation invariant comparison with fixed mapping (no permutation)."""
    if a.symbols != b.symbols or a.atom_ids != b.atom_ids:
        return float('inf')
    differences = [(distance(a, x, y)-distance(b, x, y))**2 for i, x in enumerate(a.atom_ids) for y in a.atom_ids[i+1:]]
    return math.sqrt(sum(differences)/len(differences)) if differences else norm(sub(a.positions[0], b.positions[0]))
