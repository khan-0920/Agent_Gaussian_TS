"""Rule-based TS evidence evaluation. LLM text is never validation evidence."""
import math
from .gaussian import TS_KINDS
from .geometry import distance, distance_rms, sub, unit, dot, norm, cross
from .models import Structure


def structure(data):
    return Structure.from_dict(data) if data else None


def same_geometry(a, b, tolerance=1e-4):
    return bool(a and b and distance_rms(structure(a), structure(b)) <= tolerance)


def usable(job):
    p = job.get('parsed', {})
    return (job.get('status') == 'completed' and not job.get('synthetic', True) and p.get('normal_termination')
            and not p.get('multi_job') and not p.get('parse_warnings') and p.get('energy_hartree') is not None)


def mode_alignment(reaction, geometry, mode):
    s = structure(geometry)
    displacements = mode.get('displacements', [])
    if not s or len(displacements) != len(s.atom_ids) or not reaction.bond_changes:
        return {'passed': False, 'reason': 'missing complete vibration or reaction coordinates'}
    derivatives, target = [], []
    for change in reaction.bond_changes:
        a, b = change['atoms']
        axis = unit(sub(s.point(b), s.point(a)))
        derivative = dot(axis, sub(displacements[s.index(b)], displacements[s.index(a)]))
        derivatives.append(derivative)
        target.append(-change['change'])
    length = norm(derivatives) * norm(target)
    cosine = abs(dot(derivatives, target))/length if length > 1e-12 else 0.
    # Eigenvector sign is arbitrary; one consistent sign must move all target bonds.
    sign = 1 if dot(derivatives, target) >= 0 else -1
    all_bonds = all(sign*d*t > 1e-4 for d, t in zip(derivatives, target))
    return {'passed': cosine >= .65 and all_bonds, 'absolute_cosine': cosine, 'bond_distance_derivatives': derivatives,
            'reason': 'target-bond projection; qualitative mode check, not mechanism discovery'}


def endpoint_matches(candidate, expected, reaction):
    actual = structure(candidate)
    if actual is None or distance_rms(actual, expected) > .25:
        return False
    return all(abs(distance(actual, *c['atoms']) - distance(expected, *c['atoms'])) <= .2 for c in reaction.bond_changes)


def vibrational_mode_count(geometry):
    s = structure(geometry)
    n = len(s.atom_ids)
    if n == 1:
        return 0
    origin = s.positions[0]
    vectors = [sub(p, origin) for p in s.positions[1:]]
    axis = unit(max(vectors, key=norm))
    linear = all(norm(cross(v, axis)) < 1e-4 for v in vectors)
    return 3*n - (5 if linear else 6)


def validate_ts(reaction, evidence, ts_id=None):
    """Evidence is engine-owned job records. Missing evidence always fails closed."""
    checks, used = [], []

    def check(name, passed, detail):
        checks.append({'name': name, 'passed': bool(passed), 'detail': detail})
        return bool(passed)

    def finish(label):
        return {'label': label, 'checks': checks, 'evidence_jobs': used, 'validated': label == 'VALIDATED_TS',
                'limitations': ['Endpoint identity uses fixed atom mapping, distance matrices and target bond distances.',
                                'No activation energy/free energy is reported; no thermochemical reference-state analysis.']}

    candidates = [j for j in evidence if j.get('kind') in TS_KINDS and (ts_id is None or j['id'] == ts_id)]
    if not candidates:
        check('unconstrained_stationary_point', False, 'no TS optimization evidence')
        return finish('INCONCLUSIVE')
    ts = candidates[-1]
    used.append(ts['id'])
    ts_geometry = ts['parsed'].get('final_structure')
    valid_ts = usable(ts) and not ts.get('constraints') and ts['parsed'].get('optimization_converged') and ts_geometry
    if not check('unconstrained_stationary_point', valid_ts, 'normal end, all four criteria, energy, geometry, unconstrained real calculation'):
        return finish('FAILED_SEARCH' if ts.get('status') in ('failed', 'timed_out') else 'INCONCLUSIVE')
    frequencies = [j for j in evidence if j.get('kind') == 'freq' and j.get('parent_job') == ts['id']
                   and j.get('level') == ts.get('level') and same_geometry(j.get('input_structure'), ts_geometry)]
    if not frequencies:
        check('frequency_same_stationary_point', False, 'need frequency at same mapped geometry and calculation level')
        return finish('INCONCLUSIVE')
    freq = frequencies[-1]
    used.append(freq['id'])
    p = freq['parsed']
    values = p.get('frequencies_cm1', [])
    n = len(reaction.reactant.atom_ids)
    valid_count = len(values) == vibrational_mode_count(p.get('final_structure')) if p.get('final_structure') else False
    if not check('frequency_same_stationary_point', usable(freq) and valid_count and same_geometry(p.get('final_structure'), ts_geometry),
                 'complete vibrational spectrum at TS; same level and geometry'):
        return finish('INCONCLUSIVE')
    negative = [(i, x) for i, x in enumerate(values) if x < 0]
    if not check('single_significant_imaginary_mode', len(negative) == 1 and negative[0][1] <= -30.,
                 {'frequencies_cm1': values, 'significance_threshold_cm1': -30., 'small additional negatives_require_review': True}):
        return finish('INCONCLUSIVE')
    index = negative[0][0]
    modes = [m for m in p.get('modes', []) if m.get('index') == index and m.get('frequency_cm1') == values[index]]
    if not modes or len(modes[0].get('displacements', [])) != n:
        check('target_mode', False, 'missing complete displacement vectors')
        return finish('INCONCLUSIVE')
    alignment = mode_alignment(reaction, p['final_structure'], modes[0])
    if not check('target_mode', alignment['passed'], alignment):
        return finish('WRONG_CHANNEL')
    endpoints = []
    for direction in ('forward', 'reverse'):
        ircs = [j for j in evidence if j.get('kind') == 'irc_' + direction and j.get('parent_job') == ts['id']
                and j.get('level') == ts.get('level') and same_geometry(j.get('input_structure'), ts_geometry)]
        irc = ircs[-1] if ircs else None
        if not check('irc_' + direction, irc and usable(irc) and irc['parsed'].get('irc_complete') and irc['parsed'].get('final_structure'),
                     'completed multi-point path from same TS and level'):
            return finish('CANDIDATE_TS')
        used.append(irc['id'])
        opts = [j for j in evidence if j.get('kind') == 'endpoint_' + direction and j.get('parent_job') == irc['id']
                and j.get('level') == ts.get('level') and same_geometry(j.get('input_structure'), irc['parsed'].get('final_structure'))]
        opt = opts[-1] if opts else None
        if not check('optimized_endpoint_' + direction, opt and usable(opt) and not opt.get('constraints') and opt['parsed'].get('optimization_converged')
                     and opt['parsed'].get('final_structure'), 'unconstrained optimization starting at this IRC endpoint'):
            return finish('CANDIDATE_TS')
        used.append(opt['id'])
        endpoints.append(opt['parsed']['final_structure'])
    paired = ((endpoint_matches(endpoints[0], reaction.reactant, reaction) and endpoint_matches(endpoints[1], reaction.product, reaction))
              or (endpoint_matches(endpoints[1], reaction.reactant, reaction) and endpoint_matches(endpoints[0], reaction.product, reaction)))
    check('reaction_endpoint_identity', paired, 'fixed mapping, pair-distance RMS <=0.25 A, target-bond deviations <=0.20 A; either direction')
    return finish('VALIDATED_TS' if paired else 'WRONG_CHANNEL')
