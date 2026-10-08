"""Conservative parser for a single Gaussian job; raw log remains authoritative."""
import hashlib
import re
from .models import Structure, ELEMENTS, finite

NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[DEde][-+]?\d+)?'


def number(text):
    return finite(float(text.replace('D', 'E').replace('d', 'e')), 'Gaussian number')


def parse_log(text, atom_ids=None):
    lines = text.splitlines()
    normal_positions = [m.start() for m in re.finditer(r'Normal termination of Gaussian', text)]
    error_positions = [m.start() for m in re.finditer(r'Error termination', text)]
    second_start = any('Entering Gaussian System' in text[pos:] for pos in normal_positions)
    multi = (len(normal_positions) > 1 or '--Link1--' in text or 'Proceeding to internal job step number' in text
             or second_start or text.count('Entering Gaussian System') > 1)
    normal = bool(normal_positions) and (not error_positions or normal_positions[-1] > error_positions[-1]) and not multi
    energies = [number(x) for x in re.findall(r'SCF Done:.*?=\s*(' + NUMBER + ')', text)]
    trajectory, parse_warnings = [], []
    for i, line in enumerate(lines):
        if not re.search(r'(Standard|Input|Z-Matrix) orientation:', line):
            continue
        j, dashes = i + 1, 0
        while j < len(lines) and dashes < 2:
            if re.match(r'\s*-{5,}', lines[j]):
                dashes += 1
            j += 1
        symbols, points, complete = [], [], False
        while j < len(lines):
            if re.match(r'\s*-{5,}', lines[j]):
                complete = True
                break
            row = lines[j].split()
            try:
                if len(row) != 6 or int(row[0]) != len(points) + 1:
                    raise ValueError('coordinate row')
                atomic_number = int(row[1])
                if not 1 <= atomic_number < len(ELEMENTS):
                    raise ValueError('unsupported dummy atom')
                symbols.append(ELEMENTS[atomic_number])
                points.append(tuple(number(x) for x in row[3:]))
            except ValueError:
                break
            j += 1
        if not complete or not symbols or not any('Angstrom' in x for x in lines[i+1:j]):
            parse_warnings.append('incomplete or unsupported orientation block')
            continue
        try:
            mapping = tuple(atom_ids) if atom_ids is not None else tuple(range(1, len(symbols)+1))
            trajectory.append(Structure(tuple(symbols), tuple(points), mapping, 'Gaussian orientation').to_dict())
        except ValueError as exc:
            parse_warnings.append(str(exc))
    convergence = {}
    labels = ('Maximum Force', 'RMS Force', 'Maximum Displacement', 'RMS Displacement')
    # Last table, never mix criteria from different optimization steps.
    tables = re.split(r'Item\s+Value\s+Threshold\s+Converged\?', text)
    if len(tables) > 1:
        table = tables[-1]
        for label in labels:
            m = re.search(re.escape(label) + r'\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s+(YES|NO)', table)
            if m:
                convergence[label] = {'value': number(m[1]), 'threshold': number(m[2]), 'passed': m[3] == 'YES' and number(m[1]) <= number(m[2])}
    converged = normal and 'Optimization completed.' in text and len(convergence) == 4 and all(x['passed'] for x in convergence.values())
    frequencies, modes = [], []
    for i, line in enumerate(lines):
        m = re.search(r'Frequencies\s+--\s+(.*)', line)
        if not m:
            continue
        try:
            values = [number(x) for x in m[1].split()]
        except ValueError:
            parse_warnings.append('invalid frequencies')
            continue
        offset = len(frequencies)
        frequencies.extend(values)
        block = [{'frequency_cm1': value, 'displacements': []} for value in values]
        j = i + 1
        while j < len(lines) and not re.search(r'Atom\s+AN\s+X\s+Y\s+Z', lines[j]):
            if 'Frequencies' in lines[j] or 'termination' in lines[j]:
                break
            j += 1
        if j < len(lines) and re.search(r'Atom\s+AN\s+X\s+Y\s+Z', lines[j]):
            j += 1
            while j < len(lines):
                row = lines[j].split()
                if len(row) != 2 + 3 * len(values):
                    break
                try:
                    if int(row[0]) != len(block[0]['displacements']) + 1:
                        break
                    for k, mode in enumerate(block):
                        mode['displacements'].append([number(x) for x in row[2+3*k:5+3*k]])
                except ValueError:
                    break
                j += 1
        for k, mode in enumerate(block):
            mode['index'] = offset + k
            modes.append(mode)
    diagnostics = []
    for code, patterns in [('SCF_FAILURE', ('Convergence failure', 'SCF has not converged')), ('OPT_NOT_CONVERGED', ('Number of steps exceeded', 'Optimization stopped')),
                           ('COORDINATE_FAILURE', ('Error imposing constraints', 'Linear angle', 'FormBX had a problem')), ('HESSIAN_FAILURE', ('Wrong number of Negative eigenvalues',))]:
        if any(p.lower() in text.lower() for p in patterns):
            diagnostics.append(code)
    if multi:
        diagnostics.append('MULTI_JOB_LOG_UNSUPPORTED')
    if not normal:
        diagnostics.append('ERROR_TERMINATION' if error_positions else 'INCOMPLETE_LOG')
    points = [int(x) for x in re.findall(r'Point Number\s*:?\s*(\d+)', text)]
    irc_complete = normal and bool(re.search(r'(Reaction path calculation complete|IRC calculation complete)', text, re.I)) and len(set(points)) >= 2
    return {'normal_termination': normal, 'multi_job': multi, 'energy_hartree': energies[-1] if energies else None,
            'energies_hartree': energies, 'optimization_converged': converged, 'convergence': convergence,
            'trajectory': trajectory, 'final_structure': trajectory[-1] if trajectory else None,
            'frequencies_cm1': frequencies, 'modes': modes, 'irc_points': points, 'irc_complete': irc_complete,
            'diagnostics': diagnostics, 'parse_warnings': parse_warnings,
            'log_sha256': hashlib.sha256(text.encode()).hexdigest()}
