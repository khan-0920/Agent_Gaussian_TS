"""Build Gaussian input using a small, validated set of scientific options."""
import re
from dataclasses import dataclass, asdict
from .models import integer, finite

KINDS = ('opt', 'frozen_opt', 'scan', 'ts', 'qst2', 'qst3', 'freq', 'irc_forward', 'irc_reverse', 'endpoint_forward', 'endpoint_reverse')
TS_KINDS = ('ts', 'qst2', 'qst3')
OPT_KINDS = ('opt', 'frozen_opt', 'scan') + TS_KINDS + ('endpoint_forward', 'endpoint_reverse')


@dataclass(frozen=True)
class JobSpec:
    kind: str
    method: str = 'HF'
    basis: str = 'STO-3G'
    cpus: int = 2
    memory_mb: int = 1024
    max_cycles: int = 100
    irc_max_points: int = 50
    irc_step_size: int = 10
    scf_xqc: bool = False
    solvent: str = ''
    constraints: tuple = ()

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError('unsupported Gaussian job kind')
        for name in ('method', 'basis'):
            value = getattr(self, name)
            if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9+*(),_\-]*', value):
                raise ValueError('invalid %s; raw route sections are prohibited' % name)
        if self.method.lower() in ('gen', 'genecp') or self.basis.lower() in ('gen', 'genecp'):
            raise ValueError('custom basis sections are not supported')
        for name in ('cpus', 'memory_mb', 'max_cycles', 'irc_max_points', 'irc_step_size'):
            integer(getattr(self, name), name, 1)
        if not isinstance(self.scf_xqc, bool):
            raise ValueError('scf_xqc must be boolean')
        if self.solvent not in ('', 'Water', 'Acetonitrile', 'Methanol', 'Ethanol', 'Toluene', 'Dichloromethane'):
            raise ValueError('unsupported SMD solvent')
        if self.constraints and self.kind not in ('frozen_opt', 'scan'):
            raise ValueError('constraints permitted only in frozen_opt/scan')
        if self.kind in ('frozen_opt', 'scan') and not self.constraints:
            raise ValueError('constrained optimization requires constraints')

    @property
    def level(self):
        return '%s/%s|SMD:%s|Integral:UltraFine' % (self.method.lower(), self.basis.lower(), self.solvent.lower())

    def to_dict(self):
        return asdict(self)


def constraint_lines(structure, spec):
    lines = []
    scan_found = False
    for c in spec.constraints:
        if not isinstance(c, dict) or set(c) - {'type', 'atoms', 'operation', 'steps', 'step_size'}:
            raise ValueError('unsupported constraint fields')
        ctype, operation = c.get('type'), c.get('operation')
        size = {'B': 2, 'A': 3, 'D': 4}.get(ctype)
        atoms = c.get('atoms', [])
        if size is None or len(atoms) != size or len(set(atoms)) != size or operation not in ('F', 'S'):
            raise ValueError('invalid constraint coordinate or operation')
        indices = [structure.index(atom) + 1 for atom in atoms]
        text = '%s %s %s' % (ctype, ' '.join(str(i) for i in indices), operation)
        if operation == 'S':
            if spec.kind != 'scan':
                raise ValueError('scan operation requires scan job')
            steps = integer(c.get('steps'), 'scan steps', 1)
            step_size = finite(c.get('step_size'), 'scan step size')
            if step_size == 0:
                raise ValueError('scan step cannot be zero')
            text += ' %d %g' % (steps, step_size)
            scan_found = True
        elif set(c).intersection({'steps', 'step_size'}):
            raise ValueError('freeze coordinate does not take scan parameters')
        lines.append(text)
    if spec.kind == 'scan' and not scan_found:
        raise ValueError('scan needs at least one scan coordinate')
    return lines


def prepare_input(reaction, structure, spec):
    if structure.symbols != reaction.reactant.symbols or structure.atom_ids != reaction.reactant.atom_ids:
        raise ValueError('structure does not belong to mapped reaction')
    constraints = constraint_lines(structure, spec)
    if spec.kind in TS_KINDS:
        option = {'ts': 'TS', 'qst2': 'QST2', 'qst3': 'QST3'}[spec.kind]
        operation = 'Opt=(%s,CalcFC,Tight,MaxCycles=%d)' % (option, spec.max_cycles)
    elif spec.kind in ('frozen_opt', 'scan'):
        operation = 'Opt=(ModRedundant,Tight,MaxCycles=%d)' % spec.max_cycles
    elif spec.kind == 'freq':
        operation = 'Freq'
    elif spec.kind.startswith('irc_'):
        direction = 'Forward' if spec.kind == 'irc_forward' else 'Reverse'
        operation = 'IRC=(%s,CalcFC,MaxPoints=%d,StepSize=%d)' % (direction, spec.irc_max_points, spec.irc_step_size)
    else:
        operation = 'Opt=(Tight,MaxCycles=%d)' % spec.max_cycles
    route = '#p %s/%s %s Integral=UltraFine NoSymm' % (spec.method, spec.basis, operation)
    if spec.scf_xqc:
        route += ' SCF=XQC'
    if spec.solvent:
        route += ' SCRF=(SMD,Solvent=%s)' % spec.solvent
    sections = [structure]
    if spec.kind == 'qst2':
        sections = [reaction.reactant, reaction.product]
    elif spec.kind == 'qst3':
        sections = [reaction.reactant, reaction.product, structure]
    lines = ['%%nprocshared=%d' % spec.cpus, '%%mem=%dMB' % spec.memory_mb, '%chk=job.chk', route, '']
    for i, s in enumerate(sections):
        lines += ['Gaussian TS research: %s structure %d' % (spec.kind, i+1), '', '%d %d' % (reaction.charge, reaction.multiplicity)]
        lines += ['%s %.10f %.10f %.10f' % (symbol, *p) for symbol, p in zip(s.symbols, s.positions)]
        lines += ['']
        if constraints:
            lines += constraints + ['']
    return '\n'.join(lines) + '\n'
