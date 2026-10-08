"""Synthetic Gaussian-format snippets; NOT chemistry benchmarks or real results."""


def orientation(xs=(0., 1.5, 3.)):
    rows = [' %d 1 0 %.6f 0.000000 0.000000' % (i+1, x) for i, x in enumerate(xs)]
    return '\n'.join([' Standard orientation:', ' ---------------------------------------------------------------------',
                     ' Center Atomic Atomic Coordinates (Angstroms)', ' Number Number Type X Y Z',
                     ' ---------------------------------------------------------------------'] + rows + [' ---------------------------------------------------------------------'])


def log(kind='ts', xs=(0., 1.5, 3.), normal=True):
    text = orientation(xs) + '\n SCF Done: E(UHF) = -1.500000000D+00 A.U. after 10 cycles\n'
    if kind in ('ts', 'opt'):
        text += ' Item Value Threshold Converged?\n Maximum Force 0.000001 0.000450 YES\n RMS Force 0.000001 0.000300 YES\n Maximum Displacement 0.000001 0.001800 YES\n RMS Displacement 0.000001 0.001200 YES\n Optimization completed.\n'
    if kind == 'freq':
        text += ' Frequencies -- -500.0 100.0 200.0\n Red. masses -- 1.0 1.0 1.0\n Atom AN X Y Z X Y Z X Y Z\n 1 1 0 0 0 0.1 0 0 0 0.1 0\n 2 1 1 0 0 0 0 0 0 0 0\n 3 1 0 0 0 -0.1 0 0 0 -0.1 0\n Frequencies -- 300.0\n Atom AN X Y Z\n 1 1 0 0 0.1\n 2 1 0 0 -0.2\n 3 1 0 0 0.1\n'
    if kind == 'irc':
        text += ' Point Number: 1 Path Number: 1\n Point Number: 2 Path Number: 1\n Reaction path calculation complete.\n'
    text += ' Normal termination of Gaussian 16\n' if normal else ' Error termination via Lnk1e\n'
    return text
