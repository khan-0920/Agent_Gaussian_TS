"""Offline process double, not Gaussian. Supports only the illustrative H3 case."""
import sys


def orientation(xs):
    rows = [' %d 1 0 %.6f 0.000000 0.000000' % (i+1, x) for i, x in enumerate(xs)]
    return '\n'.join([' Standard orientation:', ' ---------------------------------------------------------------------',
                     ' Center Atomic Atomic Coordinates (Angstroms)', ' Number Number Type X Y Z',
                     ' ---------------------------------------------------------------------'] + rows + [' ---------------------------------------------------------------------'])


def main():
    text = sys.stdin.read()
    route = next(line for line in text.splitlines() if line.startswith('#')).lower()
    xs = [float(row.split()[1]) for row in text.splitlines() if row.startswith('H ')]
    if len(xs) != 3:
        raise ValueError('synthetic demo only accepts H3')
    print('SYNTHETIC DEMO DATA -- NOT A GAUSSIAN CALCULATION')
    if 'opt=(ts,' in route and abs(xs[1] - 1.5) > .01:
        print(orientation(xs))
        print(' Wrong number of Negative eigenvalues')
        print(' Error termination via Lnk1e')
        return 1
    if 'irc=(forward' in route:
        xs = [0, 2.25, 3]
    elif 'irc=(reverse' in route:
        xs = [0, .75, 3]
    print(orientation(xs))
    print(' SCF Done: E(UHF) = -1.500000000 A.U. after 10 cycles')
    if 'opt=' in route:
        print(' Item Value Threshold Converged?')
        print(' Maximum Force 0.000001 0.000450 YES')
        print(' RMS Force 0.000001 0.000300 YES')
        print(' Maximum Displacement 0.000001 0.001800 YES')
        print(' RMS Displacement 0.000001 0.001200 YES')
        print(' Optimization completed.')
    elif ' freq ' in route:
        print(' Frequencies -- -500.0 100.0 200.0')
        print(' Atom AN X Y Z X Y Z X Y Z')
        print(' 1 1 0 0 0 0.1 0 0 0 0.1 0')
        print(' 2 1 1 0 0 0 0 0 0 0 0')
        print(' 3 1 0 0 0 -0.1 0 0 0 -0.1 0')
        print(' Frequencies -- 300.0')
        print(' Atom AN X Y Z')
        print(' 1 1 0 0 0.1')
        print(' 2 1 0 0 -0.2')
        print(' 3 1 0 0 0.1')
    elif 'irc=' in route:
        print(' Point Number: 1 Path Number: 1')
        print(' Point Number: 2 Path Number: 1')
        print(' Reaction path calculation complete.')
    print(' Normal termination of Gaussian 16')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
