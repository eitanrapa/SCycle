#!/usr/bin/env python3
"""Compare SCycle regression outputs by their physics rather than bit for bit.

Builds that differ in compiler, optimization or CPU give results that differ in the last digits,
and adaptive time stepping turns that into different step counts, which h5diff cannot compare.
This compares what should still agree, case by case: the earthquakes (the times at which the
fault's maximum slip rate rises through 1e-3 m/s) and the fault's slip, shear stress and state
variable at the last output.

Usage:
  tools/compare_runs.py REF OTHER [OTHER ...]

REF and each OTHER are directories written by `tools/regress.sh baseline`, with one subdirectory
per case holding its data_1D.h5. For reference, on the Mac the optimized clang and GCC builds
match the debug baseline with the same event counts, onsets within 0.021 yr over 4753 yr (ex1),
1e-6 yr (ex2) and 2.2e-4 yr (ex4s, ex4g). ex4s and ex4g end after a number of steps rather than
at a time, so builds compare their last outputs at slightly different times.
"""
import os
import sys

import h5py
import numpy as np

YEAR = 3.15576e7   # s
V_EVENT = 1e-3     # m/s


def load(d, case):
    with h5py.File(os.path.join(d, case, 'data_1D.h5'), 'r') as f:
        t = f['/time/time1D'][:, 0, 0]
        vmax = np.abs(f['/fault/slipVel'][:, :, 0]).max(axis=1)
        last = {k: f['/fault/' + k][-1, :, 0] for k in ('slip', 'tau', 'psi')}
    onsets = t[np.where((vmax[1:] >= V_EVENT) & (vmax[:-1] < V_EVENT))[0] + 1]
    return t, onsets, last


def reldiff(a, b):
    s = np.max(np.abs(a))
    return np.max(np.abs(a - b)) / s if s > 0 else np.max(np.abs(a - b))


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    ref, others = sys.argv[1], sys.argv[2:]
    cases = sorted(c for c in os.listdir(ref) if os.path.isfile(os.path.join(ref, c, 'data_1D.h5')))
    if not cases:
        sys.exit('no <case>/data_1D.h5 under %s' % ref)
    for case in cases:
        t0, e0, last0 = load(ref, case)
        print('%s: %s has %d outputs to t = %.6g yr, %d events' % (case, ref, t0.size, t0[-1] / YEAR, e0.size))
        for d in others:
            if not os.path.isfile(os.path.join(d, case, 'data_1D.h5')):
                print('  %s: no output' % d)
                continue
            t, e, last = load(d, case)
            line = '  %s: %d outputs to t = %.6g yr, %d events' % (d, t.size, t[-1] / YEAR, e.size)
            if e.size != e0.size:
                line += ' (DIFFERENT NUMBER OF EVENTS)'
            elif e.size:
                line += ', onsets within %.3g yr' % (np.max(np.abs(e - e0)) / YEAR)
            line += '; last output: slip %.1e, tau %.1e, psi %.1e relative difference' % tuple(
                reldiff(last0[k], last[k]) for k in ('slip', 'tau', 'psi'))
            print(line)


if __name__ == '__main__':
    main()
