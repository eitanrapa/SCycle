#!/usr/bin/env python3
"""Stage 5 lock-in search: can one fault take over for good, so that the other never ruptures again?

No run of the first two batches or the matrix locks in. While one fault dominates, its root takes
82-95% of the deep creep (batch 2), and the quiet fault fails again within one to three thousand
years. If the dominant fault took all the motion, the quiet fault's stress would have no secular
trend and it would never fail. The two-root model with the implemented laws
(docs/REVERSIBLE_STRENGTH_PLAN.md, 2.2) says how far the fabric can take it, with
Xf = e_full tau_c/gamma_f and e_full = vL/w = 5e-13 1/s. The second batch's fabric (beta_f 10, Xf 0.5)
locks in at a share of 0.985, so its quiet root keeps 1.5% of the creep, and its 30-yr memory fades
within each earthquake cycle. Steady states (n = 3; equal sharing is unstable when N_f > 1):

  beta_f   Xf   N_f    quiet root's share when locked
  10       0.5  1.60   1.5%
  10       2    1.25   0.26%
  10       5    0.75   0.15%, but equal sharing is stable too (the separatrix is at a share of 0.89)
  10       50   0.10   none: the quiet root's own creep keeps it weak (a 3000-yr memory at gamma_f 9.45e-4)
  20       2    1.36   0.04%

So Xf near 2 deepens the lock six-fold and keeps the instability, and gamma_f lets the memory tau_c
outlast a cycle at any Xf. Runs, on the second batch's geometry:
  k1_fab20_tau300_x2    both, with a stronger fabric
  k2_fab10_tau300_x2    both
  k3_fab10_tau300_x05   the memory alone (the second batch's steady state)
  k4_fab10_tau30_x2     the depth alone (the second batch's memory)
k1 and k2 started on 2026-10-24. The inputs of k3 and k4 were moved to OUTDIR_b (their outputs stay
in OUTDIR), where tools/fill_slots.py starts them as job slots free.

Usage:
  python3 examples/two_faults/stage5_lockin.py OUTDIR
  SCYCLE_BIN=<optimized build> tools/batch.sh start OUTDIR -j 2
  python3 tools/two_fault.py OUTDIR/k1_fab20_tau300_x2/ --zref 10 --window 1000
  python3 tools/trading_stats.py OUTDIR/k*/
"""
import os, subprocess, sys

out = sys.argv[1] if len(sys.argv) > 1 else 'data/stage5_lockin'
here = os.path.dirname(os.path.abspath(__file__))
common = ['--rheology', 'powerlaw', '--lock-depth', '20', '--weak-band', '2', '--weak-factor', '1000', '--weak-bottom', '30',
          '--series-depth', '25', '--series-width', '1', '--maxTime', '3.8e11', '--out', out]
EFULL = 5e-13   # vL/w: strain rate of a 2-km band that takes the whole plate rate (1/s)
YR = 3.15e7     # the year of stage5_batch2.py's fabric_tauC0 = 9.45e8 (30 yr)


def FAB(beta, tau_yr, Xf):
    tau = tau_yr * YR
    g = EFULL * tau / Xf
    return ['fabric_type=transient', 'fabric_gammaFVals=[%.6g %.6g]' % (g, g), 'fabric_gammaFDepths=[0 500]',
            'fabric_betaF=%g' % beta, 'fabric_tauC0=%.6g' % tau]


cases = [
    ('k1_fab20_tau300_x2', FAB(20, 300, 2)),
    ('k2_fab10_tau300_x2', FAB(10, 300, 2)),
    ('k3_fab10_tau300_x05', FAB(10, 300, 0.5)),
    ('k4_fab10_tau30_x2', FAB(10, 30, 2)),
]
for name, keys in cases:
    args = [sys.executable, os.path.join(here, 'make_inputs.py')] + common + ['--name', name]
    for k in keys:
        args += ['--set', k]
    subprocess.run(args, check=True)
    p = os.path.join(out, name + '.in')
    s = open(p).read().replace('strideChkpt = 500', 'strideChkpt = 2000').replace('stride2D = 0', 'stride2D = 20000')
    open(p, 'w').write(s)
