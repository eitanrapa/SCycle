#!/usr/bin/env python3
"""First batch of the stage 5 run matrix (docs/REVERSIBLE_STRENGTH_PLAN.md, section 5).

Two faults 20 km apart (the stage 4 power-law baseline), locked below 20 km, each with its own
viscous root: a band 2 km wide from 20 to 30 km depth whose dislocation-creep prefactor is 1000 times
larger (make_inputs.py --weak-band). Measured in that model, the roots strain at 2e-14 to 2e-13 /s
between events, the crust between and beyond them at 1e-19 to 4e-15 /s; at the full plate rate a
2 km root strains at e_full = vL/(2 km) = 5e-13 /s.

The states pair a slow strengthening memory (strain hardening, 4.1) with a fast weakening feedback
(the fabric state, 4.5), mapped onto the 0D model of section 2 (tools/zero_d_trading.py, n = 3) with
the time unit t0 = gamma'/e_full = 160 yr (gamma' = gamma_h/a, the strain that adds 1 to the
strength): R = gamma'/gamma_f = 10, T_h = tau_r/t0 = 10, T_D = tau_c/t0 and beta = beta_f. These are
cell-scale (root-averaged) strains, far below laboratory values: a core of width w_c inside the root
strains faster by (2 km)/w_c, so gamma_h = 1.25e-3 here stands for 0.06 in a 40 m core.

  null    no states                                    (the null experiment)
  hard    hardening alone                              (memory alone: equal sharing expected)
  fab     fabric alone, beta 2, T_D 0.1 (N_loc = 3)    (feedback alone: lock-in expected)
  alt     both, beta 2, T_D 0.1                        (0D: alternation, period 24.7 t0 = 4000 yr)
  lockin  both, beta 2, T_D 0.3 (N_loc = 9)            (0D: lock-in)
  share   both, beta 0.5, T_D 0.1 (N_loc = 0.75)       (0D: equal sharing)

Each runs 3.8e11 s (12,000 yr), with a checkpoint every 2000 steps and a 2D snapshot every 20000.
Usage:
  python3 examples/two_faults/stage5_batch1.py OUTDIR
  for n in null hard fab alt lockin share; do ./source/main OUTDIR/$n.in > OUTDIR/$n/run.log & done
  python3 tools/two_fault.py OUTDIR/alt/ --zref 10 --window 1000
"""
import os, subprocess, sys

out = sys.argv[1] if len(sys.argv) > 1 else 'data/stage5_batch1'
here = os.path.dirname(os.path.abspath(__file__))
common = ['--rheology', 'powerlaw', '--lock-depth', '20', '--weak-band', '2', '--weak-factor', '1000', '--weak-bottom', '30',
          '--series-depth', '25', '--maxTime', '3.8e11', '--out', out]
HARD = ['hard_type=transient', 'hard_a=0.5', 'hard_gammaHVals=[1.25e-3 1.25e-3]', 'hard_gammaHDepths=[0 500]',
        'hard_tauR0=5.05e10', 'hard_SRef=0']
FAB = lambda beta, tauc: ['fabric_type=transient', 'fabric_gammaFVals=[2.5e-4 2.5e-4]', 'fabric_gammaFDepths=[0 500]',
                          'fabric_betaF=%g' % beta, 'fabric_tauC0=%g' % tauc]
cases = {'null': [], 'hard': HARD, 'fab': FAB(2, 5.05e8), 'alt': HARD + FAB(2, 5.05e8),
         'lockin': HARD + FAB(2, 1.515e9), 'share': HARD + FAB(0.5, 5.05e8)}
for name, keys in cases.items():
    args = [sys.executable, os.path.join(here, 'make_inputs.py')] + common + ['--name', name]
    for k in keys: args += ['--set', k]
    subprocess.run(args, check=True)
    p = os.path.join(out, name + '.in')
    s = open(p).read().replace('strideChkpt = 500', 'strideChkpt = 2000').replace('stride2D = 0', 'stride2D = 20000')
    open(p, 'w').write(s)
