#!/usr/bin/env python3
"""Second batch of the stage 5 run matrix (docs/REVERSIBLE_STRENGTH_PLAN.md, sections 2.2 and 5).

The first batch (stage5_batch1.py) shared slip equally in all six runs. Screened with the laws as
implemented (tools/zero_d_trading.py --laws implemented), that is what each of them predicts: the
fabric state weakens the strength by 1/(1 + beta Phi) and saturates, so its gain at equal sharing,
N_f = n beta u/((1 + u)(1 + u + beta u)) with u = e tau_c/gamma_f, never exceeds 1 for beta <= 3
(it was 0.80, 0.65 and 0.29 there), and the hardening memory was saturated in both roots
(e_full tau_r/gamma_h = 20). This batch keeps the first one's geometry (two faults 20 km apart,
locked below 20 km, a weak band 2 km wide from 20 to 30 km beneath each, the stage 4 power law) and
takes its state parameters from the screen, with e_full = vL/(2 km) = 5e-13 /s (inside the bands
the first batch's roots strained at 1.7e-13 to 2.8e-13 /s at equal sharing, as assumed):

  memory  tau_r = 3.15e10 s (1000 yr), gamma_h = 6.3e-3: Xh = e_full tau_r/gamma_h = 2.5
  fabric  tau_c = 9.45e8 s (30 yr), gamma_f = 9.45e-4: Xf = e_full tau_c/gamma_f = 0.5, r = 0.03

            beta_f  hard_a   N_f   screen (two roots in series, implemented laws)
  fab5        5       -     1.20   lock-in, f = 0.89                  fabric alone: gate 5.5b,
  fab8        8       -     1.48   lock-in, f = 0.97                  with the first batch's fab
  fab10      10       -     1.60   lock-in, f = 0.98                  (N_f = 0.80: sharing)
  alt5        5      0.5    1.20   alternation, period 1740 yr, f 0.13-0.87
  alt8        8      0.7    1.48   alternation, period 3330 yr, f 0.03-0.97
  alt10      10      1.0    1.60   alternation, period 3360 yr, f 0.01-0.99
  lock8       8      0.5    1.48   lock-in, f = 0.89 (the memory cannot overturn the fabric)
  lock10     10      0.7    1.60   lock-in, f = 0.92

For each beta the memory amplitude straddles the screen's alternation boundary (b_max = b_c), so
the batch tests that boundary as well as the lock-in threshold. hard_a = 0.5-0.7 is a hardened to
recovered stress ratio of 1.5-1.7 at most, the range the hardening literature supports; 1.0 (alt10)
is beyond it. The predictions hold for e_full 20% lower or higher. The 2D model adds what the
screen lacks (the deep substrate beneath the bands, elastic coupling, the earthquake cycle), so its
thresholds may lie higher; gate 5.5b allows a factor of 2.

Each runs 3.8e11 s (12,000 yr), with a checkpoint every 2000 steps and a 2D snapshot every 20000;
the root probes in faultSeries.txt average over the band (seriesWidth 1 km, the first batch's 2 km
diluted them with the slow rock beside the band).
Usage:
  python3 examples/two_faults/stage5_batch2.py OUTDIR
  python3 tools/zero_d_trading.py --laws implemented --input OUTDIR/alt8.in --efull 5e-13
  tools/batch.sh start OUTDIR            (detached; tools/batch.sh status OUTDIR follows it)
  python3 tools/two_fault.py OUTDIR/alt8/ --zref 10 --window 1000
"""
import os, subprocess, sys

out = sys.argv[1] if len(sys.argv) > 1 else 'data/stage5_batch2'
here = os.path.dirname(os.path.abspath(__file__))
common = ['--rheology', 'powerlaw', '--lock-depth', '20', '--weak-band', '2', '--weak-factor', '1000', '--weak-bottom', '30',
          '--series-depth', '25', '--series-width', '1', '--maxTime', '3.8e11', '--out', out]
HARD = lambda a: ['hard_type=transient', 'hard_a=%g' % a, 'hard_gammaHVals=[6.3e-3 6.3e-3]', 'hard_gammaHDepths=[0 500]',
                  'hard_tauR0=3.15e10', 'hard_SRef=0']
FAB = lambda beta: ['fabric_type=transient', 'fabric_gammaFVals=[9.45e-4 9.45e-4]', 'fabric_gammaFDepths=[0 500]',
                    'fabric_betaF=%g' % beta, 'fabric_tauC0=9.45e8']
cases = {'fab5': FAB(5), 'fab8': FAB(8), 'fab10': FAB(10),
         'alt5': HARD(0.5) + FAB(5), 'alt8': HARD(0.7) + FAB(8), 'alt10': HARD(1.0) + FAB(10),
         'lock8': HARD(0.5) + FAB(8), 'lock10': HARD(0.7) + FAB(10)}
for name, keys in cases.items():
    args = [sys.executable, os.path.join(here, 'make_inputs.py')] + common + ['--name', name]
    for k in keys: args += ['--set', k]
    subprocess.run(args, check=True)
    p = os.path.join(out, name + '.in')
    s = open(p).read().replace('strideChkpt = 500', 'strideChkpt = 2000').replace('stride2D = 0', 'stride2D = 20000')
    open(p, 'w').write(s)
