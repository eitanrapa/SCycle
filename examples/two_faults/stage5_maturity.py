#!/usr/bin/env python3
"""Stage 5 maturity tests: does a fault whose root is permanently weaker, by inheritance or by its own
strain, end the trading?

No run of the lock-in search (stage5_lockin.py) locked in: every quiet fault kept loading at 3-8% of
the active fault's rate. During the quiet fault's last silence the active fault's band took 87-99% of
the viscous shear at every depth to the base at 58 km (k5, k1; the base is traction-free), the quiet
fault's band 0.8-1.5%, and 5-10% spread outside both below 30 km. Both faults have the same permanent
weak band (ic/disl_A 1000 x the prefactor, +-1 km, 20-30 km), so the model starts with two equally
mature faults and only the fabric decides who leads. These runs vary k2_fab10_tau300_x2 (beta_f 10,
tau_c 300 yr, Xf = e_full tau_c/gamma_f = 2, e_full = 5e-13 1/s):

  y1_fab10_tau300_x2_B1       a young second fault: no band under f2 (--weak-factor-f2 1)
  y2_fab10_tau300_x2_B100     f2's band 100 x instead of 1000 x
  y3_fab10_noheal_noband      no bands, and the fabric never heals (tau_c 1e8 yr)
  y4_fab10_tau300_x2_noband   no bands, k2's fabric as it is
  y5_fab10_tau300_x2_mbase60  the base moving with the plates across 60 km, centred between the
                              faults (--moving-base 60): the deep shear forced under both faults
  y6_fab10_noheal_k3          the bands, with the fabric of k3_fab10_tau300_x05 (gamma_f 4 x k2's)
                              but never healing; k3 itself (P(same) 0.17) is the healing twin

In y1 and y2 the young fault f2 has the nucleation patch (--trigger f2), so the mature fault gets no
head start. The fabric factor divides the creep strength: at a given stress a root creeps
(1 + beta_f Phi)^n faster, so beta_f 10 at Phi = 1 is 1331 x, the band's 1000 x; a larger cap would
also weaken the hot mantle, where the shortest Maxwell time caps the step (k1 and k2 at 45-60 km under
the active fault, Phi 0.6-0.75).

A race of two roots in series with the code's laws (shares P_i (1 + beta Phi_i)^n, P_i the band factor;
an earthquake on a root adds the strain s/w, 7 m over a 2-km band), leaving out the quiet fault's extra
loading (the 3-8% above), which is what the runs measure:
- y2: f2's band at 100 x loses to the first fault's 1000 x although f2 breaks first: its share is 97% at
  100 yr, 0.4% at 1000 yr and 0.02% from 4000 yr. With 1000 x (k2) the fault that breaks first keeps 99.7%.
- y6: one earthquake gives its root Phi 0.31 for good, and the leader passes Phi 0.9 after about 44 m of
  slip (1400 yr at the plate rate), when it takes 99.9% against a quiet root at Phi 0. The quiet fault's
  own earthquakes bring its root to 0.31, 0.52, 0.67, 0.77 and its share to 5%, 15%, 25%, 33%: the fabric
  saturates, so the quiet fault comes back in the end, and the run says how long that takes. Healed
  over 300 yr (k3) a root loses its Phi before its next turn, and the faults alternate.

Without bands there are no two roots to begin with (plan, section 6, item 10). Measured 384 yr after the
first earthquake in a trial of y3 with gamma_f 5 x k2's (Phi below 0.02): all the plate motion flows below
30 km, nearly uniformly across the domain, at 8e-15 1/s under the fault that slipped, 6e-15 between the
faults and under the other fault and 4.6e-15 at y = 160 km (vL/Ly = 5e-15), so the windows 85-95 and
105-115 km take 8% and 6%. Above that the rock barely creeps: under the fault that slipped 9e-17 1/s at
22 km (half-width 2 km) and 2e-15 at 25 km (4 km). So y3 and y4 ask whether the fabric can make roots by
itself and whether healing prevents it. At 6-8e-15 1/s k2's gamma_f has an e-folding time of about
10000 yr (5 x k2's would leave Phi below 0.25 at the end of the run). Healed over 300 yr Phi stays near
0.03 and the two-root gain N_f (plan 2.2) at that rate is about 0.7, below 1, so y4 should stay broad.
Unhealed, Phi keeps growing, and in the two-root model differences between the roots grow with an
e-folding time of a few hundred years until Phi reaches 0.7. The fabric acts at every depth as in the
lock-in runs.

Usage:
  python3 examples/two_faults/stage5_maturity.py OUTDIR [NAME ...]   (NAME: write only those cases)
  SCYCLE_BIN=<optimized build> tools/batch.sh start OUTDIR -j 6
  python3 tools/two_fault.py OUTDIR/y1_fab10_tau300_x2_B1/ --zref 10 --window 1000
  python3 tools/trading_stats.py OUTDIR/y*/
"""
import os, subprocess, sys

out = sys.argv[1] if len(sys.argv) > 1 else 'data/stage5_maturity'
only = sys.argv[2:]
here = os.path.dirname(os.path.abspath(__file__))
common = ['--rheology', 'powerlaw', '--lock-depth', '20', '--series-depth', '25', '--series-width', '1',
          '--maxTime', '3.8e11', '--out', out]
BAND = ['--weak-band', '2', '--weak-factor', '1000', '--weak-bottom', '30']
EFULL = 5e-13   # vL/w: strain rate of a 2-km band that takes the whole plate rate (1/s)
YR = 3.15e7     # the year of stage5_lockin.py
GK2 = EFULL * 300 * YR / 2   # k2's gamma_f (Xf 2 at tau_c 300 yr)


def FAB(beta, gamma, tau_yr):
    return ['fabric_type=transient', 'fabric_gammaFVals=[%.6g %.6g]' % (gamma, gamma), 'fabric_gammaFDepths=[0 500]',
            'fabric_betaF=%g' % beta, 'fabric_tauC0=%.6g' % (tau_yr * YR)]


cases = [
    ('y1_fab10_tau300_x2_B1', BAND + ['--weak-factor-f2', '1', '--trigger', 'f2'], FAB(10, GK2, 300)),
    ('y2_fab10_tau300_x2_B100', BAND + ['--weak-factor-f2', '100', '--trigger', 'f2'], FAB(10, GK2, 300)),
    ('y3_fab10_noheal_noband', [], FAB(10, GK2, 1e8)),
    ('y4_fab10_tau300_x2_noband', [], FAB(10, GK2, 300)),
    ('y5_fab10_tau300_x2_mbase60', BAND + ['--moving-base', '60'], FAB(10, GK2, 300)),
    ('y6_fab10_noheal_k3', BAND, FAB(10, 4 * GK2, 1e8)),
]
unknown = set(only) - {name for name, _, _ in cases}
if unknown:
    sys.exit('no such case: ' + ' '.join(sorted(unknown)))
for name, opts, keys in cases:
    if only and name not in only:
        continue
    args = [sys.executable, os.path.join(here, 'make_inputs.py')] + common + opts + ['--name', name]
    for k in keys:
        args += ['--set', k]
    subprocess.run(args, check=True)
    p = os.path.join(out, name + '.in')
    s = open(p).read().replace('strideChkpt = 500', 'strideChkpt = 2000').replace('stride2D = 0', 'stride2D = 20000')
    open(p, 'w').write(s)
