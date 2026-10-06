#!/usr/bin/env python3
"""Stage 5 run matrix, phase 1 (docs/REVERSIBLE_STRENGTH_PLAN.md, section 5): every mechanism that acts on
the ductile roots, alone and together, at two strengths where it matters.

The geometry is the second batch's (stage5_batch2.py): two faults 20 km apart, locked below 20 km, a weak
band 2 km wide from 20 to 30 km beneath each (1000 times the dislocation-creep prefactor), the stage 4
power law with coupled heat (W18, W19 on in every run), 12,000 yr. The second batch gives the references:
fab8, fab10 (fabric alone), alt8, alt10, lock8, lock10 (fabric with hardening).

Memories are scaled as the second batch's hardening, e_full tau/gamma = 2.5 with e_full = vL/(2 km) =
5e-13 /s and tau = 1000 yr (3.15e10 s):
  HARD(a)    strain hardening (S7) with recovery (W14), gamma_h 6.3e-3; a = 0.7 as alt8
  WATER      drying with strain (S1, S2) against recharge (W2, W3) over tau_hyd, gamma_dry 6.3e-3,
             chi_min 0.3; with --band-wet-dry R the band is the wet end-member, the dry one R times
             slower (R = 10, 30: the plan's 3-30)
  SEG(b)     phase segregation (S10) slowing pressure solution 1 + b Xi; the band's pressure solution
             (--band-dp 0.5: half the band's weakening at 2 MPa) strains at half e_full, so gamma_s is
             3.15e-3 for the same 2.5
Feedbacks:
  FAB(beta)        fabric (W11, W13), tau_c 30 yr, gamma_f 9.45e-4 (Xf 0.5), as the second batch
  FAB(beta, True)  the same, directional (W12): fault-parallel shear only

  m01_aniso8       W12 alone (beta 8): directional against fab8
  m02_wet10        hydration against dehydration (the paper's first pair) alone, R = 10
  m03_dpnull       the pressure-solution band without states (reference for m04, m06, m07)
  m04_dpseg4       S10 alone, beta_s 4
  m05_fab8wet10    fabric with the water memory
  m06_dpfab8       fabric in the pressure-solution band (reference for m07)
  m07_dpfab8seg4   fabric with the segregation memory
  m08_aniso8hard7  directional fabric with hardening (against alt8)
  m09_all          directional fabric, hardening, water and segregation in a wet-dry pressure-solution band
  m10_aniso10      W12 alone, beta 10 (against fab10)
  m11_wet30        the water memory alone, R = 30
  m12_fab10wet30   fabric with the water memory at the larger strengths
  m13_fab8hard7wet10  fabric with two memories

Not in phase 1, since they cannot reach the roots in this geometry: the cement state (S6, S9), the
pseudotachylite products (W7) and the cataclastic sink (W6) take the faults' work, which is done above
20 km where the faults slip, not in the bands below; pore-pressure pulses (W4) and cohesion healing act on
the faults' friction; grain size (W8, S3, W10) needs a band whose weakness is its fine grain. Phase 2.

Usage:
  python3 examples/two_faults/stage5_matrix.py OUTDIR
  SCYCLE_BIN=<optimized build> tools/batch.sh start OUTDIR -j 2
  python3 tools/two_fault.py OUTDIR/m02_wet10/ --zref 10 --window 1000
"""
import os, subprocess, sys

out = sys.argv[1] if len(sys.argv) > 1 else 'data/stage5_matrix'
here = os.path.dirname(os.path.abspath(__file__))
common = ['--rheology', 'powerlaw', '--lock-depth', '20', '--weak-band', '2', '--weak-factor', '1000', '--weak-bottom', '30',
          '--series-depth', '25', '--series-width', '1', '--maxTime', '3.8e11', '--out', out]
HARD = lambda a: ['hard_type=transient', 'hard_a=%g' % a, 'hard_gammaHVals=[6.3e-3 6.3e-3]', 'hard_gammaHDepths=[0 500]',
                  'hard_tauR0=3.15e10', 'hard_SRef=0']
FAB = lambda beta, aniso=False: (['fabric_type=transient', 'fabric_gammaFVals=[9.45e-4 9.45e-4]', 'fabric_gammaFDepths=[0 500]',
                                  'fabric_betaF=%g' % beta, 'fabric_tauC0=9.45e8'] + (['fabric_anisotropic=1'] if aniso else []))
WATER = ['water_type=transient', 'water_chiSupVals=[1 1]', 'water_chiSupDepths=[0 500]', 'water_chiMin=0.3',
         'water_tauHyd=3.15e10', 'water_gammaDry=6.3e-3', 'water_mix=log']
SEG = lambda b: ['seg_type=transient', 'seg_betaS=%g' % b, 'seg_gammaSVals=[3.15e-3 3.15e-3]', 'seg_gammaSDepths=[0 500]',
                 'seg_tauM0=3.15e10']
WD = lambda R: ['--band-wet-dry', str(R)]
DP = ['--band-dp', '0.5']
cases = [
    ('m01_aniso8', [], FAB(8, True)),
    ('m02_wet10', WD(10), WATER),
    ('m03_dpnull', DP, []),
    ('m04_dpseg4', DP, SEG(4)),
    ('m05_fab8wet10', WD(10), WATER + FAB(8)),
    ('m06_dpfab8', DP, FAB(8)),
    ('m07_dpfab8seg4', DP, FAB(8) + SEG(4)),
    ('m08_aniso8hard7', [], FAB(8, True) + HARD(0.7)),
    ('m09_all', DP + WD(10), FAB(8, True) + HARD(0.7) + WATER + SEG(4)),
    ('m10_aniso10', [], FAB(10, True)),
    ('m11_wet30', WD(30), WATER),
    ('m12_fab10wet30', WD(30), WATER + FAB(10)),
    ('m13_fab8hard7wet10', WD(10), WATER + FAB(8) + HARD(0.7)),
]
for name, extra, keys in cases:
    args = [sys.executable, os.path.join(here, 'make_inputs.py')] + common + extra + ['--name', name]
    for k in keys: args += ['--set', k]
    subprocess.run(args, check=True)
    p = os.path.join(out, name + '.in')
    s = open(p).read().replace('strideChkpt = 500', 'strideChkpt = 2000').replace('stride2D = 0', 'stride2D = 20000')
    open(p, 'w').write(s)
