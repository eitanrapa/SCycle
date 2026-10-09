#!/usr/bin/env python3
"""Inputs that replicate Allison & Dunham (2018, Tectonophysics 733, 232-256) with this code.

The paper published no inputs. This reconstruction takes
- from the paper: the domain (Ly = 120 km, Lz = 40 km, a half-space), the mapping
  y = Ly sinh(10 q)/sinh(10) (appendix A), the depth spacing of its Fig. 21, dc = 16 mm, f0 = 0.6,
  V0 = 1e-6 m/s, mu = 30 GPa, rho = 3 g/cm^3, the flow laws of its Table 2 (quartz-diorite crust, wet
  olivine mantle) mixed with the law of Ji et al. (2003) over the volume fractions of its Fig. 5,
  and linear geotherms of 20, 25 and 30 K/km;
- from upstream's input max_cremeBrulee_cycle.in of June 2017 (commit dd0c35c), which matches the
  paper in domain, dc, friction and normal stress: a and b (Blanpied et al. 1991 and 1995, set with
  a 30 K/km geotherm), the effective normal stress [5.668 5.668 1336] MPa at [0 0.25 60] km
  (22.3 MPa/km), 4th-order operators, RK32 with proportional control and a tolerance of 1e-8 on slip
  and state, minDeltaT 1e-3 s, maxDeltaT 2e6 s, and Ny = 201;
- inferred: a surface temperature of 293 K (the June 2017 input; the paper's Fig. 4 and 6b agree),
  the volume fraction of olivine 0.5 [1 + tanh((z - 30 km)/3 km)] (a fit to Fig. 5 within 0.005),
  and Nz from the spacing of Fig. 21.
The paper spun up with a four-step procedure. These runs start from the steady state of the code's
guessSteadyStateICs (with computeSSMomBal = 1) and run long enough to reach the limit cycle.

Usage: make_ad2018.py --out <batch dir> [--prefix r4] [--years 15000] [--cases le g20 g25 g30]
Writes <out>/<prefix><i>_ad18_<case>.in and its ic/ folder (grid, temperature, flow law). Paths
are relative to the repository root, where tools/batch.sh runs the inputs.
"""
import argparse
import os

import numpy as np
from scipy.interpolate import PchipInterpolator

YEAR = 3.15576e7
LY, LZ, NY = 120.0, 40.0, 201
R = 8.314                              # J/(mol K), as in upstream's inputs (Q/R = 2.6341e4 K)
QD = dict(A=1.3e-3, n=2.4, Q=219e3)    # quartz-diorite, Table 2
OL = dict(A=3.6e5, n=3.5, Q=480e3)     # wet olivine, Table 2
# Fig. 21a: depth spacing (km) against depth (km), digitized
DZ_Z = [0, 4, 6, 8, 10, 12, 18, 20, 22, 24, 26, 28, 30, 34, 40]
DZ = [0.130, 0.116, 0.099, 0.073, 0.044, 0.030, 0.030, 0.092, 0.170, 0.222, 0.255, 0.273, 0.284, 0.295, 0.299]


def write_vec(path, v):
    """PETSc binary Vec: classid 1211214, length, big-endian doubles."""
    with open(path, 'wb') as fh:
        np.array([1211214, v.size], dtype='>i4').tofile(fh)
        np.asarray(v, dtype='>f8').tofile(fh)


def depth_grid():
    """Nodes whose spacing follows Fig. 21: march with the midpoint spacing, then stretch to Lz."""
    dz = PchipInterpolator(DZ_Z, DZ)
    z = [0.0]
    while z[-1] < LZ:
        h = float(dz(min(z[-1] + 0.5 * float(dz(z[-1])), LZ)))
        z.append(z[-1] + h)
    z = np.array(z)
    if z[-1] - LZ > 0.5 * (z[-1] - z[-2]):   # nearer the node before: drop the overshoot
        z = z[:-1]
    return z * LZ / z[-1]


def mixed_flow_law(z):
    """Ji et al. (2003) for layered phases: n, A (MPa^-n s^-1) and Q/R (K) at each depth."""
    p2 = 0.5 * (1.0 + np.tanh((z - 30.0) / 3.0))
    p1 = 1.0 - p2
    n = 1.0 / (p1 / QD['n'] + p2 / OL['n'])
    A = QD['A'] ** (n * p1 / QD['n']) * OL['A'] ** (n * p2 / OL['n'])
    Q = n * (p1 * QD['Q'] / QD['n'] + p2 * OL['Q'] / OL['n'])
    return n, A, Q / R


TEMPLATE = """# Allison & Dunham (2018), {title}: written by examples/replication/make_ad2018.py
inputDir = {ic}/
outputDir = {out}/
order = 4
Ny = {Ny}
Nz = {Nz}
Ly = {Ly:g}
Lz = {Lz:g}
bCoordTrans = -1 # y and z come from the files y and z in inputDir
vL = 1e-9
momentumBalanceType = quasidynamic
bulkDeformationType = {bulk}
guessSteadyStateICs = 1
computeSSMomBal = 1
linSolverSS = MUMPSCHOLESKY
linSolverTrans = MUMPSCHOLESKY
muVals = [30 30]
muDepths = [0 500]
rhoVals = [3 3]
rhoDepths = [0 500]
{rheology}
stateLaw = agingLaw
f0 = 0.6
v0 = 1e-6
DcVals = [16e-3 16e-3]
DcDepths = [0 60]
aVals = [0.0135 0.0300 0.0700 0.2652]
aDepths = [0 11.0000 19.3333 60]
bVals = [0.0230 0.0379 0.0350 0.0375 0.0497]
bDepths = [0 10 11 19.3333 60]
sNVals = [5.668 5.668 1336]
sNDepths = [0 0.25 60]
thermalCoupling = no
heatEquationType = steadyState
withViscShearHeating = no
withFrictionalHeating = no
timeIntegrator = RK32
timeControlType = P
normType = L2_absolute
timeStepTol = 1e-8
timeIntInds = [slip psi]
initDeltaT = 1e-3
minDeltaT = 1e-3
maxDeltaT = 2e6
maxStepCount = 1e8
maxTime = {maxTime:.6g}
stride1D = 10
stride2D = 5000
strideChkpt = 2000
restartFromChkpt = 1
"""

POWERLAW = """# power law: dislocation creep, with A, Q/R and n of the mixed crust and mantle from the files
# disl_A, disl_QR and disl_n in inputDir (the values below are placeholders the files replace);
# the temperature comes from the file T in inputDir
wDislCreep = yes
wDiffCreep = no
wLinearMaxwell = no
wPlasticity = no
disl_AVals = [1.3e-3 1.3e-3]
disl_ADepths = [0 500]
disl_QRVals = [2.6341e4 2.6341e4]
disl_QRDepths = [0 500]
disl_nVals = [2.4 2.4]
disl_nDepths = [0 500]
maxEffVisc = 1e17
grainSizeEv_grainSizeEvType = no
grainSizeEv_grainSizeEvTypeSS = no
grainSizeEvCoupling = no
grainSizeEvCouplingSS = no
# the power-law model always builds the heat equation, which needs these (no heat equation is
# solved: thermalCoupling = no; k and c from the June 2017 input); TVals is the geotherm of the
# file T, which replaces it
kVals = [1.89e-9 1.89e-9]
kDepths = [0 60]
cVals = [900 900]
cDepths = [0 60]
withRadioHeatGeneration = no
TVals = [293 {T40:g}]
TDepths = [0 40]"""


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', required=True)
    p.add_argument('--prefix', default='r4')
    p.add_argument('--years', type=float, default=15000.0)
    p.add_argument('--cases', nargs='+', default=['le', 'g20', 'g25', 'g30'])
    a = p.parse_args()

    root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    rel = lambda path: os.path.relpath(os.path.abspath(path), root)
    q = np.linspace(0.0, 1.0, NY)
    y = LY * np.sinh(10.0 * q) / np.sinh(10.0)
    z = depth_grid()
    Nz = z.size
    Y, Z = np.repeat(y, Nz), np.tile(z, NY)        # body fields, index iy*Nz + iz
    n, A, QR = mixed_flow_law(Z)
    print(f'grid {NY} x {Nz}: dy {1e3 * (y[1] - y[0]):.3f} m at the fault, dz {1e3 * np.diff(z).min():.1f} to '
          f'{1e3 * np.diff(z).max():.1f} m; mixed law at 30 km: n {np.interp(30, z, n[:Nz]):.3f}, '
          f'A {np.interp(30, z, A[:Nz]):.3g}, Q {np.interp(30, z, QR[:Nz]) * R / 1e3:.0f} kJ/mol')
    for i, case in enumerate(a.cases):
        tag = f'{a.prefix}{chr(ord("a") + i)}_ad18_{case}'
        d = os.path.join(a.out, tag)
        os.makedirs(os.path.join(d, 'ic'), exist_ok=True)
        write_vec(os.path.join(d, 'ic', 'y'), Y)
        write_vec(os.path.join(d, 'ic', 'z'), Z)
        if case == 'le':
            title, bulk, rheology = 'linear elastic', 'linearElastic', 'wDislCreep = no'
        else:
            G = float(case[1:])
            title, bulk, rheology = f'power law, {G:g} K/km', 'powerLaw', POWERLAW.format(T40=293.0 + 40.0 * G)
            write_vec(os.path.join(d, 'ic', 'T'), 293.0 + G * Z)
            write_vec(os.path.join(d, 'ic', 'disl_A'), A)
            write_vec(os.path.join(d, 'ic', 'disl_QR'), QR)
            write_vec(os.path.join(d, 'ic', 'disl_n'), n)
        with open(os.path.join(a.out, tag + '.in'), 'w') as f:
            f.write(TEMPLATE.format(title=title, ic=rel(os.path.join(d, 'ic')), out=rel(d), Ny=NY, Nz=Nz,
                                    Ly=LY, Lz=LZ, bulk=bulk, rheology=rheology, maxTime=a.years * YEAR))
        print(os.path.join(a.out, tag + '.in'))


if __name__ == '__main__':
    main()
