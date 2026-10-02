#!/usr/bin/env python3
"""Create the grid and initial-condition files for SEAS benchmark BP1 (numpy only).

Python version of createICs.m, which needs constructCoord_constL1 (never committed)
and wrote to a hard-coded path. Writes PETSc binary Vecs y, z (size Ny*Nz, z fastest)
and psi, prestress, tau, tauQS (size Nz, along the fault) into ICs/ next to this file,
and prints the Ny and Nz that BP1.in must use.

Grid: uniform spacing dz (25 m in the benchmark) to Hy = 5 km from the fault and to
Hz = 55 km depth, then 200 cells growing geometrically to Ly = Lz = 500 km.
The growth law of the original constructCoord_constL1 is unknown; any smooth
stretching outside the seismogenic zone serves the benchmark.

Usage: SEAS_benchmarks/BP1/createICs.py [--dz 0.025] [--out DIR]
"""
import argparse, os
import numpy as np

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--dz', type=float, default=0.025, help='uniform grid spacing near the fault (km), BP1: 0.025')
p.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ICs'))
args = p.parse_args()

# BP1 parameters (Erickson et al. 2020)
Ly = Lz = 500.0       # km
Hy, Hz = 5.0, 55.0    # km, uniform spacing to these distances
N2 = 201              # points in each stretched part
rho, cs = 2.670, 3.464  # g/cm^3, km/s
mu = rho * cs**2      # GPa (32.04)
eta = mu * 1e3 / (2 * cs * 1e3)  # radiation damping mu/(2 cs), MPa s/m (4.62)
sN = 50.0             # MPa
H, h = 15.0, 3.0      # km
v0, f0, Dc = 1e-6, 0.6, 0.008
Vp = Vinit = 1e-9     # m/s
a0, amax, b0 = 0.010, 0.025, 0.015


def stretched_coord(L, Hc, d1, n2):
    """Uniform spacing d1 on [0, Hc], then n2-1 cells growing by a constant ratio to L."""
    n1 = int(round(Hc / d1)) + 1
    x1 = np.linspace(0.0, Hc, n1)
    m, target = n2 - 1, (L - Hc) / d1
    lo, hi = 1.0 + 1e-12, 2.0  # solve sum_{k=1..m} r^k = target for r
    for _ in range(200):
        r = 0.5 * (lo + hi)
        s = r * (r**m - 1.0) / (r - 1.0)
        lo, hi = (r, hi) if s < target else (lo, r)
    x2 = Hc + np.cumsum(d1 * r ** np.arange(1, m + 1))
    x2[-1] = L
    return np.concatenate([x1, x2]), r


y, ry = stretched_coord(Ly, Hy, args.dz, N2)
z, rz = stretched_coord(Lz, Hz, args.dz, N2)
Ny, Nz = y.size, z.size

a = np.interp(z, [0, H, H + h, Lz], [a0, a0, amax, amax])
f = amax * np.arcsinh((0.5 * Vinit / v0) * np.exp((f0 + b0 * np.log(v0 / Vinit)) / amax))
tau0 = sN * f + eta * Vinit + 0 * z   # quasi-static shear stress (prestress)
tau = sN * f + 0 * z                  # frictional strength
theta = (Dc / v0) * np.exp((a / b0) * np.log((2 * v0 / Vinit) * np.sinh((tau0 - eta * Vinit) / (a * sN))) - f0 / b0)
psi = f0 + b0 * np.log(theta * v0 / Dc)


def write_vec(path, v):
    """PETSc binary Vec: classid 1211214, length, then big-endian doubles."""
    with open(path, 'wb') as fh:
        np.array([1211214, v.size], dtype='>i4').tofile(fh)
        np.asarray(v, dtype='>f8').tofile(fh)


os.makedirs(args.out, exist_ok=True)
Y = np.repeat(y, Nz)   # index iy*Nz + iz: y constant within each column of Nz values
Z = np.tile(z, Ny)
for name, v in [('y', Y), ('z', Z), ('psi', psi), ('prestress', tau0), ('tau', tau), ('tauQS', tau0)]:
    write_vec(os.path.join(args.out, name), v)

print('wrote y, z, psi, prestress, tau, tauQS to %s' % args.out)
print('Ny = %d\nNz = %d' % (Ny, Nz))
print('(growth ratio of the stretched cells: y %.4f, z %.4f; tau0 = %.4f MPa, psi(0) = %.6f)' % (ry, rz, tau0[0], psi[0]))
