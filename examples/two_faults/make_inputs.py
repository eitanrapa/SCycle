#!/usr/bin/env python3
"""Two parallel strike-slip faults in a full domain loaded by remote plate motion (numpy only).

Writes an SCycle input with two interior faults a distance d apart, centred in y in [0, Ly], with
the left side moving at -vL/2 and the right at +vL/2. The first fault, named fault, reads the plain
keys; the second, named f2, inherits them and overrides some with f2_ keys (see CLAUDE.md).

  --f2 rs      f2 has the same rate-and-state friction as fault (examples/ex2.in): two seismic faults
  --f2 vs      f2 is velocity strengthening everywhere (b = a - 0.004): it creeps
  --f2 locked  f2 never slips (lockedVals = 1, constantState): the single-fault limit
  --single     no f2 at all, on the same grid (the reference for --f2 locked)

Both faults start at steady sliding at vL. A smooth patch of fault (or f2, --trigger) near 8 km
depth starts at 1e-6 m/s, so the first event nucleates within hours, the same way in every run;
without it the first event grows out of round-off and its timing is arbitrary.

The grid is uniform with spacing --h0 within --band of each fault, each fault midway between two
rows, and grows geometrically (ratio --growth) toward the boundaries and between the faults.

Run from the repository root:
  python3 examples/two_faults/make_inputs.py [--d 20] [--f2 rs]
  ./source/main data/two_faults/rs_d20.in
  python3 tools/two_fault.py data/two_faults/rs_d20/
"""
import argparse, math, os
import numpy as np


def geometric_cells(D, h0, growth):
    """Cell sizes h0 r, h0 r^2, ... h0 r^n that add up to D exactly, with r close to growth."""
    if D <= 1.5*h0:
        return np.array([D])
    n = max(1, math.ceil(math.log(1 + D*(growth - 1)/(h0*growth))/math.log(growth)))
    while n > 1 and n*h0 > D:
        n -= 1
    total = lambda r: h0*n if abs(r - 1) < 1e-14 else h0*r*(r**n - 1)/(r - 1)
    lo, hi = 1.0, growth
    while total(hi) < D:
        hi *= 1.5
    for _ in range(200):
        mid = 0.5*(lo + hi)
        lo, hi = (mid, hi) if total(mid) < D else (lo, mid)
    r = 0.5*(lo + hi)
    cells = h0*r**np.arange(1, n + 1)
    return cells*(D/cells.sum())


def two_fault_grid(Ly, faults, h0, band, growth):
    """Grid rows in [0, Ly]: spacing h0 within band of each fault (each midway between two rows),
    geometric growth elsewhere. Faults closer than 2*band share one uniform band, so their
    separation must then be a multiple of h0 (h0 is adjusted down to make it one)."""
    ys = sorted(faults)
    m = max(1, int(round(band/h0)))
    # uniform bands: one per fault, or one shared by faults whose bands overlap
    groups = [[ys[0]]]
    for y in ys[1:]:
        if y - groups[-1][-1] < 2*m*h0:
            groups[-1].append(y)
        else:
            groups.append([y])
    for g in groups:  # (for two faults at most one shared band, so one h0 fits)
        if len(g) > 1:
            h0 = (g[-1] - g[0])/math.ceil((g[-1] - g[0])/h0 - 1e-9)
    rows = []
    edges = []
    for g in groups:
        lo, hi = g[0] - (m - 0.5)*h0, g[-1] + (m - 0.5)*h0
        n = int(round((hi - lo)/h0))
        rows += list(lo + h0*np.arange(n + 1))
        edges.append((lo, hi))
    if edges[0][0] <= 0 or edges[-1][1] >= Ly:
        raise SystemExit('the faults and their bands must lie inside (0, %g) km' % Ly)
    # growth toward the left and right boundaries
    c = geometric_cells(edges[0][0], h0, growth)
    rows += list(edges[0][0] - np.cumsum(c))
    c = geometric_cells(Ly - edges[-1][1], h0, growth)
    rows += list(edges[-1][1] + np.cumsum(c))
    # between bands: growth from both sides to the midpoint
    for (a0, a1), (b0, b1) in zip(edges[:-1], edges[1:]):
        c = geometric_cells(0.5*(b0 - a1), h0, growth)
        rows += list(a1 + np.cumsum(c)) + list(b0 - np.cumsum(c))
    y = np.unique(np.round(np.array(rows), 12))
    y[0], y[-1] = 0.0, Ly
    return y, h0


def write_vec(path, v):
    """PETSc binary Vec: classid 1211214, length, big-endian doubles."""
    with open(path, 'wb') as fh:
        np.array([1211214, v.size], dtype='>i4').tofile(fh)
        np.asarray(v, dtype='>f8').tofile(fh)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--d', type=float, default=20.0, help='fault separation (km)')
    p.add_argument('--Ly', type=float, default=200.0, help='domain width (km)')
    p.add_argument('--Nz', type=int, default=153, help='grid points in depth (Lz = 30 km)')
    p.add_argument('--h0', type=float, default=0.05, help='grid spacing near the faults (km)')
    p.add_argument('--band', type=float, default=4.0, help='half-width of the uniform band around each fault (km)')
    p.add_argument('--growth', type=float, default=1.05, help='spacing ratio of neighbouring cells outside the bands')
    p.add_argument('--f2', choices=['rs', 'vs', 'locked'], default='rs', help='friction of the second fault')
    p.add_argument('--sN2', type=float, default=50.0, help='effective normal stress of f2 (MPa)')
    p.add_argument('--single', action='store_true', help='only the first fault, on the same grid')
    p.add_argument('--trigger', choices=['fault', 'f2', 'none'], default='fault', help='fault with the nucleation patch')
    p.add_argument('--maxTime', type=float, default=3.0e10, help='final time (s); 3e10 s is about 950 years')
    p.add_argument('--stride1D', type=int, default=10, help='steps between 1D outputs (strideSeries is 1)')
    p.add_argument('--name', default=None, help='run name (default from the options)')
    p.add_argument('--out', default='data/two_faults', help='directory for inputs, initial conditions and output')
    args = p.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    name = args.name or ('single' if args.single else '%s_d%g' % (args.f2, args.d))
    Ly, Lz, Nz = args.Ly, 30.0, args.Nz
    y1, y2 = 0.5*(Ly - args.d), 0.5*(Ly + args.d)
    y, h0 = two_fault_grid(Ly, [y1, y2], args.h0, args.band, args.growth)

    # friction and steady state at vL (examples/ex2.in, Fault::computePsiSS)
    z = np.linspace(0, Lz, Nz)
    aD, aV = [0, 11.0, 19.3333, 60], [0.0135, 0.0300, 0.0700, 0.2652]
    bD, bV = [0, 10, 11, 19.3333, 60], [0.0230, 0.0379, 0.0350, 0.0375, 0.0497]
    f0, v0, vL, V1 = 0.6, 1e-6, 1e-9, 1e-6
    eta = 0.5*np.sqrt(30.0*3.0)                 # radiation damping sqrt(mu*rho)/2 of an interior fault
    def steady(a, b, sN):
        psi = f0 - b*np.log(vL/v0)
        return psi, a*sN*np.arcsinh(vL/(2*v0)*np.exp(psi/a)) + eta*vL
    a = np.interp(z, aD, aV)
    b1 = np.interp(z, bD, bV)
    b2 = np.interp(z, aD, np.array(aV) - 0.004) if args.f2 == 'vs' else b1
    psi1, tau1 = steady(a, b1, 50.0)
    psi2, tau2 = steady(a, b2, args.sN2)
    w = np.exp(-((z - 8.0)/1.5)**2)             # trigger patch: V = vL (V1/vL)^w at the same stress
    if args.trigger == 'fault': psi1 = psi1 - a*w*np.log(V1/vL)
    if args.trigger == 'f2' and not args.single: psi2 = psi2 - a*w*np.log(V1/vL)

    d = os.path.join(args.out, name)
    os.makedirs(os.path.join(d, 'ic'), exist_ok=True)
    write_vec(os.path.join(d, 'ic', 'y'), np.repeat(y, Nz))  # body field, index iy*Nz + iz
    write_vec(os.path.join(d, 'ic', 'psi'), psi1)
    write_vec(os.path.join(d, 'ic', 'prestress'), tau1)
    base = [l for l in open(os.path.join(root, 'examples/ex2.in')).read().splitlines() if not l.startswith((
        'outputDir', 'inputDir', 'guessSteadyStateICs', 'computeSSMomBal', 'maxTime', 'stride1D', 'stride2D',
        'restartFromChkpt', 'momBal_bc', 'Ny =', 'Nz =', 'Ly =', 'bCoordTrans', 'maxStepCount', 'strideChkpt'))]
    lines = base + [
        'Ny = %d' % y.size, 'Nz = %d' % Nz, 'Ly = %g' % Ly,
        'bCoordTrans = -1 # the grid comes from the file y in inputDir',
        'outputDir = %s/' % d, 'inputDir = %s/ic/' % d,
        'guessSteadyStateICs = 0 # interior faults take their initial state from files or keys',
        'computeSSMomBal = 0', 'maxTime = %g' % args.maxTime, 'maxStepCount = 1e8',
        'stride1D = %d' % args.stride1D, 'stride2D = 0', 'strideSeries = 1', 'strideChkpt = 500',
        'restartFromChkpt = 1 # continue from data in outputDir if a run stopped',
        'momBal_bcL_qd = remoteLoading # full domain: the left side moves at -vL/2',
        'momBal_bcR_qd = remoteLoading', 'momBal_bcT_qd = freeSurface', 'momBal_bcB_qd = freeSurface']
    if args.single:
        lines += ['interiorFaults = [fault]', 'fault_y = %.12g # (km)' % y1]
    else:
        write_vec(os.path.join(d, 'ic', 'f2_psi'), psi2)
        write_vec(os.path.join(d, 'ic', 'f2_prestress'), tau2)
        lines += ['interiorFaults = [fault f2]', 'fault_y = %.12g # (km)' % y1, 'f2_y = %.12g # (km)' % y2]
        if args.f2 == 'vs':
            lines += ['f2_bVals = [%s] # velocity strengthening: a - b = 0.004' % ' '.join('%.4f' % v for v in np.array(aV) - 0.004),
                      'f2_bDepths = [%s]' % ' '.join('%g' % v for v in aD)]
        if args.f2 == 'locked':
            lines += ['f2_lockedVals = [1 1] # f2 never slips', 'f2_lockedDepths = [0 60]',
                      'f2_stateLaw = constantState # and its state does not evolve']
        if args.sN2 != 50.0:
            lines += ['f2_sNVals = [%g %g] # (MPa)' % (args.sN2, args.sN2), 'f2_sNDepths = [0 60]']
    open(os.path.join(args.out, name + '.in'), 'w').write('\n'.join(lines) + '\n')
    dy = np.diff(y)
    print('wrote %s.in: Ny = %d (spacing %.3g km at the faults, %.3g km at most), faults at y = %g and %g km'
          % (os.path.join(args.out, name), y.size, h0, dy.max(), y1, y2))


if __name__ == '__main__':
    main()
