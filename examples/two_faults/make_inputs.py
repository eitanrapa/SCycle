#!/usr/bin/env python3
"""Two parallel strike-slip faults in a full domain loaded by remote plate motion (numpy only).

Writes an SCycle input with two interior faults a distance d apart, centred in y in [0, Ly], with
the left side moving at -vL/2 and the right at +vL/2. The first fault, named fault, reads the plain
keys; the second, named f2, inherits them and overrides some with f2_ keys (see CLAUDE.md).

  --f2 rs      f2 has the same rate-and-state friction as fault (examples/ex2.in): two seismic faults
  --f2 vs      f2 is velocity strengthening everywhere (b = a - 0.004): it creeps
  --f2 locked  f2 never slips (lockedVals = 1, constantState): the single-fault limit
  --single     no f2 at all, on the same grid (the reference for --f2 locked)

  --rheology powerlaw  the stage 4 baseline: ex4's dislocation creep and geotherm, the domain
               deepened to --Lz (60 km) on a depth grid that is uniform (--dz) down to --zfine and
               grows geometrically below, temperature evolving with frictional heat (a Gaussian of
               width --w metres on each fault) and viscous heating, and the implicit-explicit
               integrator. --grainsize adds the wattmeter grain-size law with the quartz constants
               of docs/REVERSIBLE_STRENGTH_PLAN.md (section 4.3), uncoupled (no grain-size-sensitive
               creep law yet).

Stage 5 (docs/REVERSIBLE_STRENGTH_PLAN.md): --lock-depth Z locks both faults below Z km, so that
the plate motion beneath them flows viscously (without it the deep fault creeps and the viscous roots
barely strain); --series-depth Z writes the root-strength probes at Z km to faultSeries.txt; and
--set KEY=VALUE (repeatable) appends any input line, for the bulk state fields and the other stage 5
keys. --weak-band W gives each fault its own viscous root: a band W km wide beneath it, below
--weak-top km and above --weak-bottom km, where the dislocation-creep prefactor is --weak-factor
times larger (the file
ic/disl_A, which the power law reads as a body field; log-linear cosine tapers of --weak-taper km at
the edges). Without one the flow beneath the faults is broad and the two roots merge.
--mantle-zone W adds a mantle shear zone W km wide beneath both bands (centred between the faults),
weak like them from --weak-bottom to --mantle-bottom and tapering into the asthenosphere below, so
that the plate motion beneath the roots stays localized instead of spreading through the mantle.
--moving-base W moves the bottom boundary with the plates instead of leaving it traction-free:
from -vL/2 to +vL/2 across a transition W km wide centred between the faults (--base-center), so
the base carries load beneath the far field (the faults must then be locked at depth).

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
    p.add_argument('--Nz', type=int, default=153, help='grid points in depth (elastic, Lz = 30 km)')
    p.add_argument('--rheology', choices=['elastic', 'powerlaw'], default='elastic', help='bulk rheology')
    p.add_argument('--Lz', type=float, default=60.0, help='depth of the power-law domain (km)')
    p.add_argument('--dz', type=float, default=0.2, help='power law: depth spacing down to --zfine (km)')
    p.add_argument('--zfine', type=float, default=25.0, help='power law: depth of the uniform part of the depth grid (km)')
    p.add_argument('--w', type=float, default=10.0, help='power law: width of the frictional heat source (m)')
    p.add_argument('--grainsize', action='store_true', help='power law: evolve the grain size (uncoupled)')
    p.add_argument('--h0', type=float, default=0.05, help='grid spacing near the faults (km)')
    p.add_argument('--band', type=float, default=4.0, help='half-width of the uniform band around each fault (km)')
    p.add_argument('--growth', type=float, default=1.05, help='spacing ratio of neighbouring cells outside the bands')
    p.add_argument('--f2', choices=['rs', 'vs', 'locked'], default='rs', help='friction of the second fault')
    p.add_argument('--sN2', type=float, default=50.0, help='effective normal stress of f2 (MPa)')
    p.add_argument('--single', action='store_true', help='only the first fault, on the same grid')
    p.add_argument('--trigger', choices=['fault', 'f2', 'none'], default='fault', help='fault with the nucleation patch')
    p.add_argument('--maxTime', type=float, default=3.0e10, help='final time (s); 3e10 s is about 950 years')
    p.add_argument('--stride1D', type=int, default=10, help='steps between 1D outputs (strideSeries is 1)')
    p.add_argument('--lock-depth', type=float, default=None, help='lock both faults below this depth (km)')
    p.add_argument('--series-depth', type=float, default=None, help='depth of the root-strength probes in faultSeries.txt (km)')
    p.add_argument('--series-width', type=float, default=2.0, help='half-width of the probe average around each fault (km)')
    p.add_argument('--set', action='append', default=[], metavar='KEY=VALUE', help='append the input line "KEY = VALUE" (repeatable)')
    p.add_argument('--weak-band', type=float, default=None, help='power law: full width (km) of a weak creep band beneath each fault')
    p.add_argument('--weak-factor', type=float, default=1000.0, help='factor on the dislocation-creep prefactor in the band')
    p.add_argument('--weak-top', type=float, default=None, help='top of the band (km; default the lock depth, else 20)')
    p.add_argument('--weak-bottom', type=float, default=30.0, help='bottom of the band (km): deeper, the hot mantle is weak already and\n'
                   'a weaker band shortens the Maxwell time that caps the steps')
    p.add_argument('--weak-taper', type=float, default=0.5, help='width of the band edges (km)')
    p.add_argument('--mantle-zone', type=float, default=None, metavar='W',
                   help='with --weak-band: a mantle shear zone W km wide below the bands, centred between the faults\n'
                   '(--base-center), weakened like the bands from --weak-bottom to --mantle-bottom')
    p.add_argument('--mantle-bottom', type=float, default=40.0, help='bottom of the fully weak mantle zone (km)')
    p.add_argument('--mantle-taper', type=float, default=5.0, help='width of the mantle zone\'s sides and bottom edge (km)')
    p.add_argument('--moving-base', type=float, default=None, metavar='W',
                   help='move the bottom boundary with the plates (momBal_bcB_qd = movingBase), from -vL/2 to +vL/2 across\n'
                   'a transition W km wide (holding 90%% of the change; 0: a step); needs --lock-depth')
    p.add_argument('--base-center', type=float, default=None, help='centre of that transition (km; default midway between the faults)')
    p.add_argument('--name', default=None, help='run name (default from the options)')
    p.add_argument('--out', default='data/two_faults', help='directory for inputs, initial conditions and output')
    args = p.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    name = args.name or ('single' if args.single else '%s_d%g' % (args.f2, args.d)) + ('_pl' if args.rheology == 'powerlaw' else '')
    Ly = args.Ly
    y1, y2 = 0.5*(Ly - args.d), 0.5*(Ly + args.d)
    y, h0 = two_fault_grid(Ly, [y1, y2], args.h0, args.band, args.growth)
    if args.rheology == 'elastic':
        Lz, Nz = 30.0, args.Nz
        z = np.linspace(0, Lz, Nz)
    else:  # uniform to zfine, then geometric growth to Lz
        Lz = args.Lz
        nf = int(round(args.zfine/args.dz))
        z = np.concatenate([args.dz*np.arange(nf + 1), args.zfine + np.cumsum(geometric_cells(Lz - args.zfine, args.dz, args.growth))])
        z[-1] = Lz
        Nz = z.size

    # friction and steady state at vL (examples/ex2.in, Fault::computePsiSS)
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
    if args.rheology == 'powerlaw':
        write_vec(os.path.join(d, 'ic', 'z'), np.tile(z, y.size))
    write_vec(os.path.join(d, 'ic', 'psi'), psi1)
    write_vec(os.path.join(d, 'ic', 'prestress'), tau1)
    base = [l for l in open(os.path.join(root, 'examples/ex2.in')).read().splitlines() if not l.startswith((
        'outputDir', 'inputDir', 'guessSteadyStateICs', 'computeSSMomBal', 'maxTime', 'stride1D', 'stride2D',
        'restartFromChkpt', 'momBal_bc', 'Ny =', 'Nz =', 'Ly =', 'Lz =', 'bCoordTrans', 'maxStepCount', 'strideChkpt'))]
    lines = base + [
        'Ny = %d' % y.size, 'Nz = %d' % Nz, 'Ly = %g' % Ly, 'Lz = %g' % Lz,
        'bCoordTrans = -1 # the grid comes from the file y in inputDir',
        'outputDir = %s/' % d, 'inputDir = %s/ic/' % d,
        'guessSteadyStateICs = 0 # interior faults take their initial state from files or keys',
        'computeSSMomBal = 0', 'maxTime = %g' % args.maxTime, 'maxStepCount = 1e8',
        'stride1D = %d' % args.stride1D, 'stride2D = 0', 'strideSeries = 1', 'strideChkpt = 500',
        'restartFromChkpt = 1 # continue from data in outputDir if a run stopped',
        'momBal_bcL_qd = remoteLoading # full domain: the left side moves at -vL/2',
        'momBal_bcR_qd = remoteLoading', 'momBal_bcT_qd = freeSurface', 'momBal_bcB_qd = freeSurface']
    if args.rheology == 'powerlaw':
        pl = [l for l in open(os.path.join(root, 'examples/ex4s.in')).read().splitlines() if l.startswith((
            'disl_', 'maxEffVisc', 'kVals', 'kDepths', 'cVals', 'cDepths', 'withViscShearHeating', 'withFrictionalHeating'))]
        lines = [l for l in lines if not l.startswith(('bulkDeformationType', 'timeIntegrator'))]
        lines += pl + [
            'bulkDeformationType = powerLaw', 'momentumBalanceType = quasidynamic', 'systemEvolutionType = transient',
            'TVals = [283.15 1488.15 1623.15] # (K) the geotherm of ex4 (its LAB at 50 km)', 'TDepths = [0 50 500]',
            'thermalCoupling = coupled', 'evolveTemperature = 1', 'heatEquationType = transient',
            'withRadioHeatGeneration = no', 'timeIntegrator = RK43_WBE',
            'wVals = [%g %g] # (m) frictional heat: a Gaussian of this width on each fault' % (args.w, args.w), 'wDepths = [0 500]',
            'bcLType_trans = Dirichlet # heat: the far left side like the right']
        if args.grainsize:
            lines += ['evolveGrainSize = 1', 'grainSizeEvCoupling = uncoupled # no grain-size-sensitive creep law yet',
                      'grainSizeEv_grainSizeEvType = transient',
                      'grainSizeEv_grainSizeVals = [1e-5 1e-5] # initial grain size', 'grainSizeEv_grainSizeDepths = [0 500]',
                      'grainSizeEv_AVals = [1e-16 1e-16] # quartz, Tokle & Hirth (2021) (REVERSIBLE_STRENGTH_PLAN.md 4.3)',
                      'grainSizeEv_ADepths = [0 500]', 'grainSizeEv_QRVals = [16100 16100]', 'grainSizeEv_QRDepths = [0 500]',
                      'grainSizeEv_pVals = [3 3]', 'grainSizeEv_pDepths = [0 500]',
                      'grainSizeEv_gammaVals = [1 1]', 'grainSizeEv_gammaDepths = [0 500]', 'grainSizeEv_c = 3.14159',
                      'grainSizeEv_fVals = [0.015 0.015]', 'grainSizeEv_fDepths = [0 500]']
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
    if args.lock_depth is not None:
        lines += ['lockedVals = [0 0 1 1] # the faults end at %g km: the plate motion beneath flows viscously' % args.lock_depth,
                  'lockedDepths = [0 %g %g 500]' % (args.lock_depth, args.lock_depth)]
    if args.moving_base is not None:
        if args.lock_depth is None: raise SystemExit('--moving-base needs --lock-depth: the base pins the faults where they meet it')
        yc = args.base_center if args.base_center is not None else (y1 if args.single else 0.5*(y1 + y2))
        lines = [l for l in lines if not l.startswith('momBal_bcB_qd')] + [
            'momBal_bcB_qd = movingBase # the base moves with the plates, -vL/2 to +vL/2 across the transition',
            'momBal_bcB_center = %.12g # (km)' % yc,
            'momBal_bcB_width = %g # (km) holds 90%% of the velocity change' % args.moving_base]
    if args.series_depth is not None:
        lines += ['seriesDepth = %g # (km) root-strength probes in faultSeries.txt' % args.series_depth, 'seriesWidth = %g # (km)' % args.series_width]
    if args.weak_band is not None:
        if args.rheology != 'powerlaw': raise SystemExit('--weak-band needs --rheology powerlaw')
        A0 = float([l for l in lines if l.startswith('disl_AVals')][0].split('[')[1].split()[0])
        top = args.weak_top if args.weak_top is not None else (args.lock_depth if args.lock_depth is not None else 20.0)
        T = args.weak_taper
        def ramp(x):  # 0 for x <= 0, 1 for x >= 1, cosine between
            x = np.clip(x, 0.0, 1.0); return 0.5*(1.0 - np.cos(np.pi*x))
        faults_y = [y1] if args.single else [y1, y2]
        my = np.max([ramp((0.5*args.weak_band + T - np.abs(y - yk))/T) for yk in faults_y], axis=0)
        mz = ramp((z - (top - T))/T)*ramp((args.weak_bottom + T - z)/T)
        m = my[:, None]*mz[None, :]
        if args.mantle_zone is not None:  # one zone below both bands: weak from where they end, tapering into the asthenosphere
            yc = args.base_center if args.base_center is not None else (y1 if args.single else 0.5*(y1 + y2))
            Tm = args.mantle_taper
            myz = ramp((0.5*args.mantle_zone + Tm - np.abs(y - yc))/Tm)
            mzz = ramp((z - (args.weak_bottom - T))/T)*ramp((args.mantle_bottom + Tm - z)/Tm)
            m = np.maximum(m, myz[:, None]*mzz[None, :])
        write_vec(os.path.join(d, 'ic', 'disl_A'), (A0*args.weak_factor**m).ravel())  # index iy*Nz + iz
        lines.append('# weak band: ic/disl_A is %g x the prefactor within %g km of each fault from %g to %g km (tapers %g km)'
                     % (args.weak_factor, 0.5*args.weak_band, top, args.weak_bottom, T))
        if args.mantle_zone is not None:
            lines.append('# mantle zone: the same factor within %g km of y = %g km from %g to %g km (sides and bottom taper %g km)'
                         % (0.5*args.mantle_zone, yc, args.weak_bottom, args.mantle_bottom, Tm))
    for kv in args.set:
        k, v = kv.split('=', 1)
        lines.append('%s = %s' % (k.strip(), v.strip()))
    open(os.path.join(args.out, name + '.in'), 'w').write('\n'.join(lines) + '\n')
    dy = np.diff(y)
    print('wrote %s.in: Ny = %d (spacing %.3g km at the faults, %.3g km at most), Nz = %d (Lz = %g km), faults at y = %g and %g km'
          % (os.path.join(args.out, name), y.size, h0, dy.max(), Nz, Lz, y1, y2))


if __name__ == '__main__':
    main()
