#!/usr/bin/env python3
"""Summarize a multi-fault SCycle run: earthquakes, slip partitioning, surface velocity.

Reads data_context.h5, data_1D.h5, faultSeries.txt and mediator.txt from a run directory (the
outputDir of StrikeSlip_LinearElastic_qd with interior faults; see examples/two_faults) and writes:

  events.csv     one row per earthquake: fault, onset and end (s and yr), duration (s), peak slip rate
                 (m/s), depth of the slip-rate maximum at onset (km), potency (m^2, slip
                 integrated over depth) and moment per unit fault length (N m/m)
  partition.csv  per 1D output: time (yr), each fault's slip at --zref (m), and each fault's
                 share of the slip at --zref over the trailing --window years
  surfvel.csv    surface position y (km), the mean interseismic surface velocity (weighted by time)
                 and the latest interseismic one (mm/yr), from outputs where every fault slips slower
                 than --vinter
  switches.csv   one row per change of the dominant fault: time (yr), the faults it passes from and to,
                 the length of the phase that ended (yr) and each fault's slip at --zref during it (m).
                 The dominant fault has the largest share; a switch needs the new one's share above
                 0.5 + --margin, so fluctuations about equal sharing are not switches.

and prints a summary with the long-term slip rates and their sum against the plate rate vL.

  python3 tools/two_fault.py data/two_faults/rs_d20/ [--zref 5] [--window 300]

Needs numpy and h5py.
"""
import argparse, os, re
import numpy as np
import h5py

YEAR = 365.25*24*3600.0


def series(run):
    """faultSeries.txt as {column: array}, keeping the last line of each step (restarts append)."""
    path = os.path.join(run, 'faultSeries.txt')
    with open(path) as f:
        header = f.readline().lstrip('#').split()
        rows = {}
        for line in f:
            if line.strip() and not line.startswith('#'):
                v = line.split()
                rows[int(v[0])] = [float(x) for x in v]
    data = np.array([rows[k] for k in sorted(rows)])
    names = [re.sub(r'\(.*\)$', '', h) for h in header]
    return {n: data[:, i] for i, n in enumerate(names)}


def events(t, maxV, zMaxV, potency, vseis, merge):
    """Intervals with maxV > vseis, merging those less than merge seconds apart."""
    above = np.flatnonzero(maxV > vseis)
    if above.size == 0:
        return []
    runs, start, prev = [], above[0], above[0]
    for i in above[1:]:
        if i != prev + 1 and t[i] - t[prev] > merge:
            runs.append((start, prev)); start = i
        prev = i
    runs.append((start, prev))
    out = []
    for i0, i1 in runs:
        before = max(i0 - 1, 0)
        out.append(dict(onset=t[i0], end=t[i1], peakV=maxV[i0:i1 + 1].max(), z=zMaxV[i0],
                        potency=potency[min(i1 + 1, t.size - 1)] - potency[before]))
    return out


def switches(t, shares, slip, names, margin):
    """Changes of the dominant fault: rows (time, from, to, phase length, slip of each fault in the phase).
    shares[n][j] is fault n's share at output j (nan before the window fills); the dominant fault is
    the first whose share exceeds 0.5 + margin, and keeps dominance until another one does."""
    rows, dom, start = [], None, None
    for j in range(t.size):
        new = None
        for n in names:
            if np.isfinite(shares[n][j]) and shares[n][j] > 0.5 + margin and n != dom:
                new = n
                break
        if new is None: continue
        if dom is not None:
            rows.append((t[j], dom, new, t[j] - t[start], {n: slip[n][j] - slip[n][start] for n in names}))
        dom, start = new, j
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('run', help='run directory (the outputDir prefix, ending in /)')
    p.add_argument('--zref', type=float, default=5.0, help='depth for slip histories and partitioning (km)')
    p.add_argument('--window', type=float, default=300.0, help='trailing window for slip shares (yr)')
    p.add_argument('--vseis', type=float, default=1e-3, help='slip rate above which a fault is in an earthquake (m/s)')
    p.add_argument('--merge', type=float, default=3600.0, help='merge seismic intervals closer than this (s)')
    p.add_argument('--vinter', type=float, default=1e-8, help='interseismic: every fault slower than this (m/s)')
    p.add_argument('--margin', type=float, default=0.1, help='share above 0.5 that a new dominant fault must reach')
    args = p.parse_args()
    run = args.run if args.run.endswith('/') else args.run + '/'

    ctx = h5py.File(run + 'data_context.h5', 'r')
    d1 = h5py.File(run + 'data_1D.h5', 'r')
    s = series(run)
    names = [k[:-len('_maxV')] for k in s if k.endswith('_maxV')]
    vL = 1e-9
    for line in open(run + 'mediator.txt'):
        if line.startswith('vL = '): vL = float(line.split()[2])

    # faults: position, depth grid, shear modulus
    info = {}
    for n in names:
        g = ctx[n]
        info[n] = dict(y=float(g.attrs['y']) if 'y' in g.attrs else 0.0, z=np.array(g['z']).ravel(),
                       mu=np.array(ctx[n + '_qd/mu']).ravel())
    t1 = np.array(d1['time/time1D']).ravel()
    print('%s: %d steps to t = %.4g yr, %d 1D outputs; vL = %g m/s' % (run, s['step'].size, s['time'][-1]/YEAR, t1.size, vL))
    for n in names:
        print('  fault %-8s y = %8.3f km' % (n, info[n]['y']))

    # earthquakes, from every step
    cat = []
    for n in names:
        mu = info[n]['mu'].mean()*1e9  # Pa (uniform in the examples)
        for e in events(s['time'], s[n + '_maxV'], s[n + '_zMaxV'], s[n + '_potency'], args.vseis, args.merge):
            e.update(fault=n, moment=mu*e['potency']); cat.append(e)
    cat.sort(key=lambda e: e['onset'])
    with open(run + 'events.csv', 'w') as f:
        f.write('fault,onset_s,end_s,onset_yr,end_yr,duration_s,peakV_m_per_s,depth_km,potency_m2,moment_Nm_per_m\n')
        for e in cat:
            f.write('%s,%.15e,%.15e,%.9g,%.9g,%.6g,%.6g,%.4g,%.6g,%.6g\n' % (e['fault'], e['onset'], e['end'], e['onset']/YEAR, e['end']/YEAR,
                    e['end'] - e['onset'], e['peakV'], e['z'], e['potency'], e['moment']))
    print('\n%d earthquakes (slip rate above %g m/s):' % (len(cat), args.vseis))
    print('  %-8s %12s %10s %9s %8s %12s' % ('fault', 'onset (yr)', 'dur (s)', 'peakV', 'z (km)', 'potency (m2)'))
    for e in cat[:40]:
        print('  %-8s %12.4f %10.1f %9.3g %8.2f %12.4g' % (e['fault'], e['onset']/YEAR, e['end'] - e['onset'], e['peakV'], e['z'], e['potency']))
    if len(cat) > 40: print('  ... (%d more in events.csv)' % (len(cat) - 40))

    # slip at zref and the shares of each fault over a trailing window, from the 1D outputs
    slip = {}
    for n in names:
        z = info[n]['z']
        S = np.array(d1[n + '/slip'])[:, :, 0]
        slip[n] = np.array([np.interp(args.zref, z, row) for row in S])
    W = args.window*YEAR
    j0 = np.searchsorted(t1, t1 - W)  # output at or after t - W
    total = sum(slip[n] - slip[n][j0] for n in names)
    with open(run + 'partition.csv', 'w') as f:
        f.write('time_yr,' + ','.join('%s_slip_m' % n for n in names) + ',' + ','.join('%s_share' % n for n in names) + '\n')
        share = {n: np.full(t1.size, np.nan) for n in names}
        for j in range(t1.size):
            shares = [(slip[n][j] - slip[n][j0[j]])/total[j] if total[j] > 0 and t1[j] - t1[j0[j]] > 0.5*W else np.nan for n in names]
            for n, x in zip(names, shares): share[n][j] = x
            f.write('%.9g,' % (t1[j]/YEAR) + ','.join('%.9g' % slip[n][j] for n in names) + ',' + ','.join('%.6g' % x for x in shares) + '\n')
    # switches of the dominant fault
    sw = switches(t1, share, slip, names, args.margin)
    with open(run + 'switches.csv', 'w') as f:
        f.write('time_yr,from,to,phase_yr,' + ','.join('%s_slip_m' % n for n in names) + '\n')
        for (tj, a, b, dur, ds) in sw:
            f.write('%.9g,%s,%s,%.9g,' % (tj/YEAR, a, b, dur/YEAR) + ','.join('%.6g' % ds[n] for n in names) + '\n')
    if sw:
        print('\n%d switches of the dominant fault (share above %.2f, window %g yr):' % (len(sw), 0.5 + args.margin, args.window))
        for (tj, a, b, dur, ds) in sw[:20]:
            print('  %10.2f yr  %s -> %s after %8.2f yr; slip in the phase: %s' % (tj/YEAR, a, b, dur/YEAR, ', '.join('%s %.2f m' % (n, ds[n]) for n in names)))
        if len(sw) >= 3:
            per = [sw[k + 2][0] - sw[k][0] for k in range(len(sw) - 2)]  # back to the same fault
            print('  alternation period %.2f yr (mean of %d); slip per phase on the dominant fault %.2f m (mean)'
                  % (np.mean(per)/YEAR, len(per), np.mean([ds[a] for (tj, a, b, dur, ds) in sw])))
    else:
        print('\nno switch of the dominant fault (share above %.2f, window %g yr)' % (0.5 + args.margin, args.window))
    # long-term rates: between the first and last onsets of the most active fault, if it has 2 or more
    # events spanning at least a quarter of the run (whole cycles), else over the last window
    ev = {n: [e for e in cat if e['fault'] == n] for n in names}
    main_fault = max(names, key=lambda n: len(ev[n]))
    span = ev[main_fault][-1]['onset'] - ev[main_fault][0]['onset'] if len(ev[main_fault]) >= 2 else 0.0
    if len(ev[main_fault]) >= 2 and span >= 0.25*(t1[-1] - t1[0]):
        ta, tb = ev[main_fault][0]['onset'], ev[main_fault][-1]['onset']
        how = 'between the first and last onsets on %s (%d cycles)' % (main_fault, len(ev[main_fault]) - 1)
    else:
        ta, tb = max(t1[-1] - W, t1[0]), t1[-1]
        how = 'over the last %.0f yr' % ((tb - ta)/YEAR)
    ia, ib = np.searchsorted(t1, ta), np.searchsorted(t1, tb)
    ib = min(ib, t1.size - 1)
    print('\nlong-term slip rates at z = %g km, %s (%.4g to %.4g yr):' % (args.zref, how, t1[ia]/YEAR, t1[ib]/YEAR))
    rates = {n: (slip[n][ib] - slip[n][ia])/(t1[ib] - t1[ia]) for n in names}
    for n in names:
        print('  %-8s %.4g m/s = %.4g vL, share %.4f' % (n, rates[n], rates[n]/vL, rates[n]/sum(rates.values())))
    print('  sum      %.4g vL' % (sum(rates.values())/vL))

    # interseismic surface velocity
    wrote = ['events.csv', 'partition.csv', 'switches.csv']
    if 'momBal/surfVel' in d1:
        V = np.array(d1['momBal/surfVel'])[:, :, 0]*1e3*YEAR  # mm/yr
        yb = np.array(ctx['domain/y']).ravel()
        Ny = V.shape[1]
        ys = yb.reshape(Ny, -1)[:, 0]
        vmax = np.max([np.array(d1[n + '/slipVel'])[:, :, 0].max(axis=1) for n in names], axis=0)
        first = cat[0]['end'] if cat else t1[0]
        sel = np.flatnonzero((vmax < args.vinter) & (t1 > first))
        if sel.size:
            # each output stands for the time halfway to its neighbours: outputs come every stride1D
            # steps, and a plain mean would weight the short steps after earthquakes far beyond their time
            edges = np.concatenate([[t1[0]], 0.5*(t1[1:] + t1[:-1]), [t1[-1]]])
            w = np.clip(np.diff(edges), 0.0, None)[sel]
            if w.sum() <= 0.0: w = np.ones(sel.size)
            with open(run + 'surfvel.csv', 'w') as f:
                f.write('y_km,mean_interseismic_mm_per_yr,latest_interseismic_mm_per_yr\n')
                mean, last = (w[:, None]*V[sel]).sum(axis=0)/w.sum(), V[sel[-1]]
                for k in range(Ny):
                    f.write('%.9g,%.6g,%.6g\n' % (ys[k], mean[k], last[k]))
            wrote.append('surfvel.csv')
            print('\ninterseismic surface velocity: %d outputs (%.4g yr) with every fault below %g m/s; far-field %.3f and %.3f mm/yr (vL/2 = %.3f)'
                  % (sel.size, w.sum()/YEAR, args.vinter, mean[0], mean[-1], 0.5*vL*1e3*YEAR))
        else:
            print('\nno interseismic outputs (every fault below %g m/s after the first earthquake)' % args.vinter)
    print('\nwrote %s in %s' % (', '.join(wrote), run))


if __name__ == '__main__':
    main()
