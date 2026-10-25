#!/usr/bin/env python3
"""How two faults trade slip: earthquake sequences, slip shares and root persistence of two-fault runs.

Reads, for each run, the events.csv and partition.csv of tools/two_fault.py (run it with --zref 10
--window 1000 first) and the root probes of faultSeries.txt (strideSeries, --series-depth). Everything
after --spin years of spin-up:

  N eq, interval   earthquakes and their mean interval (yr)
  P(same)          fraction of consecutive earthquakes on the same fault: 0 for strict turns, 0.5 for
                   independent faults, 1 for one fault only
  run mean / max   consecutive earthquakes on one fault
  >75% 500 / 1000  fraction of the time one fault takes more than 75% of the slip at --zref over the
                   trailing 500 or 1000 yr (std: the spread of that share)
  phases, length   dominance phases with hysteresis on the 500-yr share (a fault dominant once above
                   0.7, the other once it falls below 0.3): completed phases and their mean length (yr)
  root kept        how long the root beneath a fault keeps most (> 50%) of the deep creep after the
                   fault's earthquake: the share of the two roots' strain at the probe depth over a
                   trailing 50 yr, sampled every 25 yr up to 3000 yr; median over earthquakes (yr)
  share            the first fault's share of all slip at --zref after the spin-up

Usage: trading_stats.py [LABEL=]RUN_DIR ... [--spin 2000] [--csv out.csv]
These are the definitions of the stage 5 second batch (docs/REVERSIBLE_STRENGTH_PLAN.md, item 10).
"""
import argparse
import csv
import os

import numpy as np

YEAR = 3.15576e7


def events(d, spin):
    ev = [(float(r['onset_yr']), r['fault']) for r in csv.DictReader(open(os.path.join(d, 'events.csv')))]
    return sorted(e for e in ev if e[0] > spin)


def partition(d):
    t, s1, s2 = [], [], []
    for r in csv.DictReader(open(os.path.join(d, 'partition.csv'))):
        t.append(float(r['time_yr'])); s1.append(float(r['fault_slip_m'])); s2.append(float(r['f2_slip_m']))
    return np.array(t), np.array(s1), np.array(s2)


def shares(t, s1, s2, w, spin):
    tc = np.arange(spin + w, t[-1], 5.0)
    d1 = np.interp(tc, t, s1) - np.interp(tc - w, t, s1)
    d2 = np.interp(tc, t, s2) - np.interp(tc - w, t, s2)
    return tc, d1 / (d1 + d2)


def root_kept(d, ev):
    """Median years the root beneath the ruptured fault keeps > 50% of the deep creep (nan without probes)."""
    p = os.path.join(d, 'faultSeries.txt')
    if not os.path.exists(p):
        return np.nan
    hdr = open(p).readline().lstrip('#').split()
    col = {h: i for i, h in enumerate(hdr)}
    keys = ['step', 'time(s)'] + [k for k in hdr if k.endswith('_dgVdev@25km(1e-3/s)')]
    if len(keys) != 4:
        return np.nan
    a = np.loadtxt(p, comments='#', usecols=[col[k] for k in keys])
    _, idx = np.unique(a[::-1, 0], return_index=True)      # restarts append: keep the last line per step
    a = a[::-1][idx]
    a = a[np.argsort(a[:, 0])]
    t, dt = a[:, 1] / YEAR, np.diff(a[:, 1])
    E = {name: np.concatenate([[0.0], np.cumsum(a[1:, i] * dt)]) for i, name in ((2, 'fault'), (3, 'f2'))}
    kept = []
    for ta, fa in ev:
        tt = np.arange(ta + 25.0, min(ta + 3000.0, t[-1]), 25.0)
        inc = {n: np.interp(tt, t, E[n]) - np.interp(tt - 50.0, t, E[n]) for n in E}
        own = inc[fa] / (inc['fault'] + inc['f2'])
        lost = np.flatnonzero(own < 0.5)
        if lost.size:
            kept.append(tt[lost[0]] - ta)
    return float(np.median(kept)) if kept else np.nan


def stats(d, spin):
    ev = events(d, spin)
    f = [e[1] for e in ev]
    same = float(np.mean([x == y for x, y in zip(f[:-1], f[1:])])) if len(f) > 1 else np.nan
    runs, n = [], 1
    for x, y in zip(f[:-1], f[1:]):
        if x == y:
            n += 1
        else:
            runs.append(n); n = 1
    runs.append(n)
    t, s1, s2 = partition(d)
    row = dict(t_end=t[-1], n_eq=len(ev), interval=float(np.mean(np.diff([e[0] for e in ev]))) if len(ev) > 1 else np.nan,
               p_same=same, run_mean=float(np.mean(runs)), run_max=max(runs))
    for w in (500.0, 1000.0):
        _, sh = shares(t, s1, s2, w, spin)
        row[f'std{w:.0f}'] = float(np.nanstd(sh))
        row[f'out{w:.0f}'] = float(np.nanmean((sh > 0.75) | (sh < 0.25)))
    tc, sh = shares(t, s1, s2, 500.0, spin)
    state, starts = None, []
    for ti, x in zip(tc, sh):
        new = 'fault' if x > 0.7 else ('f2' if x < 0.3 else state)
        if new != state and new is not None:
            starts.append(ti); state = new
    lens = np.diff(starts)
    row['phases'], row['phase_yr'] = len(lens), float(np.mean(lens)) if len(lens) else np.nan
    row['root_kept'] = root_kept(d, ev)
    i0 = np.searchsorted(t, spin)
    d1, d2 = s1[-1] - s1[i0], s2[-1] - s2[i0]
    row['share'] = d1 / (d1 + d2)
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('runs', nargs='+')
    p.add_argument('--spin', type=float, default=2000.0)
    p.add_argument('--csv')
    a = p.parse_args()
    rows = []
    print('%-16s %5s %5s %6s %8s %5s | %4s %5s %4s %5s | %6s %6s | %6s %5s' % (
        'run', 'N eq', 'dt', 'Psame', 'run mean', 'max', 'std5', '>75%', 'std1', '>75%', 'phases', 'length', 'kept', 'share'))
    for r in a.runs:
        label, d = r.split('=', 1) if '=' in r else (os.path.basename(r.rstrip('/')), r)
        if not all(os.path.exists(os.path.join(d, f)) for f in ('events.csv', 'partition.csv')):
            print('%-16s no events.csv or partition.csv (run tools/two_fault.py)' % label)
            continue
        s = stats(d, a.spin)
        s['run'] = label
        rows.append(s)
        print('%-16s %5d %5.0f %6.2f %8.1f %5d | %4.2f %5.2f %4.2f %5.2f | %6d %6.0f | %6.0f %5.2f%s' % (
            label, s['n_eq'], s['interval'], s['p_same'], s['run_mean'], s['run_max'], s['std500'], s['out500'],
            s['std1000'], s['out1000'], s['phases'], s['phase_yr'], s['root_kept'], s['share'],
            '' if s['t_end'] > 12000 else '   (to %.0f yr)' % s['t_end']))
    if a.csv and rows:
        with open(a.csv, 'w') as fh:
            w = csv.DictWriter(fh, fieldnames=['run'] + [k for k in rows[0] if k != 'run'])
            w.writeheader(); w.writerows(rows)
        print('wrote', a.csv)


if __name__ == '__main__':
    main()
