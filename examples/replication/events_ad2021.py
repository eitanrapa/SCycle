#!/usr/bin/env python3
"""Earthquakes of a run and of the published runs of Allison & Dunham (2021), measured the same way.

For each earthquake (the maximum slip rate rising through 1 mm/s), as defined in the paper and its
summaryData.xlsx: onset time, nucleation depth (where the slip rate first reaches 1 mm/s), coseismic
slip (from the last output before the onset to the first after the slip rate falls back below 1 mm/s
everywhere), its maximum, the down-dip limit (the depth above which 98% of the coseismic potency
lies), the moment per unit length log10(mu * integral of slip dz) and the recurrence interval.

Usage: events_ad2021.py RUN_OR_MAT [RUN_OR_MAT ...] [--mu 32.4]
  RUN_OR_MAT is a run's outputDir (with data_1D.h5) or one of the paper's .mat files (data_eqs/*.mat,
  or a *_onFault.mat).
"""
import argparse
import os

import h5py
import numpy as np

YEAR = 3.15576e7
V_EQ = 1e-3   # m/s


def load(path):
    """Time (s), depth (km), slip (m) and slip rate (m/s) as (time, depth) arrays."""
    if os.path.isdir(path):
        with h5py.File(os.path.join(path, 'data_1D.h5'), 'r') as f:
            t = f['/time/time1D'][:, 0, 0]
            slip = f['/fault/slip'][:, :, 0]
            V = f['/fault/slipVel'][:, :, 0]
        with h5py.File(os.path.join(path, 'data_context.h5'), 'r') as f:
            z = np.array(f['/fault/z']).ravel() if '/fault/z' in f else None
        if z is None:
            with h5py.File(os.path.join(path, 'checkpoint.h5'), 'r') as f:
                z = np.array(f['/fault_qd/z']).ravel()
    else:
        with h5py.File(path, 'r') as f:
            t = np.array(f['time']).ravel()
            slip, V = np.array(f['slip']), np.array(f['slipVel'])
            z = np.array(f['z'])[0] if np.array(f['z']).ndim == 2 else np.array(f['z']).ravel()
    return t, z, slip, V


def events(t, z, slip, V, mu, merge=YEAR):
    """Earthquakes; rises through 1 mm/s less than `merge` seconds after the previous one join it."""
    vmax = np.abs(V).max(axis=1)
    on = np.where((vmax[1:] >= V_EQ) & (vmax[:-1] < V_EQ))[0] + 1
    groups = []
    for i in on:
        if groups and t[i] - t[groups[-1][-1]] < merge:
            groups[-1].append(i)
        else:
            groups.append([i])
    out = []
    for g in groups:
        i, j = g[0], g[-1]
        while j < t.size and vmax[j] >= V_EQ:
            j += 1
        if j >= t.size:
            break                                   # still slipping at the end of the record
        d = slip[j] - slip[i - 1]                   # coseismic slip
        zz = z * 1e3
        pot = np.concatenate([[0.0], np.cumsum(0.5 * (d[1:] + d[:-1]) * np.diff(zz))])
        downdip = np.interp(0.98 * pot[-1], pot, z)
        nuc = z[np.argmax(np.abs(V[i]))]
        out.append(dict(t=t[i], nuc=nuc, dmax=d.max(), downdip=downdip,
                        logM=np.log10(mu * 1e9 * pot[-1]), dsurf=d[0]))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('runs', nargs='+')
    p.add_argument('--mu', type=float, default=32.4, help='shear modulus (GPa)')
    a = p.parse_args()
    for path in a.runs:
        t, z, slip, V = load(path)
        ev = events(t, z, slip, V, a.mu)
        print(f'== {path}: {t.size} outputs, t = {t[0] / YEAR:.2f} .. {t[-1] / YEAR:.2f} yr, {len(ev)} earthquakes')
        prev = None
        for e in ev:
            rec = f'{(e["t"] - prev) / YEAR:8.2f}' if prev is not None else '       -'
            print(f'   onset {e["t"] / YEAR:9.2f} yr  recurrence {rec}  nucleation {e["nuc"]:6.2f} km  '
                  f'max slip {e["dmax"]:6.3f} m (surface {e["dsurf"]:6.3f})  down-dip {e["downdip"]:6.2f} km  log10 M {e["logM"]:.4f}')
            prev = e['t']


if __name__ == '__main__':
    main()
