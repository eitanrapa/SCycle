#!/usr/bin/env python3
"""Turn a run's checkpoint into initial-condition files that start a new run from its state.

A run stopped at a chosen moment (tools/batch.sh stop, or SIGTERM: it writes checkpoint.h5 at that
step) can seed a run with other settings, such as the transient-heating cycle of Allison & Dunham
(2021) started just after an earthquake at the end of a spin-up with the thermal anomaly held fixed.
The new run starts at time 0 from the checkpoint's fault slip and state, displacement, viscous
strains and stresses, effective viscosity and boundary displacements; momBal_bcR is the right
boundary's displacement at the checkpoint, which the new run continues with vL t/2.

Usage: ckpt_to_ics.py CHECKPOINT OUTPREFIX [--base BASEPREFIX]
  First copies the files BASEPREFIX* (the original run's inputDir files: grid, flow law,
  temperatures) to OUTPREFIX*, then writes OUTPREFIX<name> for each field of the checkpoint. Run the
  new input with inputDir = OUTPREFIX, guessSteadyStateICs = 0 and computeSSMomBal = 0.
Only the boundary fault is converted (interior faults have their own groups and prefixes).
"""
import argparse
import glob
import os
import shutil
import sys

import h5py
import numpy as np

YEAR = 3.15576e7
FIELDS = {                      # file read by the code: checkpoint dataset
    'psi': '/fault_qd/psi',
    'slip': '/fault_qd/slip',
    'momBal_bcL': '/momBal/bcL',
    'momBal_bcR': '/momBal/bcR',
    'u': '/momBal/u',
    'GVxy': '/momBal/gVxy',
    'GVxz': '/momBal/gVxz',
    'Sxy': '/momBal/sxy',
    'Sxz': '/momBal/sxz',
    'effVisc': '/momBal/effVisc',
}


def write_vec(path, v):
    """PETSc binary Vec: classid 1211214, length, big-endian doubles."""
    with open(path, 'wb') as fh:
        np.array([1211214, v.size], dtype='>i4').tofile(fh)
        np.asarray(v, dtype='>f8').tofile(fh)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('checkpoint')
    p.add_argument('outprefix')
    p.add_argument('--base', help='prefix of the files to copy first (the original inputDir)')
    a = p.parse_args()

    outdir = os.path.dirname(a.outprefix) or '.'
    os.makedirs(outdir, exist_ok=True)
    if a.base:
        files = [f for f in glob.glob(a.base + '*') if os.path.isfile(f)]
        if not files:
            sys.exit(f'no files match {a.base}*')
        for f in files:
            shutil.copyfile(f, a.outprefix + f[len(a.base):])
        print(f'copied {len(files)} files from {a.base}*')
    with h5py.File(a.checkpoint, 'r') as h:
        t = float(np.array(h['/time1D/time1D']).ravel()[0]) if '/time1D/time1D' in h else float('nan')
        V = np.abs(np.array(h['/fault_qd/slipVel']).ravel()) if '/fault_qd/slipVel' in h else None
        print(f'checkpoint at t = {t / YEAR:.3f} yr' + (f', max slip rate {V.max():.3g} m/s' if V is not None else ''))
        for name, key in FIELDS.items():
            if key not in h:
                print(f'  {name}: not in the checkpoint, skipped')
                continue
            v = np.array(h[key], dtype=float).ravel()
            write_vec(a.outprefix + name, v)
            print(f'  {name} <- {key} ({v.size} values)')


if __name__ == '__main__':
    main()
