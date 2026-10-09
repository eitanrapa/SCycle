#!/usr/bin/env python3
"""Runs of Allison & Dunham (2021) that start from the paper's own limit-cycle states.

The paper's runs reached their limit cycles after 100 or more earthquake cycles (the published time
series start with about 950 to 2250 m of accumulated slip at 5 km), which would take weeks per case
here. Its published results (osf.io/nzmw2, simulationData.zip) hold complete states to start from:

- heating cases (LAB 50 km, lambda 0.37): the transient run with w = 0.3 m stores 2D snapshots; its
  first, 158 s after its start, gives the viscous strains, stresses and effective viscosity, and its
  fault series the slip, state and boundary displacements (the w = 0.3 and w = 1 m runs start from the
  same state, at the end of the paper's spin-up with the anomaly held at its steady state);
- elastic case: the first output of eqs_le_l0.37 gives slip, state and the boundary displacements.

The new runs start at time 0 from those states (guessSteadyStateICs = 0) and are compared with the
paper's next cycle: its ssT cycle (recurrence 354 yr), its w = 1 m transient cycle (earthquake at
352.84 yr) and its elastic cycle (earthquake at 266.81 yr). The temperature starts as the paper's
steady-state anomaly (the T file of its initial conditions), as in its transient runs.

Usage: seed_ad2021.py --osf <osf dir with inputFiles/ and outputFiles/data_eqs/> --out <batch dir>
Writes <out>/seed_heat_ic/, <out>/seed_le_ic/ and the inputs q1 (ssT cycle), q2 (elastic cycle),
q3 (transient heating, w = 1 m).
"""
import argparse
import glob
import os
import shutil

import h5py
import numpy as np

from make_ad2021 import SN, DC, YEAR, settings

CASE = 'zlab50_l0.37'


def write_vec(path, v):
    """PETSc binary Vec: classid 1211214, length, big-endian doubles."""
    with open(path, 'wb') as fh:
        np.array([1211214, v.size], dtype='>i4').tofile(fh)
        np.asarray(v, dtype='>f8').tofile(fh)


def seed(base, outprefix, fields):
    """Copy the files base*, except those a field replaces (compared without case: macOS), then write fields."""
    os.makedirs(os.path.dirname(outprefix), exist_ok=True)
    replaced = {k.lower() for k in fields}
    for f in glob.glob(base + '*'):
        name = f[len(base):]
        if os.path.isfile(f) and name.lower() not in replaced:
            shutil.copyfile(f, outprefix + name)
    for name, v in fields.items():
        write_vec(outprefix + name, np.ascontiguousarray(v, dtype=float).ravel())


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--osf', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    rel = lambda path: os.path.relpath(os.path.abspath(path), root)
    ics = os.path.join(a.osf, 'inputFiles', 'IC_files_eqs')
    src = os.path.join(a.osf, 'inputFiles', 'input_files_scycle')
    eqs = os.path.join(a.osf, 'outputFiles', 'data_eqs')

    # heating cases: the first 2D snapshot of the w = 0.3 m run and its first fault output
    pre = os.path.join(eqs, f'eqs_wfvsh_w0.3_{CASE}_')
    two_d = {}
    for name, (fname, key) in {'GVxy': ('GVxy', 'GVxy'), 'GVxz': ('GVxz', 'GVxz'), 'Sxy': ('Sxy', 'Sxy'),
                               'Sxz': ('Sxz', 'Sxz'), 'effVisc': ('EffVisc', 'EffVisc')}.items():
        with h5py.File(pre + fname + '.mat', 'r') as f:
            two_d[name] = np.array(f[key][0])          # (Ny, Nz): flattens z fastest, as the code orders
    with h5py.File(pre + 'onFault.mat', 'r') as f:
        fault = {'psi': np.array(f['psi'])[0], 'slip': np.array(f['slip'])[0],
                 'momBal_bcL': np.array(f['momBal/bcL'])[0], 'momBal_bcR': np.array(f['momBal/bcR'])[0]}
    seed(os.path.join(ics, f'eqs_wfvsh_{CASE}_'), os.path.join(a.out, 'seed_heat_ic', 'w0.3_t0_'), {**two_d, **fault})

    # elastic case: the first output of the paper's elastic run
    with h5py.File(os.path.join(eqs, 'eqs_le_l0.37.mat'), 'r') as f:
        le = {'psi': np.array(f['psi'])[0], 'slip': np.array(f['slip'])[0],
              'momBal_bcL': np.array(f['slip'])[0] / 2.0, 'momBal_bcR': np.array(f['momBal/bcR'])[0]}
    seed(os.path.join(ics, 'eqs_elastic_l0.37_'), os.path.join(a.out, 'seed_le_ic', 'le_t0_'), le)

    common = {'sNVals': SN['0.37'], 'DcVals': DC['0.37'], 'normType': 'L2_absolute', 'timeStepTol': '1e-7',
              'guessSteadyStateICs': '0', 'computeSSMomBal': '0', 'stride1D': '20', 'stride2D': '20000',
              'strideChkpt': '2000', 'restartFromChkpt': '1', 'initTime': '0'}
    runs = [
        ('q1_ad21_seed_ssT', 'viscoelastic_eq_sims.in', 'seed_heat_ic/w0.3_t0_',
         {'thermalCoupling': 'no', 'heatEquationType': 'steadyState', 'withViscShearHeating': 'no',
          'withFrictionalHeating': 'no', 'timeIntegrator': 'RK43', 'minDeltaT': '1e-3', 'maxDeltaT': '1e11',
          'maxTime': f'{800 * YEAR:.6g}'}),
        ('q2_ad21_seed_le', 'elastic_eq_sims.in', 'seed_le_ic/le_t0_',
         {'timeIntegrator': 'RK43', 'minDeltaT': '1.78e-5', 'maxDeltaT': '1e11', 'maxTime': f'{600 * YEAR:.6g}'}),
        ('q3_ad21_seed_w1', 'viscoelastic_eq_sims.in', 'seed_heat_ic/w0.3_t0_',
         {'thermalCoupling': 'coupled', 'heatEquationType': 'transient', 'withViscShearHeating': 'yes',
          'withFrictionalHeating': 'yes', 'wVals': '[1 1]', 'evolveTemperature': '1', 'timeIntegrator': 'RK43_WBE',
          'minDeltaT': '1e-3', 'maxDeltaT': '2.5e5', 'stride2D': '200', 'maxTime': f'{370 * YEAR:.6g}'}),
    ]
    for tag, base, icprefix, extra in runs:
        changes = dict(common)
        changes.update(extra)
        changes['inputDir'] = rel(os.path.join(a.out, icprefix))
        changes['outputDir'] = rel(os.path.join(a.out, tag)) + '/'
        os.makedirs(os.path.join(a.out, tag), exist_ok=True)
        with open(os.path.join(a.out, tag + '.in'), 'w') as f:
            f.write(settings(open(os.path.join(src, base)).read(), changes))
        print(os.path.join(a.out, tag + '.in'))


if __name__ == '__main__':
    main()
