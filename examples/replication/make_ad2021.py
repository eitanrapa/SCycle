#!/usr/bin/env python3
"""Inputs that replicate Allison & Dunham (2021, JGR Solid Earth 126, e2020JB021394) with this code.

The paper's input files, grids and steady-state initial conditions are public at osf.io/nzmw2
(inputFiles.zip). This script takes its three input files and writes inputs the current code runs,
with the settings the paper's own output files record (their /med and /he groups), which differ in
places from the posted inputs:

- normType = L2_absolute and timeStepTol = 1e-7 (the posted files say max_relative and 1e-8; the
  implicit-explicit integrator accepts only L2 norms);
- minDeltaT = 1e-3 s for the power-law runs and 1.78e-5 s, one shear-wave transit of the smallest
  cell, for the elastic run; maxDeltaT = 1e11 s, as in the paper's no-heating and ssT runs. The
  paper's elastic and transient-heating runs capped the step at 2.5e5 s, which only adds output
  between earthquakes (error control at 1e-7 sets the accuracy) and here would make an elastic cycle
  take 18 hours; the elastic run below does not cap it.

Changes the current code needs:

- evolveTemperature = 1 with transient heating (the 2020 code inferred it; since 2022 a coupled run
  must say so);
- computeSSMomBal = 1, so that guessSteadyStateICs computes the steady-state momentum balance as the
  2020 code always did;
- the "ssT" runs, the spin-up with the temperature anomaly held at its steady state, use
  thermalCoupling = no: the power law then keeps the temperature of the case's T file, which holds the
  ambient geotherm plus the steady-state anomaly. The paper's ssT runs were coupled but used the
  explicit RK43, which in the 2020 code left the temperature at its initial value; the current code
  refuses that pairing.

Cases: elastic (le), viscoelastic without shear heating (nsh), the ssT spin-up (ssT). The run with
transient shear heating starts from the end of the ssT spin-up and is written separately.

Usage: make_ad2021.py --osf <dir with input_files_scycle and IC_files_eqs> --out <batch dir>
       [--zlab 50] [--lam 0.37] [--prefix r] [--years-le 2000] [--years-ve 7500]
Paths are written relative to the repository root, where tools/batch.sh runs the inputs.
"""
import argparse
import os
import re

YEAR = 3.15576e7
SN = {'0.37': '[0.00 6128.6]', '0.6': '[0.00 3891.18]', '0.8': '[0.00 1945.59]'}
DC = {'0.37': '[0.05 0.05 0.01 0.01]', '0.6': '[0.0317 0.0317 0.0063 0.0063]',
      '0.8': '[0.0159 0.0159 0.0032 0.0032]'}


def settings(text, changes):
    """Set each key of changes (replacing an active line, or appending one); drop commented copies."""
    lines, done = [], set()
    for line in text.splitlines():
        m = re.match(r'^\s*(#?)\s*([A-Za-z_0-9]+) = ', line)
        if m and m.group(2) in changes:
            if m.group(1):          # a commented alternative
                continue
            if m.group(2) not in done:
                lines.append(f'{m.group(2)} = {changes[m.group(2)]}')
                done.add(m.group(2))
            continue
        lines.append(line)
    lines += ['', '# set by examples/replication/make_ad2021.py']
    lines += [f'{k} = {v}' for k, v in changes.items() if k not in done]
    return '\n'.join(lines) + '\n'


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--osf', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--zlab', default='50')
    p.add_argument('--lam', default='0.37', choices=sorted(SN))
    p.add_argument('--prefix', default='r')
    p.add_argument('--years-le', type=float, default=2000.0)
    p.add_argument('--years-ve', type=float, default=7500.0)
    a = p.parse_args()

    root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    rel = lambda path: os.path.relpath(os.path.abspath(path), root)
    src = os.path.join(a.osf, 'input_files_scycle')
    ics = os.path.join(a.osf, 'IC_files_eqs')
    os.makedirs(a.out, exist_ok=True)
    case = f'zlab{a.zlab}_l{a.lam}'
    common = {'sNVals': SN[a.lam], 'DcVals': DC[a.lam], 'normType': 'L2_absolute',
              'timeStepTol': '1e-7', 'computeSSMomBal': '1', 'stride1D': '20', 'stride2D': '20000',
              'strideChkpt': '2000', 'restartFromChkpt': '1'}
    noheat = {'thermalCoupling': 'no', 'heatEquationType': 'steadyState',
              'withViscShearHeating': 'no', 'withFrictionalHeating': 'no', 'timeIntegrator': 'RK43',
              'minDeltaT': '1e-3', 'maxDeltaT': '1e11', 'maxTime': f'{a.years_ve * YEAR:.6g}'}
    runs = [
        ('1', f'ad21_ssT_{case}', 'viscoelastic_eq_sims.in', f'eqs_wfvsh_{case}_', noheat),
        ('3', f'ad21_nsh_{case}', 'viscoelastic_eq_sims.in', f'eqs_nsh_{case}_', noheat),
        ('5', f'ad21_le_l{a.lam}', 'elastic_eq_sims.in', f'eqs_elastic_l{a.lam}_',
         {'timeIntegrator': 'RK43', 'minDeltaT': '1.78e-5', 'maxDeltaT': '1e11',
          'maxTime': f'{a.years_le * YEAR:.6g}'}),
    ]
    for order, name, base, icprefix, extra in runs:
        tag = f'{a.prefix}{order}_{name}'
        changes = dict(common)
        changes.update(extra)
        changes['inputDir'] = rel(os.path.join(ics, icprefix))
        changes['outputDir'] = rel(os.path.join(a.out, tag)) + '/'
        text = settings(open(os.path.join(src, base)).read(), changes)
        os.makedirs(os.path.join(a.out, tag), exist_ok=True)
        with open(os.path.join(a.out, tag + '.in'), 'w') as f:
            f.write(text)
        print(os.path.join(a.out, tag + '.in'))


if __name__ == '__main__':
    main()
