#!/usr/bin/env python3
"""Write an interior-fault example and its half-space twin (numpy only).

Both models use the friction, material and depth grid of examples/ex2.in and start at steady
sliding at vL, except for a smooth patch near 8 km depth that starts at 1e-6 m/s, so the first
event nucleates within about ten hours, the same way in every run.

- half: the usual model, the fault on the left boundary (symmFault), y in [0, 100] km with
  ex2's sinh stretching, remote loading on the right.
- full: the full domain y in [0, 200] km, the fault inside it at y = 100 km (interiorFaults),
  remote loading on both sides, and the grid mirrored about the fault so the spacing there
  matches the half-space. The fault lies midway between grid rows Ny/2 - 1 and Ny/2.

By antisymmetry the two describe the same physics: the full domain's right half is the half
space. With --ny 301 (and 602 rows in the full domain) the first event's onset, peak slip rate
and final slip agree within 0.14%, 0.3% and 0.03% (docs/TWO_FAULT_DESIGN.md, stage 2).

Run from the repository root:
  python3 examples/interior_fault/make_inputs.py [--ny 151]
  ./source/main data/interior_fault/full.in
  ./source/main data/interior_fault/half.in
"""
import argparse, os
import numpy as np

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--ny', type=int, default=151, help='rows of the half-space grid (the full domain gets 2*(ny-1)+2)')
p.add_argument('--out', default='data/interior_fault', help='directory for inputs, initial conditions and output')
p.add_argument('--maxTime', type=float, default=3.0e8, help='final time (s)')
args = p.parse_args()

root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
Nz, Lz, Ly, bct = 153, 30.0, 100.0, 5.0
z = np.linspace(0, Lz, Nz)
a = np.interp(z, [0, 11.0, 19.3333, 60], [0.0135, 0.0300, 0.0700, 0.2652])
b = np.interp(z, [0, 10, 11, 19.3333, 60], [0.0230, 0.0379, 0.0350, 0.0375, 0.0497])
f0, v0, vL, sN, V1 = 0.6, 1e-6, 1e-9, 50.0, 1e-6
eta = 0.5*np.sqrt(30.0*3.0)                      # radiation damping sqrt(mu*rho)/2 (mu 30 GPa, rho 3 g/cm^3)
psi_ss = f0 - b*np.log(vL/v0)                     # steady state at vL (Fault::computePsiSS)
prestress = a*sN*np.arcsinh(vL/(2*v0)*np.exp(psi_ss/a)) + eta*vL
w = np.exp(-((z - 8.0)/1.5)**2)                   # trigger patch
psi = psi_ss - a*w*np.log(V1/vL)                  # initial V = vL*(V1/vL)^w at the same stress

def write_vec(path, v):
    """PETSc binary Vec: classid 1211214, length, big-endian doubles."""
    with open(path, 'wb') as fh:
        np.array([1211214, v.size], dtype='>i4').tofile(fh)
        np.asarray(v, dtype='>f8').tofile(fh)

NyH = args.ny
NyF = 2*(NyH - 1) + 2
q = np.arange(NyF)/(NyF - 1.0)
yF = Ly + Ly*np.sinh(bct*(2*q - 1))/np.sinh(bct)  # ex2's stretching mirrored about the fault

base = [l for l in open(os.path.join(root, 'examples/ex2.in')).read().splitlines() if not l.startswith((
    'outputDir', 'inputDir', 'guessSteadyStateICs', 'computeSSMomBal', 'maxTime', 'stride1D', 'stride2D',
    'restartFromChkpt', 'momBal_bc', 'Ny =', 'Ly =', 'bCoordTrans', 'maxStepCount'))]
for kind in ['half', 'full']:
    d = os.path.join(args.out, kind)
    os.makedirs(os.path.join(d, 'ic'), exist_ok=True)
    write_vec(os.path.join(d, 'ic', 'prestress'), prestress)
    write_vec(os.path.join(d, 'ic', 'psi'), psi)
    lines = base + ['outputDir = %s/' % d, 'inputDir = %s/ic/' % d,
                    'guessSteadyStateICs = 0 # interior faults take their initial state from files or keys',
                    'computeSSMomBal = 0', 'maxTime = %g' % args.maxTime, 'stride1D = 10', 'stride2D = 0',
                    'restartFromChkpt = 0', 'maxStepCount = 1e8',
                    'momBal_bcR_qd = remoteLoading', 'momBal_bcT_qd = freeSurface', 'momBal_bcB_qd = freeSurface']
    if kind == 'half':
        lines += ['Ny = %d' % NyH, 'Ly = %g' % Ly, 'bCoordTrans = %g' % bct, 'momBal_bcL_qd = symmFault']
    else:
        write_vec(os.path.join(d, 'ic', 'y'), np.repeat(yF, Nz))  # body field, index iy*Nz + iz
        lines += ['Ny = %d' % NyF, 'Ly = %g' % (2*Ly),
                  'momBal_bcL_qd = remoteLoading # full domain: the left side moves at -vL/2',
                  'interiorFaults = [fault] # named "fault", so it reads the plain keys and files like the boundary fault',
                  'fault_y = %g # (km) placed midway between the nearest grid rows' % Ly]
    open(os.path.join(args.out, kind + '.in'), 'w').write('\n'.join(lines) + '\n')
print('wrote %s/half.in (Ny = %d) and %s/full.in (Ny = %d), fault spacing %.4f km'
      % (args.out, NyH, args.out, NyF, yF[NyF//2] - yF[NyF//2 - 1]))
