#!/usr/bin/env python3
"""Two shear-zone roots in series: do they share, lock in, or alternate?

Both roots carry the same far-field stress and their slip rates add to the plate rate, so the
partition is f_1 = S_1^-n / (S_1^-n + S_2^-n). Each root has a slow strengthening memory h (grows
with strain, relaxes in time) and a fast weakening feedback D (grows with strain, heals in time):

    S = (1 + h) exp(-beta D),   dh/dt = f - h/Th,   dD/dt = R f - D/TD

in time units of w*gamma_h/vL (w shear-zone width, gamma_h hardening strain). Th and TD are the
recovery and healing times in those units, R = gamma_h/gamma_D. Equal sharing is unstable when
N_loc = n beta R TD/2 > 1; lock-in is permanent when N_loc > (n/2) ln(1 + Th), roughly.

  python3 tools/zero_d_trading.py --n 3 --beta 0.4 --Th 10 --TD 0.3 --R 10
  python3 tools/zero_d_trading.py --scan      # the regime table of docs/REVERSIBLE_STRENGTH_PLAN.md
"""
import argparse, itertools
import numpy as np


def integrate(n, beta, Th, TD, R, T=None, dt=None, eps=1e-3):
    """RK4 from the symmetric steady state with a relative perturbation eps of D. Returns t, f_1."""
    n, beta, Th, TD, R = (np.atleast_1d(np.asarray(v, float)) for v in (n, beta, Th, TD, R))
    M = max(v.size for v in (n, beta, Th, TD, R))
    n, beta, Th, TD, R = (np.broadcast_to(v, M).copy() for v in (n, beta, Th, TD, R))
    h = np.stack([0.5*Th, 0.5*Th]); D = np.stack([0.5*R*TD*(1 + eps), 0.5*R*TD*(1 - eps)])

    def rhs(h, D):
        S = (1 + h)*np.exp(-beta*D); w = S**-n; f = w/w.sum(axis=0)
        return f - h/Th, R*f - D/TD, f[0]

    T = T or 30*Th.max(); dt = dt or 0.01*min(1.0, TD.min(), 1.0/R.max())
    N = int(T/dt); every = max(1, N//20000); ts, fs = [], []
    for k in range(N):
        k1h, k1D, f1 = rhs(h, D); k2h, k2D, _ = rhs(h + 0.5*dt*k1h, D + 0.5*dt*k1D)
        k3h, k3D, _ = rhs(h + 0.5*dt*k2h, D + 0.5*dt*k2D); k4h, k4D, _ = rhs(h + dt*k3h, D + dt*k3D)
        h = h + dt*(k1h + 2*k2h + 2*k3h + k4h)/6; D = D + dt*(k1D + 2*k2D + 2*k3D + k4D)/6
        if k % every == 0:
            ts.append(k*dt); fs.append(f1.copy())
    return np.array(ts), np.array(fs)


def classify(t, f, Th):
    """Regime from the second half of the series: sharing, lock-in or alternation (with period)."""
    half = t > 0.5*t[-1]; f = f[half]; t = t[half]
    amp = f.max() - f.min(); cross = np.flatnonzero(np.diff(np.sign(f - 0.5)) != 0)
    if amp < 0.05 and abs(f[-1] - 0.5) < 0.05: return 'sharing', amp, float('nan'), f.max(), f.min()
    if amp < 0.05: return 'lock-in', amp, float('nan'), f.max(), f.min()
    per = np.mean(np.diff(t[cross]))*2 if cross.size >= 3 else float('nan')
    return ('alternation' if cross.size >= 2 else 'transient'), amp, per, f.max(), f.min()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--n', type=float, default=3.0); p.add_argument('--beta', type=float, default=0.4)
    p.add_argument('--Th', type=float, default=10.0); p.add_argument('--TD', type=float, default=0.3)
    p.add_argument('--R', type=float, default=10.0); p.add_argument('--T', type=float, default=None)
    p.add_argument('--scan', action='store_true', help='regime table over N_loc and Th for n = 2, 3, 5')
    p.add_argument('--csv', default=None, help='write the f_1(t) series of a single run to this file')
    p.add_argument('--scan-csv', default=None, help='with --scan, also write the regime table to this file')
    a = p.parse_args()
    if a.scan:
        cases = list(itertools.product([2., 3., 5.], [1.2, 1.5, 2., 3., 4., 6.], [3., 10., 30., 100.]))
        n, Nloc, Th = (np.array(c) for c in zip(*cases)); TD, R = 0.3, 10.0; beta = 2*Nloc/(n*R*TD)
        t, f = integrate(n, beta, Th, TD, R, T=30*Th.max(), dt=0.003)
        print('N_loc = n beta R TD/2 with TD = %g, R = %g; columns Th = 3, 10, 30, 100; A = alternation (period/Th)' % (TD, R))
        if a.scan_csv:
            with open(a.scan_csv, 'w') as fh:
                fh.write('n,N_loc,beta,Th,TD,R,regime,period_over_Th,f1_min,f1_max\n')
                for j, (nn, c, th) in enumerate(cases):
                    reg, amp, per, fmax, fmin = classify(t, f[:, j], th)
                    fh.write('%g,%g,%.6g,%g,%g,%g,%s,%.4g,%.4f,%.4f\n' % (nn, c, beta[j], th, TD, R, reg, per/th, fmin, fmax))
        for nn in [2., 3., 5.]:
            print('n = %g' % nn)
            for c in [1.2, 1.5, 2., 3., 4., 6.]:
                row = []
                for th in [3., 10., 30., 100.]:
                    j = cases.index((nn, c, th)); reg, amp, per, fmax, fmin = classify(t, f[:, j], th)
                    row.append({'sharing': '   S    ', 'lock-in': '   L    '}.get(reg, 'A(%5.2f)' % (per/th)))
                print('  N_loc %4.1f: %s' % (c, ' '.join(row)))
        return
    t, f = integrate(a.n, a.beta, a.Th, a.TD, a.R, T=a.T)
    reg, amp, per, fmax, fmin = classify(t, f[:, 0], a.Th)
    print('N_loc = %.2f, upper bound (n/2) ln(1+Th) = %.2f' % (a.n*a.beta*a.R*a.TD/2, a.n/2*np.log(1 + a.Th)))
    print('regime: %s; swing of f_1: %.2f .. %.2f; period: %.3g (%.2f Th)' % (reg, fmin, fmax, per, per/a.Th))
    if a.csv:
        np.savetxt(a.csv, np.column_stack([t, f[:, 0]]), header='t f1', comments='')


if __name__ == '__main__':
    main()
