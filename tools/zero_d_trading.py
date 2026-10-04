#!/usr/bin/env python3
"""Two shear-zone roots in series: do they share, lock in, or alternate?

Both roots carry the same far-field stress and their slip rates add to the plate rate, so the
partition is f_1 = S_1^-n / (S_1^-n + S_2^-n). Each root has a slow strengthening memory h (grows
with strain, relaxes in time) and a fast weakening feedback D (grows with strain, heals in time):

    S = (1 + h) exp(-beta D),   dh/dt = f - h/Th,   dD/dt = R f - D/TD

in time units of w*gamma_h/vL (w shear-zone width, gamma_h hardening strain). Th and TD are the
recovery and healing times in those units, R = gamma_h/gamma_D. Equal sharing is unstable when
N_loc = n beta R TD/2 > 1; lock-in is permanent when N_loc > (n/2) ln(1 + Th), roughly.

--laws implemented integrates the laws of SCycle's hardening and fabric states instead
(hardeningState.cpp, fabricState.cpp, PowerLaw::applyStrengthFactors), which saturate:

    root i creeps at [(1 + beta Phi_i) / H_i]^n times the common rate,  H_i = 1 + a (S_i - Sref)
    dS_i/dt   = Xh f_i (1 - S_i) - S_i                  (time in units of tau_r)
    dPhi_i/dt = (Xf f_i (1 - Phi_i) - Phi_i) / r

with Xh = e_full tau_r/gamma_h, Xf = e_full tau_c/gamma_f, r = tau_c/tau_r and e_full = vL/w the
strain rate of a root that takes the whole plate rate. The fabric's gain at equal sharing is
N_f = n beta u/((1 + u)(1 + u + beta u)), u = Xf/2, at most n beta/(1 + sqrt(1 + beta))^2: below 1
for every Xf unless beta > 3 (n = 3). Above 1 the fabric locks one root in and withstands a
hardening bias b = n ln(H_fast/H_slow) up to b_c; the memory alternates the roots when the bias it
builds there exceeds b_c, given a fast fabric (r of 0.03 or less).

  python3 tools/zero_d_trading.py --n 3 --beta 0.4 --Th 10 --TD 0.3 --R 10
  python3 tools/zero_d_trading.py --scan      # the regime table of docs/REVERSIBLE_STRENGTH_PLAN.md
  python3 tools/zero_d_trading.py --laws implemented --beta 8 --Xf 0.5 --a 0.7 --Xh 2.5 --r 0.03
  python3 tools/zero_d_trading.py --laws implemented --input data/stage5_batch2/alt8.in --efull 5e-13
  python3 tools/zero_d_trading.py --laws implemented --scan    # the table of section 2.2
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


def integrate_implemented(n, a, Sref, Xh, beta, Xf, r, T=40.0, eps=1e-3):
    """The laws as implemented (module docstring), RK4 in time units of tau_r from the symmetric
    steady state with a relative perturbation eps of Phi. a = 0 is no memory, beta = 0 no fabric.
    Returns t, f_1."""
    P = [np.atleast_1d(np.asarray(v, float)) for v in (n, a, Sref, Xh, beta, Xf, r)]
    M = max(p.size for p in P)
    n, a, Sref, Xh, beta, Xf, r = (np.broadcast_to(p, M).copy() for p in P)
    S0 = 0.5*Xh/(1 + 0.5*Xh); P0 = 0.5*Xf/(1 + 0.5*Xf)
    S = np.stack([S0, S0]); Ph = np.stack([P0*(1 + eps), P0*(1 - eps)])

    def rhs(S, Ph):
        lw = n*(np.log1p(beta*Ph) - np.log(1 + a*(S - Sref))); w = np.exp(lw - lw.max(axis=0))
        f = w/w.sum(axis=0)
        return Xh*f*(1 - S) - S, (Xf*f*(1 - Ph) - Ph)/r, f[0]

    dt = min(0.02*np.min(r/(1 + Xf)), 0.02/np.max(1 + Xh))
    N = int(T/dt); every = max(1, N//20000); ts, fs = [], []
    for k in range(N):
        k1s, k1p, f1 = rhs(S, Ph); k2s, k2p, _ = rhs(S + 0.5*dt*k1s, Ph + 0.5*dt*k1p)
        k3s, k3p, _ = rhs(S + 0.5*dt*k2s, Ph + 0.5*dt*k2p); k4s, k4p, _ = rhs(S + dt*k3s, Ph + dt*k3p)
        S = S + dt*(k1s + 2*k2s + 2*k3s + k4s)/6; Ph = Ph + dt*(k1p + 2*k2p + 2*k3p + k4p)/6
        if k % every == 0:
            ts.append(k*dt); fs.append(f1.copy())
    return np.array(ts), np.array(fs)


def fabric_gain(n, beta, Xf):
    """N_f, the fabric's localization gain at equal sharing (equal sharing is unstable above 1)."""
    u = 0.5*Xf
    return n*beta*u/((1 + u)*(1 + u + beta*u))


def fabric_hold(n, beta, Xf):
    """The largest hardening bias b_c that the fabric's lock-in branch withstands, the share f* of
    the fast root where the branch ends, and the lock-in share f_L without a memory; b_c = 0 and
    f* = 0.5 when equal sharing is stable."""
    f = np.linspace(0.5, 1 - 1e-6, 200001); Phi = lambda g: Xf*g/(1 + Xf*g)
    b = n*np.log((1 + beta*Phi(f))/(1 + beta*Phi(1 - f))) - np.log(f/(1 - f))
    j = int(np.argmax(b)); down = np.flatnonzero((b[:-1] > 0) & (b[1:] <= 0))
    return b[j], f[j], (f[down[-1]] if down.size else float('nan'))


def memory_push(n, a, Sref, Xh, f):
    """The bias b = n ln(H_fast/H_slow) of the hardening memory between roots at shares f and 1 - f,
    each in its steady state S = Xh f/(1 + Xh f)."""
    Sss = lambda g: Xh*g/(1 + Xh*g)
    return n*np.log((1 + a*(Sss(f) - Sref))/(1 + a*(Sss(1 - f) - Sref)))


def keys_from_input(path, efull):
    """n, a, Sref, Xh, beta, Xf, r and tau_r (s) from the hardening and fabric keys of an SCycle
    input file (first values of the depth profiles), with e_full (1/s) the strain rate of a root
    that takes the whole plate rate. Without hardening a = 0 and tau_r = 100 tau_c (time unit only);
    without fabric beta = 0 and tau_c = tau_r/100; without either the time unit is 1 s."""
    k = {}
    for line in open(path):
        if ' = ' in line:
            key, val = line.split(' = ', 1); val = val.split('#')[0].strip()
            k[key.strip()] = val.strip('[]').split()[0] if val.startswith('[') else val.split()[0]
    hard = k.get('hard_type', 'off') != 'off'; fab = k.get('fabric_type', 'off') != 'off'
    n = float(k.get('disl_nVals', 3))
    tauC = float(k['fabric_tauC0']) if fab else None
    tauR = float(k['hard_tauR0']) if hard else (100*tauC if fab else 1.0)
    tauC = tauC or tauR/100
    if hard and 'hard_SRef' not in k: raise SystemExit('hard_eRef: give hard_SRef to screen this input')
    a, Sref, Xh = (float(k.get('hard_a', 0.5)), float(k['hard_SRef']), efull*tauR/float(k['hard_gammaHVals'])) if hard else (0.0, 0.0, 0.0)
    beta, Xf = (float(k.get('fabric_betaF', 4.0)), efull*tauC/float(k['fabric_gammaFVals'])) if fab else (0.0, 0.0)
    return n, a, Sref, Xh, beta, Xf, tauC/tauR, tauR


def classify(t, f, Th):
    """Regime from the second half of the series: sharing, lock-in or alternation (with period)."""
    half = t > 0.5*t[-1]; f = f[half]; t = t[half]
    amp = f.max() - f.min(); cross = np.flatnonzero(np.diff(np.sign(f - 0.5)) != 0)
    if amp < 0.05 and abs(f[-1] - 0.5) < 0.05: return 'sharing', amp, float('nan'), f.max(), f.min()
    if amp < 0.05: return 'lock-in', amp, float('nan'), f.max(), f.min()
    per = np.mean(np.diff(t[cross]))*2 if cross.size >= 3 else float('nan')
    return ('alternation' if cross.size >= 2 else 'transient'), amp, per, f.max(), f.min()


def main_implemented(a):
    if a.scan:
        cases = []
        for beta in [4., 5., 6., 8., 10.]:
            for Xf in [0.5, round(2/np.sqrt(1 + beta), 3), 1.0]:
                cases += [(beta, Xf, aa, Xh, r) for aa in [0.5, 0.7, 1.0] for Xh in [1.5, 2.5, 4.0] for r in [0.01, 0.03, 0.1]]
        B, XF, A, XH, R = (np.array(c) for c in zip(*cases))
        t, f = integrate_implemented(a.n, A, a.Sref, XH, B, XF, R)
        out = open(a.scan_csv, 'w') if a.scan_csv else None
        if out: out.write('n,beta,Xf,a,Xh,r,Sref,N_f,b_c,f_star,b_max,regime,period_over_tauR,f1_min,f1_max\n')
        print('n = %g, Sref = %g; per row the regime for r = tau_c/tau_r = 0.01, 0.03, 0.1 (S sharing, L lock-in,' % (a.n, a.Sref))
        print('A alternation with the period in tau_r and the swing of f_1)')
        print('beta    Xf    a   Xh   N_f    b_c  b_max | r = 0.01            r = 0.03            r = 0.1')
        for j0 in range(0, len(cases), 3):
            beta, Xf, aa, Xh, _ = cases[j0]
            bc, fs, fL = fabric_hold(a.n, beta, Xf); bm = memory_push(a.n, aa, a.Sref, Xh, fs)
            cells = []
            for j in range(j0, j0 + 3):
                reg, amp, per, fmax, fmin = classify(t, f[:, j], None)
                cells.append('%s %5s %.2f-%.2f' % (reg[0].upper(), '' if np.isnan(per) else '%.2f' % per, fmin, fmax))
                if out: out.write('%g,%g,%g,%g,%g,%g,%g,%.4f,%.4f,%.4f,%.4f,%s,%.4g,%.4f,%.4f\n' % (a.n, beta, Xf, aa, Xh, cases[j][4], a.Sref,
                                  fabric_gain(a.n, beta, Xf), bc, fs, bm, reg, per, fmin, fmax))
            print('%4g %5.3g %4g %4g %5.3f %6.3f %6.3f | %s' % (beta, Xf, aa, Xh, fabric_gain(a.n, beta, Xf), bc, bm, ' | '.join(cells)))
        return
    if a.input:
        if a.efull is None: raise SystemExit('--input needs --efull (vL/w, 1/s)')
        n, aa, Sref, Xh, beta, Xf, r, tauR = keys_from_input(a.input, a.efull)
    else:
        n, aa, Sref, Xh, beta, Xf, r, tauR = a.n, a.a, a.Sref, a.Xh, a.beta, a.Xf, a.r, None
    print('n = %g; memory a = %g, Sref = %g, Xh = %.4g; fabric beta = %g, Xf = %.4g; r = tau_c/tau_r = %.4g'
          % (n, aa, Sref, Xh, beta, Xf, r))
    bc, fs, fL = fabric_hold(n, beta, Xf); Nf = fabric_gain(n, beta, Xf)
    print('fabric gain N_f = %.3f (largest over Xf: %.3f); ' % (Nf, n*beta/(1 + np.sqrt(1 + beta))**2), end='')
    if Nf <= 1: print('equal sharing is stable, the screen predicts sharing')
    else:
        bm = memory_push(n, aa, Sref, Xh, fs) if aa > 0 else 0.0
        print('lock-in at f = %.3f without memory, hold b_c = %.3f (branch ends at f = %.3f);' % (fL, bc, fs))
        print('memory bias there b_max = %.3f: the screen predicts %s' % (bm, 'alternation' if bm > bc else 'lock-in'))
    t, f = integrate_implemented(n, aa, Sref, Xh, beta, Xf, r, T=a.T or 40.0)
    reg, amp, per, fmax, fmin = classify(t, f[:, 0], None)
    yrs = '' if tauR is None or np.isnan(per) else ' = %.0f yr' % (per*tauR/3.15576e7)
    print('two-root ODE: %s; swing of f_1: %.2f .. %.2f; period: %.3g tau_r%s' % (reg, fmin, fmax, per, yrs))
    if a.csv:
        np.savetxt(a.csv, np.column_stack([t, f[:, 0]]), header='t_over_tauR f1', comments='')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--laws', choices=['linear', 'implemented'], default='linear',
                   help='linear: the model of section 2; implemented: the saturating laws of the code')
    p.add_argument('--n', type=float, default=3.0); p.add_argument('--beta', type=float, default=0.4)
    p.add_argument('--Th', type=float, default=10.0); p.add_argument('--TD', type=float, default=0.3)
    p.add_argument('--R', type=float, default=10.0); p.add_argument('--T', type=float, default=None)
    p.add_argument('--a', type=float, default=0.5, help='implemented: hardening amplitude (hard_a; 0: no memory)')
    p.add_argument('--Sref', type=float, default=0.0, help='implemented: hard_SRef')
    p.add_argument('--Xh', type=float, default=2.5, help='implemented: e_full tau_r/gamma_h')
    p.add_argument('--Xf', type=float, default=0.5, help='implemented: e_full tau_c/gamma_f')
    p.add_argument('--r', type=float, default=0.03, help='implemented: tau_c/tau_r')
    p.add_argument('--input', default=None, help='implemented: take the parameters from the hard_ and fabric_ keys of this input')
    p.add_argument('--efull', type=float, default=None, help='with --input: e_full = vL/w (1/s)')
    p.add_argument('--scan', action='store_true', help='regime table (linear: over N_loc and Th for n = 2, 3, 5)')
    p.add_argument('--csv', default=None, help='write the f_1(t) series of a single run to this file')
    p.add_argument('--scan-csv', default=None, help='with --scan, also write the regime table to this file')
    a = p.parse_args()
    if a.laws == 'implemented':
        return main_implemented(a)
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
