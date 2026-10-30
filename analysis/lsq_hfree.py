#!/usr/bin/env python3
"""
lsq_hfree.py — closed-form (s2, s4) recovery from VELOCITIES ONLY.
h is never observed and never used, anywhere.

WHY THIS EXISTS
---------------
The `--velocity-only` flag in lsq_from_data.py drops the continuity equation's
rows but still uses the true h inside eq1/eq2's pressure terms, so it only ever
proved that eq5 is REDUNDANT — not that h is dispensable. On solar data h is not
observable at all, so that gap had to be closed before any solar claim.

It is closed analytically. For a single azimuthal mode, d/dphi -> i*m makes h
appear ALGEBRAICALLY in eq2 (verified: eq2 contains no derivative of h), so eq2
inverts for h in terms of the velocities. Substituting into eq1 gives a relation
in the velocities alone.

WHAT THE SYMBOLIC DERIVATION ESTABLISHED (sympy 2026-08-13; nu4 added 2026-09-30)
--------------------------------------------------------------------------------
  * h itself is EXACTLY LINEAR in (s2, s4)  -- monomials (0,0),(0,1),(1,0)
  * h-free eq1 is EXACTLY LINEAR in (s2, s4) -- so the closed-form least squares
    survives elimination:   R = c0 + s2*c2 + s4*c4, velocities only.
  * h-free eq5 is QUADRATIC (it contains s2^2, s2*s4, s4^2), because eq5 holds
    Omega_d*h and h already carries s2/s4. eq5 therefore CANNOT be used in the
    closed form and is not implemented here.
  * G CANCELS EXACTLY. eq2 gives h ~ 1/G; eq1 multiplies it back by G. The
    h-free relation therefore contains NO G -- the recovery needs no knowledge
    of the layer depth H0. (Note hfree_coeffs() takes no G argument.) This
    survives the addition of hyperdiffusion below.

HYPERDIFFUSION (the nu4 terms) — added 2026-09-30
-------------------------------------------------
sw_model.py integrates eq1 and eq2 WITH a mu-biharmonic term:

    eq1:  d_t u_th + Om*im*u_th - f_th*u_ph - G*s_th*d_mu h + nu4*d_mu^4 u_th = 0
    eq2:  d_t u_ph + Om*im*u_ph + f_ph*u_th + G*varsigma*im*h + nu4*d_mu^4 u_ph = 0

The original elimination dropped both nu4 terms. That omission is a MODEL
MISMATCH, not a numerical error: it biases (s2, s4) by a fixed amount that is
independent of the Chebyshev order and of the noise level -- which is exactly
the signature measured on 2026-09-30 (1.85% / 8.90%, flat across cheb=20..60
and unchanged by rank-1 SVD truncation even at noise 0.005).

Re-deriving with nu4 retained (verified symbolically by --selftest):
  * still EXACTLY LINEAR in (s2, s4)
  * G still cancels
  * nu4 enters c0 ONLY -- c2 and c4 are untouched, so cond() does not change
  * the correction has the same functional form as the (Ut, Vt, Vmt) group with
    (d_mu^4 u_th, d_mu^4 u_ph, d_mu^5 u_ph) substituted, because nu4*d_mu^4 sits
    in eq1/eq2 in precisely the slot d_t occupies:

      c0 += nu4*(m*U4 + i*mu*V4 - i*V5 + i*mu^2*V5 - i*eps*sqrt(1-mu^2)*V5)/m

THE COST
--------
Inviscid h-free eq1 needs a MIXED SECOND derivative d2(u_phi)/dmu dt. With nu4
retained it additionally needs d_mu^5 u_phi -- a FIFTH mu-derivative. These are
spectral (truncated-Chebyshev) derivatives, not finite differences, so they are
computable on smooth data, but they amplify noise hard. Hence --no-visc: the Sun
has no biharmonic hyperdiffusion (it is a numerical device), so the inviscid
coefficients are the physically relevant ones there, and the nu4 terms exist to
remove the bias when validating against THIS simulation.

    --no-visc  -> inviscid coefficients (original behaviour; the solar limit)
    default    -> nu4 from input.json (removes the simulation's own bias)

Running both on the same data measures the hyperdiffusion bias directly.
nu4_h does NOT appear: only eq1 and eq2 enter the elimination, and nu4_h lives
in eq5.

USAGE
-----
    # on clean simulation snapshots (start here -- truth is known)
    python3 lsq_hfree.py --input <data-dir>/input.json --tmin 10 --scan

    # same, with the nu4 terms dropped -- the difference IS the nu4 bias
    python3 lsq_hfree.py --input <data-dir>/input.json --tmin 10 --scan --no-visc

    # on a PINN checkpoint (PINN-smoothed field, then h-free recovery)
    python3 lsq_hfree.py --ckpt checkpoints --input <data-dir>/input.json --scan

    # verify the hand-transcribed coefficients against sympy (needs sympy)
    python3 lsq_hfree.py --selftest

Flags: --cheb / --svd for a single setting, --scan to sweep them,
       --tmin/--tmax for the window, --field pred|true with --ckpt,
       --nu4 to override, --no-visc to drop the hyperdiffusion terms.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import h5py

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--input", default=None, help="input.json (m, G, s2, s4, eps_reg, nu4)")
ap.add_argument("--ckpt", default=None,
                help="a PINN checkpoint dir holding pinn_fields.h5; if given, the "
                     "field comes from there instead of raw snapshots")
ap.add_argument("--field", choices=["pred", "true"], default="pred",
                help="with --ckpt: PINN field (pred) or the field it was given (true)")
ap.add_argument("--tmin", type=float, default=-np.inf)
ap.add_argument("--tmax", type=float, default=np.inf)
ap.add_argument("--edge-trim", type=int, default=None)
ap.add_argument("--cheb", type=int, default=30)
ap.add_argument("--svd", type=int, default=0)
ap.add_argument("--scan", action="store_true", help="sweep cheb x svd")
ap.add_argument("--nu4", type=float, default=None,
                help="override the hyperdiffusion coefficient (default: from input.json)")
ap.add_argument("--no-visc", action="store_true",
                help="drop the nu4 terms from the h-free relation. This is the "
                     "SOLAR limit (no biharmonic hyperdiffusion in the Sun) and "
                     "reproduces the original inviscid coefficients.")
ap.add_argument("--band-deg", type=float, default=None,
                help="use ONLY the latitude band |lat| <= BAND (deg): the mu "
                     "derivatives come from a least-squares Chebyshev fit on that "
                     "band alone, and the equation is applied only there. Needed "
                     "whenever the data stop short of the poles (Stage 3 output, "
                     "ring-diagram flows): zero-filled caps put a jump at the band "
                     "edge that the full-range derivatives turn into a 17%% s4 bias")
ap.add_argument("--band-trim", type=int, default=4,
                help="with --band-deg: drop this many mu points at each band edge "
                     "from the equation rows (the fit is least accurate there)")
ap.add_argument("--selftest", action="store_true",
                help="check the numpy coefficients against a fresh sympy "
                     "derivation on random data, then exit")


# ── the h-free eq1 coefficients ──────────────────────────────────────────────
# R = c0 + s2*c2 + s4*c4, from eliminating h between eq1 and eq2.
# Generated from sympy rather than hand-copied; --selftest re-derives
# symbolically and compares, so a transcription slip cannot go unnoticed.
#   U=u_theta  V=u_phi  Um=d_mu u_theta  Ut=d_t u_theta
#   Vm=d_mu u_phi  Vt=d_t u_phi  Vmt=d_mu d_t u_phi
#   U4=d_mu^4 u_theta  V4=d_mu^4 u_phi  V5=d_mu^5 u_phi   (nu4 terms only)
# nu4 enters c0 ONLY. c2 and c4 are identical to the inviscid derivation.
def hfree_coeffs(mu, m, eps, nu4, U, V, Um, Ut, Vm, Vt, Vmt, U4, V4, V5):
    sq = np.sqrt(np.maximum(1 - mu**2, 0.0))
    j = 1j
    c0 = (-2*j*U*eps*sq + 4*j*U*mu**2 - 2*j*U
          - 2*j*Um*eps*mu*sq + 2*j*Um*mu**3 - 2*j*Um*mu
          - j*Vmt*eps*sq + j*Vmt*mu**2 - j*Vmt
          + j*Vt*mu + m*(Ut - 2*V*mu)) / m
    # hyperdiffusion: same shape as the (Ut, Vt, Vmt) group above, because
    # nu4*d_mu^4 occupies exactly the slot d_t occupies in eq1/eq2.
    c0 = c0 + nu4*(-j*V5*eps*sq + j*V5*mu**2 - j*V5
                   + j*V4*mu + m*U4) / m
    c2 = (12*j*U*eps*mu**2*sq - 2*j*U*eps*sq
          - 16*j*U*mu**4 + 16*j*U*mu**2 - 2*j*U
          + 4*j*Um*eps*mu**3*sq - 2*j*Um*eps*mu*sq
          - 4*j*Um*mu**5 + 6*j*Um*mu**3 - 2*j*Um*mu
          - m*mu*(j*U*m*mu + 2*V*eps*sq - 5*V*mu**2 + 2*V
                  + Vm*eps*mu*sq - Vm*mu**3 + Vm*mu)) / m
    c4 = mu**2*(30*j*U*eps*mu**2*sq - 12*j*U*eps*sq
                - 36*j*U*mu**4 + 46*j*U*mu**2 - 12*j*U
                + 6*j*Um*eps*mu**3*sq - 4*j*Um*eps*mu*sq
                - 6*j*Um*mu**5 + 10*j*Um*mu**3 - 4*j*Um*mu
                - m*mu*(j*U*m*mu + 4*V*eps*sq - 7*V*mu**2 + 4*V
                        + Vm*eps*mu*sq - Vm*mu**3 + Vm*mu)) / m
    return c0, c2, c4


def selftest():
    """Re-derive symbolically and compare to the numpy transcription."""
    import sympy as sp
    mu_s, t_s, m_s, G_s, eps_s, s2_s, s4_s, nu4_s = sp.symbols(
        'mu t m G eps s2 s4 nu4', real=True)
    Uth = sp.Function('Uth')(mu_s, t_s)
    Uph = sp.Function('Uph')(mu_s, t_s)
    Hh = sp.Function('H')(mu_s, t_s)
    Om = -s2_s*mu_s**2 - s4_s*mu_s**4
    st = sp.sqrt(1-mu_s**2)
    ist = 1/(st+eps_s)
    corf_th = 2*(1+Om)*mu_s
    corf_phi = 2*(1+Om)*mu_s + (1-mu_s**2)*(2*s2_s*mu_s + 4*s4_s*mu_s**3)
    # eq1 / eq2 EXACTLY as sw_model.py integrates them, hyperdiffusion included.
    eq1 = (sp.diff(Uth, t_s) + Om*sp.I*m_s*Uth - corf_th*Uph
           - G_s*st*sp.diff(Hh, mu_s) + nu4_s*sp.diff(Uth, mu_s, 4))
    eq2 = (sp.diff(Uph, t_s) + Om*sp.I*m_s*Uph + corf_phi*Uth
           + G_s*ist*sp.I*m_s*Hh + nu4_s*sp.diff(Uph, mu_s, 4))
    Hsol = sp.solve(eq2, Hh)[0]
    e = sp.expand(eq1.subs(sp.Derivative(Hh, mu_s), sp.diff(Hsol, mu_s)).subs(Hh, Hsol))
    P = sp.Poly(e, s2_s, s4_s)
    mons = sorted(P.monoms())
    print("h-free eq1 monomials in (s2,s4):", mons)
    assert all(sum(k) <= 1 for k in mons), "NOT LINEAR — derivation broken"
    print("  -> LINEAR confirmed")
    assert G_s not in e.free_symbols, "G did not cancel — derivation broken"
    print("  -> G CANCELS confirmed (recovery needs no layer depth H0)")

    U_, V_, Um_, Ut_, Vm_, Vt_, Vmt_, U4_, V4_, V5_ = sp.symbols(
        'U V Um Ut Vm Vt Vmt U4 V4 V5')
    # highest derivatives first so nothing is shadowed
    sub = [(sp.Derivative(Uph, (mu_s, 5)), V5_),
           (sp.Derivative(Uph, (mu_s, 4)), V4_),
           (sp.Derivative(Uth, (mu_s, 4)), U4_),
           (sp.Derivative(Uph, mu_s, t_s), Vmt_),
           (sp.Derivative(Uph, t_s, mu_s), Vmt_),
           (sp.Derivative(Uth, mu_s), Um_), (sp.Derivative(Uth, t_s), Ut_),
           (sp.Derivative(Uph, mu_s), Vm_), (sp.Derivative(Uph, t_s), Vt_),
           (Uth, U_), (Uph, V_)]
    syms = [mu_s, m_s, eps_s, nu4_s,
            U_, V_, Um_, Ut_, Vm_, Vt_, Vmt_, U4_, V4_, V5_]
    fns = []
    for k in (sp.Integer(1), s2_s, s4_s):
        c = P.coeff_monomial(k)
        for old, new in sub:
            c = c.subs(old, new)
        leftover = c.atoms(sp.Derivative)
        assert not leftover, f"unsubstituted derivative(s): {leftover}"
        fns.append(sp.lambdify(syms, c, "numpy"))

    rng = np.random.default_rng(0)
    # nu4 deliberately O(1) here so the hyperdiffusion terms are exercised
    vals = [rng.uniform(-0.9, 0.9), 5.0, 0.05, 0.37] + \
           [rng.standard_normal() + 1j*rng.standard_normal() for _ in range(10)]
    sym = [f(*vals) for f in fns]
    num = hfree_coeffs(*vals)
    ok = True
    for nm, a_, b_ in zip(("c0", "c2", "c4"), sym, num):
        d = abs(a_-b_)/max(abs(a_), 1e-30)
        print(f"  {nm}: sympy={a_:.6e}  numpy={b_:.6e}  rel.diff={d:.2e}")
        ok &= d < 1e-10
    # and check that nu4=0 reproduces the inviscid coefficients exactly
    v0 = list(vals); v0[3] = 0.0
    ok0 = True
    for a_, b_ in zip([f(*v0) for f in fns], hfree_coeffs(*v0)):
        ok0 &= abs(a_-b_)/max(abs(a_), 1e-30) < 1e-10
    print(f"  nu4=0 limit reproduces inviscid coefficients: {'yes' if ok0 else 'NO'}")
    ok &= ok0
    print("SELFTEST", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


a = ap.parse_args()
if a.selftest:
    sys.exit(selftest())
if not a.input:
    sys.exit("--input is required (needs m, G, s2, s4, eps_reg, nu4)")

P = json.load(open(a.input))
m = int(P["m"])
S2_TRUE, S4_TRUE = float(P["s2"]), float(P["s4"])
EPS = float(P.get("eps_reg", 0.05))
EDGE = a.edge_trim if a.edge_trim is not None else int(P.get("edge_trim", 8))

# ── hyperdiffusion coefficient ───────────────────────────────────────────────
if a.no_visc:
    NU4 = 0.0
    VISC_NOTE = "DROPPED (--no-visc) — this is the solar limit"
elif a.nu4 is not None:
    NU4 = float(a.nu4)
    VISC_NOTE = f"nu4={NU4:.3e} (from --nu4)"
elif "nu4" in P:
    NU4 = float(P["nu4"])
    VISC_NOTE = f"nu4={NU4:.3e} (from input.json)"
else:
    NU4 = 0.0
    VISC_NOTE = ("nu4 ABSENT from input.json -> treated as 0. If the run used "
                 "hyperdiffusion, pass --nu4 explicitly.")

# ── load the field ───────────────────────────────────────────────────────────
if a.ckpt:
    fp = os.path.join(a.ckpt, "pinn_fields.h5")
    with h5py.File(fp, "r") as f:
        mu = f["mu"][:]; t = f["t"][:]
        k = a.field
        Uth = f["u_theta"][f"{k}_r"][:] + 1j*f["u_theta"][f"{k}_i"][:]
        Uph = f["u_phi"][f"{k}_r"][:] + 1j*f["u_phi"][f"{k}_i"][:]
    src = f"{fp} [{a.field}]"
else:
    d = os.path.dirname(os.path.abspath(a.input))
    sd = P.get("snapshots_dir", "snapshots")
    sd = sd if os.path.isabs(sd) else os.path.normpath(os.path.join(d, sd))
    files = sorted(glob.glob(os.path.join(sd, "snapshots_s*.h5")) or
                   glob.glob(os.path.join(sd, "*.h5")))
    if not files:
        sys.exit(f"no snapshots in {sd}")
    tc, Uc, Vc = [], [], []
    for fpath in files:
        with h5py.File(fpath, "r") as f:
            tt = f["scales/sim_time"][:]
            for nm in f["scales"]:
                if "mu" in nm.lower():
                    arr = np.array(f["scales"][nm][:])
                    if arr.ndim == 1:
                        mu = arr.ravel()
            nphi = f["tasks/u_theta"].shape[1]
            n = min(len(tt), f["tasks/u_theta"].shape[0])
            msk = (tt[:n] >= a.tmin) & (tt[:n] <= a.tmax)
            if not msk.any():
                continue
            mm = lambda nm: np.fft.fft(f[f"tasks/{nm}"][:n], axis=1)[:, m, :]/nphi
            tc.append(tt[:n][msk]); Uc.append(mm("u_theta")[msk]); Vc.append(mm("u_phi")[msk])
    t = np.concatenate(tc); Uth = np.concatenate(Uc); Uph = np.concatenate(Vc)
    o = np.argsort(t); t, Uth, Uph = t[o], Uth[o], Uph[o]
    src = sd

if a.band_deg is not None:
    # keep only the observed band; everything else is never touched
    MU_B = float(np.sin(np.deg2rad(a.band_deg)))
    keep = np.abs(mu) <= MU_B + 1e-12
    mu = mu[keep]; Uth = Uth[:, keep]; Uph = Uph[:, keep]
N_t, N_mu = len(t), len(mu)
mus = np.sort(mu)
if a.band_deg is not None:
    EDGE = a.band_trim
MU_MAX = float(mus[-EDGE-1])
order = np.argsort(mu); inv = np.argsort(order); msrt = mu[order]
dt_u = float(np.median(np.diff(t)))

print("="*74)
print(f"lsq_hfree  |  {src}")
print(f"  m={m}  true s2={S2_TRUE}  s4={S4_TRUE}  eps_reg={EPS}  edge_trim={EDGE}")
print(f"  N_t={N_t}  N_mu={N_mu}  t in [{t.min():.1f}, {t.max():.1f}]")
print(f"  hyperdiffusion: {VISC_NOTE}")
if a.band_deg is not None:
    print(f"  BAND: |lat| <= {a.band_deg:g} deg (|mu| <= {MU_B:.4f}), {N_mu} mu points; "
          f"band-limited Chebyshev fit, {a.band_trim} edge points trimmed")
print("  VELOCITIES ONLY — h is never observed and never used")
print("  G-FREE — the elimination cancels G, so no layer depth H0 is needed")
print("  (one linear equation family: h-free eq1. eq5 is quadratic after")
print("   elimination and cannot be used in closed form.)")
if NU4:
    print("  NOTE: the nu4 terms require d_mu^5 u_phi. Expect faster degradation")
    print("        at high cheb on noisy data — compare against --no-visc.")
print("="*74)


def cheb_D(x, K):
    x = np.asarray(x, float); N = len(x)
    th = np.arccos(np.clip(x, -1, 1)); kk = np.arange(N)
    C = (2.0/N)*np.cos(np.outer(kk, th)); C[0] *= 0.5
    Dm = np.zeros((N, N))
    for k in range(N-1):
        Dm[k, k+1] = 2*(k+1)
    for k in range(N-3, -1, -1):
        Dm[k, :] += Dm[k+2, :]
    Dm[0, :] *= 0.5
    Tr = np.diag((np.arange(N) < K).astype(float))
    return np.cos(np.outer(th, kk)) @ Dm @ Tr @ C


def cheb_D_band(x, K, xb):
    """d/dmu on an ARBITRARY node set inside [-xb, xb]: least-squares fit of
    Chebyshev T_0..T_{K-1} in y = x/xb, differentiated exactly."""
    y = np.asarray(x, float)/xb
    K = min(K, len(y))
    V = np.polynomial.chebyshev.chebvander(y, K - 1)          # (N, K)
    dV = np.zeros_like(V)
    for k in range(K):
        e = np.zeros(K); e[k] = 1.0
        dV[:, k] = np.polynomial.chebyshev.chebval(y, np.polynomial.chebyshev.chebder(e))
    return (dV @ np.linalg.pinv(V))/xb


def svd_trunc(A, r):
    if r is None or r <= 0 or r >= min(A.shape):
        return A
    U_, S_, Vh_ = np.linalg.svd(A, full_matrices=False)
    return (U_[:, :r]*S_[:r]) @ Vh_[:r, :]


def dt4(A):
    d = np.zeros_like(A)
    d[2:-2] = (-A[4:] + 8*A[3:-1] - 8*A[1:-3] + A[:-4])/(12*dt_u)
    return d


def solve(K, r):
    U_, V_ = svd_trunc(Uth, r), svd_trunc(Uph, r)
    D = cheb_D(msrt, K) if a.band_deg is None else cheb_D_band(msrt, K, MU_B)
    dmu = lambda F: (F[:, order] @ D.T)[:, inv]
    Um, Ut = dmu(U_), dt4(U_)
    Vm, Vt = dmu(V_), dt4(V_)
    Vmt = dmu(dt4(V_))                      # the mixed second derivative
    if NU4:
        U4 = dmu(dmu(dmu(dmu(U_))))         # d_mu^4 u_theta
        V4 = dmu(dmu(dmu(dmu(V_))))         # d_mu^4 u_phi
        V5 = dmu(V4)                        # d_mu^5 u_phi  (the new cost)
    else:
        U4 = V4 = V5 = np.zeros_like(U_)
    mu2 = mu[None, :]
    c0, c2, c4 = hfree_coeffs(mu2, m, EPS, NU4, U_, V_, Um, Ut, Vm, Vt, Vmt,
                              U4, V4, V5)
    msk = np.zeros((N_t, N_mu), bool)
    msk[:, np.abs(mu) <= MU_MAX] = True
    msk[:2, :] = False; msk[-2:, :] = False   # dt4 stencil edges
    S = lambda x: np.concatenate([x[msk].real, x[msk].imag])
    b = S(c0)
    A = np.column_stack([S(c2), S(c4)])
    sol, *_ = np.linalg.lstsq(A, -b, rcond=None)
    cond = float(np.linalg.cond(A))
    res = float(np.linalg.norm(A@sol + b)/(np.linalg.norm(b)+1e-30))
    return float(sol[0]), float(sol[1]), cond, res


def report(K, r):
    s2, s4, cond, res = solve(K, r)
    e2 = abs(s2-S2_TRUE)/abs(S2_TRUE)*100
    e4 = abs(s4-S4_TRUE)/abs(S4_TRUE)*100
    print(f"  cheb={K:<3d} svd={r:<2d}  s2={s2:+.5f} ({e2:6.2f}%)  "
          f"s4={s4:+.5f} ({e4:7.2f}%)  cond={cond:7.2f}  resid={res:.3e}")
    return e2+e4, (K, r), (s2, s4, cond)


if a.scan:
    print("  setting        s2                 s4                 cond      resid")
    best = None
    for K in (20, 30, 40, 60):
        for r in (0, 1, 2, 4):
            try:
                sc, cfg, vals = report(K, r)
                if best is None or sc < best[0]:
                    best = (sc, cfg, vals)
            except Exception as ex:
                print(f"  cheb={K} svd={r}: failed ({ex})")
    print("="*74)
    print(f"BEST (truth-informed, for validation only): cheb={best[1][0]}, svd={best[1][1]}")
    print(f"   s2={best[2][0]:.6f}  s4={best[2][1]:.6f}  cond={best[2][2]:.2f}")
    print("   NOTE: 'best of scan' is a MINIMUM over many noisy estimates and")
    print("   flatters itself as noise rises — for a degradation curve use a")
    print("   FIXED setting or a truth-free selector instead.")
else:
    report(a.cheb, a.svd)
print("="*74)