#!/usr/bin/env python3
"""
fit_growth.py — extract exponential growth rate from spherical-2d-hydro logs.

KE ~ exp(2*sigma*t), so amplitude rate sigma = 0.5 * d ln(KE)/dt.

Two modes:
  manual window :  fit ln(KE) over [--tmin, --tmax]
  auto tail     :  --tail  detects the onset of exponential growth and fits
                   only the growing tail (ignores the neutral phase before it)

Usage:
  python3 fit_growth.py logger_info.txt --tmin 6 --tmax 36
  python3 fit_growth.py run_Nmu96.log --tail
  python3 fit_growth.py run_*.log --tail --log growth_scan.log
  python3 fit_growth.py run.log --check          # stable / growing / decaying verdict
"""
import argparse, re, sys
import numpy as np

LINE = re.compile(r"t=([-\d.eE+]+).*?KE=([-\d.eE+]+)")


def parse(path):
    t, ke = [], []
    with open(path) as f:
        for ln in f:
            m = LINE.search(ln)
            if m:
                t.append(float(m.group(1))); ke.append(float(m.group(2)))
    return np.array(t), np.array(ke)


def linfit(x, y):
    """least-squares line; return (slope, R^2)."""
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intc), res, *_ = np.linalg.lstsq(A, y, rcond=None)
    ss_res = res[0] if len(res) else float(np.sum((y - (slope*x + intc))**2))
    ss_tot = float(np.sum((y - y.mean())**2)) + 1e-300
    return float(slope), 1.0 - ss_res/ss_tot


def envelope_peaks(t, ke):
    """
    Local-maxima envelope of KE(t): indices where ke[i] is a strict local max
    (ke[i] > both neighbours). Works on the raw sampled series -- no
    scipy dependency. Requires at least 3 points around each candidate.
    """
    idx = [i for i in range(1, len(ke)-1) if ke[i] > ke[i-1] and ke[i] > ke[i+1]]
    return np.array(idx, dtype=int)


def stability_check(t, ke, emit, growing_tol, decaying_tol):
    """
    Decide GROWING / DECAYING / NEUTRAL (STABLE) from the envelope (peak-to-peak)
    trend of KE, which removes the oscillation and isolates the secular trend.
    sigma_env is defined the same way as the tail-fit: KE_envelope ~ exp(2*sigma*t).
    """
    pk = envelope_peaks(t, ke)
    if pk.size < 4:
        emit("   not enough oscillation peaks detected to fit an envelope "
             "(need >=4) -- try a longer run or check --tmin/--tmax.")
        return None
    tp, yp = t[pk], np.log(ke[pk])
    slope, r2 = linfit(tp, yp)
    sigma_env = 0.5 * slope

    if sigma_env > growing_tol:
        verdict = "GROWING (unstable)"
    elif sigma_env < -decaying_tol:
        verdict = "DECAYING (damped)"
    else:
        verdict = "NEUTRAL / STABLE (bounded oscillation)"

    emit(f"   envelope fit: {pk.size} peaks over t=[{tp[0]:.2f},{tp[-1]:.2f}]")
    emit(f"   sigma_envelope = {sigma_env:+.5f}   R2={r2:.4f}")
    emit(f"   thresholds: growing > {growing_tol:+.5f}, "
         f"decaying < {-decaying_tol:+.5f}")
    emit(f"   VERDICT: {verdict}")
    return sigma_env, verdict


def find_tail(t, y, r2min, min_efolds, min_pts):
    """
    y = ln(KE). Scan onset index from early to late; return the EARLIEST i such
    that the linear fit of y[i:] has slope>0, R^2>=r2min, and the tail spans at
    least min_efolds in ln(KE). That is the longest clean exponential tail.
    Returns (i, slope, r2) or None.
    """
    N = len(t)
    for i in range(0, N - min_pts):
        seg = y[i:]
        if (seg.max() - seg.min()) < min_efolds:   # not enough growth left
            continue
        slope, r2 = linfit(t[i:], seg)
        if slope > 0 and r2 >= r2min:
            return i, slope, r2
    return None


def analyse(path, args, emit):
    t, ke = parse(path)
    if t.size == 0:
        emit(f"{path}: no KE lines found"); return
    good = ke > 0
    t, ke = t[good], ke[good]
    y = np.log(ke)

    if args.check:
        emit(f"{path}")
        return stability_check(t, ke, emit, args.growing_tol, args.decaying_tol)

    if args.tail:
        res = find_tail(t, y, args.r2min, args.min_efolds, args.min_pts)
        if res is None:
            slope, r2 = linfit(t, y)
            emit(f"{path}")
            emit(f"   no growing tail found (r2min={args.r2min}, "
                 f"min_efolds={args.min_efolds}) -> looks NEUTRAL/decaying")
            emit(f"   full-series fit: sigma={0.5*slope:+.5f}  R2={r2:.4f}")
            return
        i, slope, r2 = res
        sigma = 0.5*slope
        efolds_KE = y[-1] - y[i]
        emit(f"{path}")
        emit(f"   TAIL fit  onset t={t[i]:.2f}  window t=[{t[i]:.2f},{t[-1]:.2f}]  "
             f"n={t.size-i}")
        emit(f"   sigma={sigma:+.5f}   (KE-rate={slope:+.5f})   R2={r2:.5f}   "
             f"KE e-folds in tail={efolds_KE:.1f}")
        return sigma, t[i]
    else:
        tmin = args.tmin if args.tmin is not None else t.min()
        tmax = args.tmax if args.tmax is not None else t.max()
        msk = (t >= tmin) & (t <= tmax)
        if msk.sum() < 3:
            emit(f"{path}: <3 points in [{tmin},{tmax}]"); return
        slope, r2 = linfit(t[msk], y[msk])
        emit(f"{path}  window t=[{tmin:g},{tmax:g}]  n={int(msk.sum())}  "
             f"sigma={0.5*slope:+.5f}  (KE-rate={slope:+.5f})  R2={r2:.5f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+")
    ap.add_argument("--tmin", type=float, default=None)
    ap.add_argument("--tmax", type=float, default=None)
    ap.add_argument("--tail", action="store_true",
                    help="auto-detect growth onset and fit only the exponential tail")
    ap.add_argument("--check", action="store_true",
                    help="fit the oscillation ENVELOPE (peaks) and report a "
                         "GROWING / DECAYING / NEUTRAL(STABLE) verdict")
    ap.add_argument("--growing-tol", type=float, default=3e-3,
                    help="envelope sigma above this counts as growing (default 3e-3)")
    ap.add_argument("--decaying-tol", type=float, default=3e-3,
                    help="envelope sigma below -this counts as decaying (default 3e-3)")
    ap.add_argument("--r2min", type=float, default=0.999,
                    help="R^2 threshold for accepting the tail fit (default 0.999)")
    ap.add_argument("--min-efolds", type=float, default=5.0,
                    help="tail must span >= this many e-folds of KE (default 5)")
    ap.add_argument("--min-pts", type=int, default=8,
                    help="minimum points in the tail (default 8)")
    ap.add_argument("--log", default="fit_growth.log",
                    help="path to save the report (default: fit_growth.log)")
    a = ap.parse_args()

    _lines = []
    def emit(s=""):
        print(s); _lines.append(s)

    import datetime
    mode = "check" if a.check else ("tail" if a.tail else "window")
    emit("="*70)
    emit(f"  fit_growth.py   {datetime.datetime.now():%Y-%m-%d %H:%M:%S}   "
         f"mode={mode}")
    emit("="*70)

    summary = []
    for p in a.logs:
        out = analyse(p, a, emit)
        emit("")
        if a.tail and isinstance(out, tuple):
            summary.append((p, out[0], out[1]))
        elif a.check and isinstance(out, tuple):
            summary.append((p, out[0], out[1]))

    if a.tail and len(summary) > 1:
        emit("-"*70)
        emit("  summary            sigma      onset_t")
        for p, sig, onset in summary:
            emit(f"    {p[-28:]:28s} {sig:+8.5f}   {onset:8.2f}")
    elif a.check and len(summary) > 1:
        emit("-"*70)
        emit("  summary            sigma_env    verdict")
        for p, sig, verdict in summary:
            emit(f"    {p[-28:]:28s} {sig:+8.5f}   {verdict}")
    emit("="*70)

    with open(a.log, "w") as f:
        f.write("\n".join(_lines) + "\n")
    print(f"[saved] {a.log}")


if __name__ == "__main__":
    main()