#!/usr/bin/env bash
# =====================================================================
# Re-run the two recovered drivers with the ARCHIVED (pinned) code:
#   run_sec44_cadence.sh  -> Sect. 4.4 / Fig. 1b, stress pairs
#   run_sec5_los.sh       -> Table 2, stress pairs
# and print the numbers next to the values in the paper.
#
# Both drivers ran on 2026-09-21 with an older lsq_from_data.py. This
# confirms the paper's stress-pair numbers with the pinned version
# (md5 176d096a...), the one archived in the repository.
#
# Run on IGAM12WS FROM THE G=0.25 SWEEP FOLDER (the one holding the
# s2=..._s4=.../ case folders; takes a few minutes, no SLURM needed):
#   ARCH=<path to scripts_for_papers_data_avail> bash rerun_stress_pairs_pinned.sh
# run_sec44_cadence.sh and run_sec5_los.sh must sit next to this script,
# or set DRIVERS=<folder holding them>.
#
# SAFE: works in its own folder (WORK, default ./rerun_pinned_<date>) with
# symlinks to the case folders, so the September logs and the LOS
# checkpoints next to the cases are never overwritten.
# =====================================================================
set -uo pipefail
SWEEP=${SWEEP:-$PWD}
ARCH=${ARCH:?set ARCH to the scripts_for_papers_data_avail folder}
LSQ=${LSQ:-$ARCH/analysis/lsq_from_data.py}
LOS=${LOS:-$ARCH/analysis/los_inversion.py}
HERE=${DRIVERS:-"$(cd "$(dirname "$0")" && pwd)"}
WORK=${WORK:-$PWD/rerun_pinned_$(date +%Y%m%d)}
PIN_LSQ=176d096ac17a3371dc06f164e662de14
PIN_LOS=320d5a4771713a42c636033b0932c0e9

# ── preflight ────────────────────────────────────────────────────────────────
for f in "$LSQ" "$LOS" "$HERE/run_sec44_cadence.sh" "$HERE/run_sec5_los.sh"; do
  [ -f "$f" ] || { echo "!!! not found: $f"; exit 1; }
done
[ "$(md5sum "$LSQ" | cut -d' ' -f1)" = "$PIN_LSQ" ] || { echo "!!! $LSQ is NOT the pinned version"; exit 1; }
[ "$(md5sum "$LOS" | cut -d' ' -f1)" = "$PIN_LOS" ] || { echo "!!! $LOS is NOT the archived version"; exit 1; }
for C in s2=0.29_s4=0.03 s2=0.03_s4=0.33; do
  [ -f "$SWEEP/$C/input.json" ] || { echo "!!! $SWEEP/$C/input.json missing"; exit 1; }
done
[ -e "$WORK" ] && { echo "!!! $WORK exists -- remove it or set WORK=..."; exit 1; }
mkdir -p "$WORK" && cd "$WORK"
ln -s "$SWEEP/s2=0.29_s4=0.03" . && ln -s "$SWEEP/s2=0.03_s4=0.33" .
echo "### work dir : $WORK"
echo "### lsq      : $LSQ  (pinned OK)"
echo "### los      : $LOS  (archived OK)"

# ── the two drivers, unchanged, with the pinned code ────────────────────────
LSQ="$LSQ" LOG=sec44_cadence_pinned.log bash "$HERE/run_sec44_cadence.sh" > /dev/null; sleep 2
LSQ="$LSQ" LOS="$LOS" LOG=sec5_los_pinned.log bash "$HERE/run_sec5_los.sh" > /dev/null; sleep 2
echo "### failed runs (should be none):"
grep -hE '^--- .*exit=[^0]' sec44_cadence_pinned.log sec5_los_pinned.log || echo "    none"

# ── summary ──────────────────────────────────────────────────────────────────
python3 - <<'PY'
import re, statistics as st
def runs(log, head):
    """yield (banner fields, s2, s4) for each recovery in a driver log."""
    cur, s2 = None, None
    for ln in open(log):
        if ln.startswith(head):
            cur = dict(kv.split("=", 1) for kv in ln.split()[1:] if "=" in kv)
            # CASE=s2=0.29_s4=0.03 keeps its own '=' after the first split
        m = re.match(r"\s+s2 = ([-\d.eE+]+)\s+true ([-\d.eE+]+)", ln)
        if m: s2 = (float(m[1]), float(m[2]))
        m = re.match(r"\s+s4 = ([-\d.eE+]+)\s+true ([-\d.eE+]+)", ln)
        if m and cur is not None and s2 is not None:
            s4 = (float(m[1]), float(m[2]))
            yield cur, s2, s4
            s2 = None
rel = lambda v: abs(v[0]-v[1])/abs(v[1])*100

# Sect. 4.4: 5-seed mean +- std per (case, tstride); paper values = Fig. 1b (make_figures.py CAD)
PAPER = {("s2=0.29_s4=0.03","1"):(0.19,7.07), ("s2=0.29_s4=0.03","2"):(0.49,10.10),
         ("s2=0.29_s4=0.03","4"):(5.59,76.29), ("s2=0.03_s4=0.33","1"):(1.51,0.43),
         ("s2=0.03_s4=0.33","2"):(1.18,0.35), ("s2=0.03_s4=0.33","4"):(18.93,1.97)}
cells = {}
for b, s2, s4 in runs("sec44_cadence_pinned.log", "=== CASE="):
    cells.setdefault((b["CASE"], b["TSTRIDE"]), []).append((rel(s2), rel(s4)))
print("\n=== Sect. 4.4 cadence (5 seeds): pinned mean +- std  |  paper (Fig. 1b)")
for k in sorted(cells):
    e2 = [x[0] for x in cells[k]]; e4 = [x[1] for x in cells[k]]
    sd = lambda x: st.stdev(x) if len(x) > 1 else 0.0
    p = PAPER.get(k, (float("nan"),)*2)
    print(f"  {k[0]} tstride {k[1]}  n={len(e2)}  s2 {st.mean(e2):6.2f} +- {sd(e2):5.2f}%  "
          f"s4 {st.mean(e4):6.2f} +- {sd(e4):5.2f}%   |  paper {p[0]:.2f} / {p[1]:.2f}%")

# Table 2: one line per inversion and recovery; compare with the stress-pair rows of Table 2
print("\n=== Table 2, stress pairs (pinned): relative % and absolute error")
print(f"  {'case':16s} {'noise':>5s} {'seed':>4s} {'field':5s} {'eq':>3s}   "
      f"{'s2 %':>7s} {'s4 %':>7s}   {'|ds2|':>9s} {'|ds4|':>9s}")
for b, s2, s4 in runs("sec5_los_pinned.log", "=== RECOVERY"):
    print(f"  {b['CASE']:16s} {b['NOISE']:>5s} {b['SEED']:>4s} {b['FIELD']:5s} {b['EQ']:>3s}   "
          f"{rel(s2):7.2f} {rel(s4):7.2f}   {abs(s2[0]-s2[1]):9.2e} {abs(s4[0]-s4[1]):9.2e}")
print("\n  noise 0.05: average the five seeds before comparing with Table 2.")
print("  field=true rows are the controls (Table 2 control rows); field=pred the end-to-end values.")
PY
echo
echo "### logs: $WORK/sec44_cadence_pinned.log  $WORK/sec5_los_pinned.log"
