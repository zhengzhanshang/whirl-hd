#!/usr/bin/env bash
# =====================================================================
# Table 1 rerun -- G=0.25, m=5, clean data, full three-equation stack,
# FIXED K=30, r=0 (the settings the Table 1 caption states), t >= 10.
#
# Why all six pairs: Table 1 was assembled in July from an older
# lsq_from_data.py, and row (0.15,0.17) was re-run later with the
# refactored script. Rows (0.29,0.03) and (0.03,0.33) disagree with the
# September log. Re-running all six with ONE script version gives a table
# that is consistent and reproducible from the archived code.
#
# Run from the G=0.25 sweep directory on IGAM12WS:
#   /home/zshang/dedalus3/HD_SW_alpha_spherical_paper/only_h_small_G=0.25_paper/
#     vary_s2s4_G=0.25_nu4-h=1e-4_paper/sweep_s2s4_pairs
#
#   bash run_table1_G025.sh
#
# --cheb 30 --svd 0 are passed EXPLICITLY: the script defaults are now the
# degraded-data values (cheb=20, svd=2), which are wrong for clean data.
#
# CHECK IN THE LOG: every case must load N_t=379-380 over t in [10.5, 200]
# (snap_dt = 0.5). The (0.15,0.17) run was written every 0.1 and is
# subsampled with --tstride 5 from t=10.5 to land on the same grid.
# =====================================================================
LSQ=${LSQ:-lsq_from_data.py}
LOG=${LOG:-table1_rerun_G025.log}
CASES=("s2=0.29_s4=0.03" "s2=0.25_s4=0.19" "s2=0.19_s4=0.23"
       "s2=0.15_s4=0.17" "s2=0.15_s4=0.15" "s2=0.03_s4=0.33")

[ -f "$LSQ" ] || { echo "!!! $LSQ not found in $(pwd)"; exit 1; }

exec > >(tee "$LOG") 2>&1

echo "#####################################################################"
echo "# TABLE 1 rerun  G=0.25 m=5  clean, full stack, K=30 r=0, tmin=10"
echo "# started : $(date -Is)"
echo "# host    : $(hostname)"
echo "# cwd     : $(pwd)"
echo "# lsq     : $LSQ   md5 $(md5sum "$LSQ" | cut -d' ' -f1)"
echo "#####################################################################"

for C in "${CASES[@]}"; do
  echo
  echo "====================================================================="
  echo "=== CASE=$C"
  echo "=== $(date -Is)"
  echo "====================================================================="
  if [ ! -f "$C/input.json" ]; then echo "!!! $C/input.json missing -- SKIPPED"; continue; fi
  # (0.15,0.17) was written every 0.1; keep every 5th snapshot from t=10.5 so
  # it is analysed on the same 10.5, 11.0, ... grid as the other five pairs
  EXTRA=(--tmin 10)
  [ "$C" = "s2=0.15_s4=0.17" ] && EXTRA=(--tmin 10.4 --tstride 5)
  set -x
  python3 "$LSQ" --input "$C/input.json" "${EXTRA[@]}" --cheb 30 --svd 0
  { rc=$?; set +x; } 2>/dev/null
  echo "--- exit=$rc  CASE=$C"
done

echo
echo "#####################################################################"
echo "# SUMMARY (s2 err, s4 err, cond, residual) -- paste into Table 1"
echo "#####################################################################"
awk '/^=== CASE=/{c=$2} /loaded N_t=/{nt=$0}
     /^ +s2 = /{s2=$NF} /^ +s4 = /{s4=$NF}
     /cond\(\[c2,c4\]\)/{print c, " s2", s2, " s4", s4, " ", $0; print "      ", nt}' "$LOG"
echo "# finished : $(date -Is)"