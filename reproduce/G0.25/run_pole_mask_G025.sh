#!/usr/bin/env bash
# =====================================================================
# Sect. 4.1 / Fig. 1a rerun -- polar masking, G=0.25, m=5, clean data,
# full three-equation stack, FIXED K=30, r=0, t >= 10 (same as Table 1).
#
# Why: the Sect. 4.1 numbers and the POLE block in make_figures.py come
# from an older lsq_from_data.py (reference pair no-mask 0.10/0.12%, while
# the new Table 1 gives 0.02/0.02%). The stress pairs were also best-of-
# K-scan. This reruns all three pairs at four cutoffs with ONE script
# version and ONE fixed setting.
#
# Run from the G=0.25 sweep directory on IGAM12WS (same as Table 1):
#   /home/zshang/dedalus3/HD_SW_alpha_spherical_paper/only_h_small_G=0.25_paper/
#     vary_s2s4_G=0.25_nu4-h=1e-4_paper/sweep_s2s4_pairs
#
#   bash run_pole_mask_G025.sh
#
# The no-mask row must reproduce Table 1 (sanity check):
#   (0.19,0.23) 0.02/0.02 cond 7.18 | (0.29,0.03) 0.03/0.37 cond 7.20
#   (0.03,0.33) 0.07/0.00 cond 6.61
# CHECK IN THE LOG: N_t=379-380 for every case.
# =====================================================================
LSQ=${LSQ:-lsq_from_data.py}
LOG=${LOG:-pole_mask_rerun_G025.log}
CASES=("s2=0.19_s4=0.23" "s2=0.29_s4=0.03" "s2=0.03_s4=0.33")
CUTS=(none 60 40 20)

[ -f "$LSQ" ] || { echo "!!! $LSQ not found in $(pwd)"; exit 1; }

exec > >(tee "$LOG") 2>&1

echo "#####################################################################"
echo "# POLAR MASK rerun  G=0.25 m=5  clean, full stack, K=30 r=0, tmin=10"
echo "# started : $(date -Is)"
echo "# host    : $(hostname)"
echo "# cwd     : $(pwd)"
echo "# lsq     : $LSQ   md5 $(md5sum "$LSQ" | cut -d' ' -f1)"
echo "#####################################################################"

for C in "${CASES[@]}"; do
  if [ ! -f "$C/input.json" ]; then echo "!!! $C/input.json missing -- SKIPPED"; continue; fi
  for M in "${CUTS[@]}"; do
    echo
    echo "====================================================================="
    echo "=== CASE=$C MASK=$M"
    echo "=== $(date -Is)"
    echo "====================================================================="
    MASK=()
    [ "$M" != none ] && MASK=(--mask-poles "$M")
    set -x
    python3 "$LSQ" --input "$C/input.json" --tmin 10 --cheb 30 --svd 0 "${MASK[@]}"
    { rc=$?; set +x; } 2>/dev/null
    echo "--- exit=$rc  CASE=$C MASK=$M"
  done
done

echo
echo "#####################################################################"
echo "# SUMMARY (s2 err, s4 err, cond) -- for Sect. 4.1 and POLE in make_figures.py"
echo "#####################################################################"
awk '/^=== CASE=/{c=$2" "$3} /loaded N_t=/{nt=$0}
     /^ +s2 = /{s2=$NF} /^ +s4 = /{s4=$NF}
     /cond\(\[c2,c4\]\)/{print c, " s2", s2, " s4", s4, " ", $0; print "      ", nt}' "$LOG"
echo "# finished : $(date -Is)"
