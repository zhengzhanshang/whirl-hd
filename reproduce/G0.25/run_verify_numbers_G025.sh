#!/usr/bin/env bash
# =====================================================================
# Re-check the Sect. 2.5 and Sect. 4.2-4.4 numbers with the CURRENT
# lsq_from_data.py (same version as the Table 1 and polar-mask reruns).
# All runs: G=0.25, m=5, reference pair (0.19,0.23), tmin=10.
# About 20 runs of a few seconds each.
#
# Run from the G=0.25 sweep directory on IGAM12WS:
#   /home/zshang/dedalus3/HD_SW_alpha_spherical_paper/only_h_small_G=0.25_paper/
#     vary_s2s4_G=0.25_nu4-h=1e-4_paper/sweep_s2s4_pairs
#   bash run_verify_numbers_G025.sh
#
# Paper value being checked is printed with each run as "EXPECT".
# =====================================================================
LSQ=${LSQ:-lsq_from_data.py}
LOG=${LOG:-verify_numbers_G025.log}
C=s2=0.19_s4=0.23
[ -f "$LSQ" ] || { echo "!!! $LSQ not found in $(pwd)"; exit 1; }
[ -f "$C/input.json" ] || { echo "!!! $C/input.json missing"; exit 1; }

exec > >(tee "$LOG") 2>&1
echo "#####################################################################"
echo "# VERIFY paper numbers  G=0.25 m=5  case $C"
echo "# started : $(date -Is)   host: $(hostname)"
echo "# lsq     : $LSQ   md5 $(md5sum "$LSQ" | cut -d' ' -f1)"
echo "#####################################################################"
echo "# --- defaults of the flags these runs rely on (from --help) ---"
python3 "$LSQ" --help 2>&1 | grep -i -A1 -- "--noise-seed\|--dropout\|--gap\|visc" | head -30
echo

run () {   # run LABEL EXPECT args...
  local L="$1" E="$2"; shift 2
  echo
  echo "====================================================================="
  echo "=== RUN=$L"
  echo "=== EXPECT: $E"
  echo "====================================================================="
  set -x
  python3 "$LSQ" --input "$C/input.json" --tmin 10 "$@"
  { rc=$?; set +x; } 2>/dev/null
  echo "--- exit=$rc  RUN=$L"
}

DEG=(--limb-weight best --mask-poles 60 --cheb 20 --svd 2)

# Sect. 2.5: K=20 on clean data (paper: residual 0.48 vs 6.2e-3 at K=30)
run clean_K20_r0   "residual 0.48 (K=30 gives 2.0e-3 in Table 1)" --cheb 20 --svd 0

# Sect. 4.2: centre-to-limb weighting + 60 deg mask, fixed K=20 r=2
run limb_n0.03     "s4 0.21%"  --noise 0.03 --noise-seed 0 "${DEG[@]}"
run limb_n0.1      "s4 0.72%"  --noise 0.1  --noise-seed 0 "${DEG[@]}"
run limb_n0.3      "s4 1.60%"  --noise 0.3  --noise-seed 0 "${DEG[@]}"
run limb_n1.0      "s4 1.70%"  --noise 1.0  --noise-seed 0 "${DEG[@]}"
run limb_n1.0_nomask   "s4 99.7%" --noise 1.0 --noise-seed 0 --limb-weight best --cheb 20 --svd 2
run limb_n0.3_nomask_r0 "s4 25.8%" --noise 0.3 --noise-seed 0 --limb-weight best --cheb 20 --svd 0
run limb_n0.3_mask_r0   "s4 16.1%" --noise 0.3 --noise-seed 0 --limb-weight best --mask-poles 60 --cheb 20 --svd 0

# Sect. 4.3: irregular sampling, noise 0.1 + limb + 60 deg mask, K=20 r=2
for S in 0 1 2 3 4; do
  run base_seed$S  "seed 0: 0.32/0.49%; 5-seed scatter +-0.3 / +-1.1%" --noise 0.1 --noise-seed $S "${DEG[@]}"
done
run dropout0.2     "0.00/0.25%"  --noise 0.1 --noise-seed 0 "${DEG[@]}" --dropout 0.2
run dropout0.4     "0.00/0.11%"  --noise 0.1 --noise-seed 0 "${DEG[@]}" --dropout 0.4
run gap80_120      "0.31/0.46%, dt spread 80x" --noise 0.1 --noise-seed 0 "${DEG[@]}" --gap 80 120

# Sect. 4.4: cadence, reference pair
run tstride2       "0.00/0.44%"  --noise 0.1 --noise-seed 0 "${DEG[@]}" --tstride 2
run tstride4       "5.79/5.99%"  --noise 0.1 --noise-seed 0 "${DEG[@]}" --tstride 4

echo
echo "#####################################################################"
echo "# SUMMARY  (s2 err, s4 err, cond, residual)  vs  EXPECT"
echo "#####################################################################"
awk '/^=== RUN=/{r=$2} /^=== EXPECT:/{sub(/^=== EXPECT: /,""); e=$0}
     /loaded N_t=/{nt=$2}
     /^ +s2 = /{s2=$NF} /^ +s4 = /{s4=$NF}
     /cond\(\[c2,c4\]\)/{printf "%-22s s2 %-8s s4 %-8s cond %-6s resid %-10s N_t %s | EXPECT %s\n", r, s2, s4, $3, $6, nt, e}' "$LOG"
echo "# finished : $(date -Is)"
