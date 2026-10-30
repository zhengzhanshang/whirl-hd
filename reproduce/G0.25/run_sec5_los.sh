#!/usr/bin/env bash
# =====================================================================
# Section 5 -- end-to-end LOS chain, two additional (s2,s4) pairs.
# For each noise level: invert synthetic Dopplergrams, then recover from
# BOTH fields -- --field true (control, must reproduce the direct m-mode
# answer) and --field pred (the actual end-to-end number).
# Config matches Table 2: B0=7.25, 60 deg mask, cheb=20, svd=2 (defaults).
# eq-subset 2 is the headline; eq-subset 12 runs alongside as the
# cross-equation redundancy check.
#
# CHECK IN THE LOG: fitted lambda must be Re ~ sigma, Im ~ |omega| for
# THIS case (from its postrun_analyse log). Near zero => the fit failed
# and every number below it is meaningless.
#
# All output goes to ONE log file.
# =====================================================================
LOS=${LOS:-los_inversion.py}
LSQ=${LSQ:-lsq_from_data.py}
LOG=${LOG:-sec5_los_newpairs.log}
B0=${B0:-7.25}
CASES=("s2=0.29_s4=0.03" "s2=0.03_s4=0.33")
# "noise seed" pairs; noise=0 is the noiseless control run
RUNS=("0:0" "0.01:0" "0.05:0" "0.05:1" "0.05:2" "0.05:3" "0.05:4" "0.10:0")
SUBSETS=(2 12)

exec > >(tee "$LOG") 2>&1

echo "#####################################################################"
echo "# SECTION 5  LOS chain, B0=$B0"
echo "# started : $(date -Is)"
echo "# host    : $(hostname)"
echo "# cwd     : $(pwd)"
echo "# los     : $LOS"
echo "# lsq     : $LSQ"
echo "#####################################################################"

for C in "${CASES[@]}"; do
  for R in "${RUNS[@]}"; do
    N="${R%%:*}"; S="${R##*:}"
    if [ "$N" = "0" ]; then
      OUT="los_${C}_clean"; NOISEFLAG=()
    else
      OUT="los_${C}_n${N}_s${S}"; NOISEFLAG=(--noise "$N" --noise-seed "$S")
    fi

    echo
    echo "#####################################################################"
    echo "### INVERSION  CASE=$C  NOISE=$N  SEED=$S  OUT=$OUT"
    echo "### $(date -Is)"
    echo "#####################################################################"
    set -x
    python3 "$LOS" --input "$C/input.json" --tmin 10 --b0 "$B0" \
      "${NOISEFLAG[@]}" --out "$OUT"
    { rc=$?; set +x; } 2>/dev/null
    echo "--- inversion exit=$rc  OUT=$OUT"
    if [ "$rc" -ne 0 ]; then
      echo "!!! inversion FAILED, skipping recoveries for $OUT"
      continue
    fi

    for F in true pred; do
      for E in "${SUBSETS[@]}"; do
        echo
        echo "====================================================================="
        echo "=== RECOVERY  CASE=$C  NOISE=$N  SEED=$S  FIELD=$F  EQ=$E"
        echo "====================================================================="
        set -x
        python3 "$LSQ" --ckpt "$OUT" --field "$F" \
          --input "$C/input.json" --tmin 10 --mask-poles 60 --eq-subset "$E"
        { rc=$?; set +x; } 2>/dev/null
        echo "--- exit=$rc  CASE=$C  NOISE=$N  SEED=$S  FIELD=$F  EQ=$E"
      done
    done
  done
done

echo
echo "#####################################################################"
echo "# SECTION 5 finished : $(date -Is)"
echo "#####################################################################"
