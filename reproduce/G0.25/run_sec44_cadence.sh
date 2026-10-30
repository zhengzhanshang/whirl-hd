#!/usr/bin/env bash
# =====================================================================
# Section 4.4 -- cadence sweep, two additional (s2,s4) pairs.
# DEGRADED config: noise 0.1 + centre-to-limb weighting + 60 deg mask,
# matching Fig. 1b. cheb=20 / svd=2 are the correct defaults here.
# Five noise seeds per cell -- single-seed numbers are not quotable.
# NOTE: tstride is NOT samples-per-period. Convert using each case's own
#       omega from its postrun_analyse log before plotting.
# All output goes to ONE log file.
# =====================================================================
LSQ=${LSQ:-lsq_from_data.py}
LOG=${LOG:-sec44_cadence_newpairs.log}
CASES=("s2=0.29_s4=0.03" "s2=0.03_s4=0.33")
SEEDS=(0 1 2 3 4)
STRIDES=(1 2 4)

exec > >(tee "$LOG") 2>&1

echo "#####################################################################"
echo "# SECTION 4.4  cadence, noise 0.1 + limb-weight best + mask 60"
echo "# started : $(date -Is)"
echo "# host    : $(hostname)"
echo "# cwd     : $(pwd)"
echo "# lsq     : $LSQ"
echo "#####################################################################"

for C in "${CASES[@]}"; do
  for T in "${STRIDES[@]}"; do
    for S in "${SEEDS[@]}"; do
      echo
      echo "====================================================================="
      echo "=== CASE=$C  TSTRIDE=$T  NOISE_SEED=$S"
      echo "=== $(date -Is)"
      echo "====================================================================="
      set -x
      python3 "$LSQ" --input "$C/input.json" --tmin 10 \
        --noise 0.1 --noise-seed "$S" \
        --limb-weight best --mask-poles 60 --tstride "$T"
      { rc=$?; set +x; } 2>/dev/null
      echo "--- exit=$rc  CASE=$C  TSTRIDE=$T  SEED=$S"
    done
  done
done

echo
echo "#####################################################################"
echo "# SECTION 4.4 finished : $(date -Is)"
echo "#####################################################################"
