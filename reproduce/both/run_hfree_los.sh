#!/usr/bin/env bash
# =====================================================================
# h-free recovery on the LOS-inversion checkpoints.
#
# Same inverted fields as stage C of run_gdep_G4.5_m5.sh, but recovered
# with lsq_hfree.py, which uses u_theta and u_phi ONLY (h eliminated
# analytically; G cancels). Tests whether the h-free route survives the
# noise levels where every h-based subset failed at G=4.5.
#
# Run it twice, once per G, from the directory holding the checkpoints:
#
#   G=4.5  (GSC, m=5_s2s4_sweep/):
#     CKPT_GLOB='gdep_los_s2=0.19_s4=0.23_*' INPUT=s2=0.19_s4=0.23/input.json \
#       TAG=G4.5 bash run_hfree_los.sh
#
#   G=0.25 (the original Table 2 checkpoints, after fetching them):
#     CKPT_GLOB='los_s2=0.19_s4=0.23_*' INPUT=s2=0.19_s4=0.23/input.json \
#       TAG=G0.25 bash run_hfree_los.sh
#
# For every checkpoint and both fields (true = control, pred = inverted),
# four recoveries:
#   fixed  cheb=20 svd=2   with nu4 terms      (matches the paper's settings)
#   fixed  cheb=20 svd=2   --no-visc           (solar limit; difference = nu4 bias)
#   --scan                 with nu4 terms      (best-of-scan, VALIDATION ONLY)
#   --scan                 --no-visc
#
# POLAR MASK: the paper fits |lat| < 60 deg. lsq_hfree.py has no
# --mask-poles; the newer version has --band-deg (band-limited Chebyshev
# fit). This script uses --band-deg 60 if the copy you run supports it,
# and otherwise warns loudly -- without it, the badly conditioned polar
# rows of the LOS inversion (cond up to ~7e3 above 75 deg) enter the fit.
# =====================================================================
HFREE=${HFREE:-lsq_hfree.py}
CKPT_GLOB=${CKPT_GLOB:?set CKPT_GLOB, e.g. 'gdep_los_s2=0.19_s4=0.23_*'}
INPUT=${INPUT:?set INPUT to the case input.json}
TAG=${TAG:-run}
BAND=${BAND:-60}
LOG=${LOG:-hfree_los_${TAG}.log}

[ -f "$HFREE" ] || { echo "!!! $HFREE not found in $(pwd)"; exit 1; }
[ -f "$INPUT" ] || { echo "!!! $INPUT not found"; exit 1; }

exec > >(tee "$LOG") 2>&1

echo "#####################################################################"
echo "# h-free recovery on LOS checkpoints   TAG=$TAG"
echo "# started : $(date -Is)"
echo "# host    : $(hostname)"
echo "# cwd     : $(pwd)"
echo "# hfree   : $HFREE"
echo "# input   : $INPUT"
echo "#####################################################################"

echo "### SELFTEST (coefficients vs sympy)"
python3 "$HFREE" --selftest 2>&1 | tail -3

if python3 "$HFREE" --help 2>&1 | grep -q -- '--band-deg'; then
  MASKFLAG=(--band-deg "$BAND")
  echo "### polar band: --band-deg $BAND"
else
  MASKFLAG=()
  echo "!!! WARNING: this lsq_hfree.py has no --band-deg. Polar rows are NOT"
  echo "!!! masked, unlike the paper's recoveries. Use the newer copy (from"
  echo "!!! m=3_s2s4_sweep/ or ring_diagram/) for numbers you intend to quote."
fi

shopt -s nullglob
CKPTS=( $CKPT_GLOB )
[ ${#CKPTS[@]} -gt 0 ] || { echo "!!! no checkpoints match '$CKPT_GLOB'"; exit 1; }
echo "### ${#CKPTS[@]} checkpoints: ${CKPTS[*]}"

for C in "${CKPTS[@]}"; do
  [ -f "$C/pinn_fields.h5" ] || { echo "!!! $C has no pinn_fields.h5 -- skipped"; continue; }
  for F in true pred; do
    for V in visc novisc; do
      if [ "$V" = "novisc" ]; then VFLAG=(--no-visc); else VFLAG=(); fi
      for MODE in fixed scan; do
        if [ "$MODE" = "scan" ]; then SFLAG=(--scan); else SFLAG=(--cheb 20 --svd 2); fi
        echo
        echo "====================================================================="
        echo "=== CKPT=$C  FIELD=$F  VISC=$V  MODE=$MODE"
        echo "====================================================================="
        set -x
        python3 "$HFREE" --ckpt "$C" --field "$F" --input "$INPUT" --tmin 10 \
          "${SFLAG[@]}" "${VFLAG[@]}" "${MASKFLAG[@]}"
        { rc=$?; set +x; } 2>/dev/null
        echo "--- exit=$rc  CKPT=$C  FIELD=$F  VISC=$V  MODE=$MODE"
      done
    done
  done
done

echo
echo "#####################################################################"
echo "# finished : $(date -Is)"
echo "#####################################################################"
