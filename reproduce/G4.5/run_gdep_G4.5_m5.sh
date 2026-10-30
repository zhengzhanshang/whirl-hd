#!/usr/bin/env bash
# =====================================================================
# G-dependence section -- re-analysis of the G=4.5, m=5 sweep with the
# SAME estimator and settings as the G=0.25 paper. No new simulations.
#
# Run from:  paper-redo_G=4.5/m=3-15/m=5_s2s4_sweep/
#
# Four parts, each writes its own log (STAGES selects which to run):
#   A  sweep     : all seven pairs, clean data, --scan      -> Table 1 at G=4.5
#   B  mask      : reference pair, cutoffs none/60/40/20    -> Fig. 1a at G=4.5
#   C  los       : reference pair, LOS chain, noise 0.001-0.1, three
#                  equation subsets, field true + pred       -> Fig. 2 / Table 2
#                  at G=4.5, and the h = w/lambda check
#
#   STAGES="A B C" bash run_gdep_G4.5_m5.sh      (default: all)
#   STAGES="A"     bash run_gdep_G4.5_m5.sh      (only the sweep)
#
# Cases without a snapshots/ directory are skipped with a warning, so this
# can run now on the six finished pairs and again for s2=0.15_s4=0.15
# once that job has finished:
#   STAGES="A" CASES_A="s2=0.15_s4=0.15" bash run_gdep_G4.5_m5.sh
#
# Settings copied from the G=0.25 scripts (run_sec41_mask.sh,
# run_sec5_los.sh): --tmin 10, B0=7.25, 60 deg mask, cheb=20 / svd=2
# defaults for degraded data, 5 noise seeds at 0.05.
# =====================================================================
LSQ=${LSQ:-lsq_from_data.py}
LOS=${LOS:-los_inversion.py}
B0=${B0:-7.25}
STAGES=${STAGES:-"A B C"}
REF=${REF:-s2=0.19_s4=0.23}

CASES_A=${CASES_A:-"s2=0.29_s4=0.03 s2=0.25_s4=0.19 s2=0.19_s4=0.23 s2=0.15_s4=0.17 s2=0.15_s4=0.15 s2=0.16_s4=0.12 s2=0.03_s4=0.33"}
MASKS=("none" 60 40 20)
# "noise:seed"; noise 0 is the noiseless control inversion
LOS_RUNS=("0:0" "0.001:0" "0.01:0" "0.02:0" "0.05:0" "0.05:1" "0.05:2" "0.05:3" "0.05:4" "0.10:0")
# "full" = no --eq-subset flag (all three equations); see the preflight
# block in the log to confirm that is what the default does.
SUBSETS=("full" 12 2)

for f in "$LSQ" "$LOS"; do
  [ -f "$f" ] || { echo "!!! $f not found in $(pwd) -- set LSQ=/LOS= to its path"; exit 1; }
done

have_snaps() { [ -d "$1/snapshots" ] && [ -n "$(ls -A "$1/snapshots" 2>/dev/null)" ]; }

header() {   # $1 = title
  echo "#####################################################################"
  echo "# $1"
  echo "# started : $(date -Is)"
  echo "# host    : $(hostname)"
  echo "# cwd     : $(pwd)"
  echo "# lsq     : $LSQ    los : $LOS"
  echo "#####################################################################"
  echo "### PREFLIGHT: equation-subset option as the script defines it"
  python3 "$LSQ" --help 2>&1 | grep -A3 -- '--eq-subset' || echo "    (no --eq-subset in --help)"
}

show_input() {  # print the case parameters so the log proves G=4.5, m=5
  python3 -c "import json,sys; d=json.load(open(sys.argv[1])); \
print('    input:', {k: d[k] for k in d if k in ('G','m','s2','s4','nu4','nu4_h','Nmu','Nphi','snap_dt','t_end','stop_sim_time')})" \
    "$1/input.json"
}

# ---------------------------------------------------------------- A
if [[ " $STAGES " == *" A "* ]]; then
(
  exec > >(tee "gdep_A_sweep.log") 2>&1
  header "G-DEP A  six/seven-pair sweep, clean data, G=4.5 m=5"
  for C in $CASES_A; do
    echo; echo "====================================================================="
    echo "=== CASE=$C  (clean, --scan)"; echo "=== $(date -Is)"
    echo "====================================================================="
    if ! have_snaps "$C"; then echo "!!! $C has no snapshots yet -- SKIPPED"; continue; fi
    show_input "$C"
    set -x
    python3 "$LSQ" --input "$C/input.json" --tmin 10 --scan
    { rc=$?; set +x; } 2>/dev/null
    echo "--- exit=$rc  CASE=$C"
  done
  echo; echo "# G-DEP A finished : $(date -Is)"
)
fi

# ---------------------------------------------------------------- B
if [[ " $STAGES " == *" B "* ]]; then
(
  exec > >(tee "gdep_B_mask.log") 2>&1
  header "G-DEP B  polar masking, clean data, $REF, G=4.5 m=5"
  if ! have_snaps "$REF"; then echo "!!! $REF has no snapshots -- stage B SKIPPED"; exit 0; fi
  show_input "$REF"
  for M in "${MASKS[@]}"; do
    if [ "$M" = "none" ]; then MASKFLAG=(); else MASKFLAG=(--mask-poles "$M"); fi
    echo; echo "====================================================================="
    echo "=== CASE=$REF  MASK=$M"; echo "=== $(date -Is)"
    echo "====================================================================="
    set -x
    python3 "$LSQ" --input "$REF/input.json" --tmin 10 --scan "${MASKFLAG[@]}"
    { rc=$?; set +x; } 2>/dev/null
    echo "--- exit=$rc  CASE=$REF  MASK=$M"
  done
  echo; echo "# G-DEP B finished : $(date -Is)"
)
fi

# ---------------------------------------------------------------- C
# CHECK IN THE LOG: the fitted lambda must match this case's EVP mode,
# Re ~ sigma (now very small: the mode barely decays at G=4.5) and
# Im ~ |omega|. A near-zero Im means the fit failed and every number
# below it is meaningless. The amplitude-threshold fit (--amp-frac) was
# tuned for the strongly damped G=0.25 mode -- note what it does here.
if [[ " $STAGES " == *" C "* ]]; then
(
  exec > >(tee "gdep_C_los.log") 2>&1
  header "G-DEP C  LOS chain, $REF, B0=$B0, G=4.5 m=5"
  if ! have_snaps "$REF"; then echo "!!! $REF has no snapshots -- stage C SKIPPED"; exit 0; fi
  show_input "$REF"
  for R in "${LOS_RUNS[@]}"; do
    N="${R%%:*}"; S="${R##*:}"
    if [ "$N" = "0" ]; then OUT="gdep_los_${REF}_clean"; NOISEFLAG=()
    else OUT="gdep_los_${REF}_n${N}_s${S}"; NOISEFLAG=(--noise "$N" --noise-seed "$S"); fi

    echo; echo "#####################################################################"
    echo "### INVERSION  CASE=$REF  NOISE=$N  SEED=$S  OUT=$OUT"
    echo "### $(date -Is)"
    echo "#####################################################################"
    set -x
    python3 "$LOS" --input "$REF/input.json" --tmin 10 --b0 "$B0" "${NOISEFLAG[@]}" --out "$OUT"
    { rc=$?; set +x; } 2>/dev/null
    echo "--- inversion exit=$rc  OUT=$OUT"
    if [ "$rc" -ne 0 ]; then echo "!!! inversion FAILED, skipping recoveries for $OUT"; continue; fi

    for F in true pred; do
      for E in "${SUBSETS[@]}"; do
        if [ "$E" = "full" ]; then EQFLAG=(); else EQFLAG=(--eq-subset "$E"); fi
        echo; echo "====================================================================="
        echo "=== RECOVERY  CASE=$REF  NOISE=$N  SEED=$S  FIELD=$F  EQ=$E"
        echo "====================================================================="
        set -x
        python3 "$LSQ" --ckpt "$OUT" --field "$F" --input "$REF/input.json" \
          --tmin 10 --mask-poles 60 "${EQFLAG[@]}"
        { rc=$?; set +x; } 2>/dev/null
        echo "--- exit=$rc  CASE=$REF  NOISE=$N  SEED=$S  FIELD=$F  EQ=$E"
      done
    done
  done
  echo; echo "# G-DEP C finished : $(date -Is)"
)
fi
