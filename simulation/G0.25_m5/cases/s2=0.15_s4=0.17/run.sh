#!/bin/bash
#SBATCH --job-name=🐸_rossby_G025
#SBATCH -p ib1
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=128G
#SBATCH -o job-%J.out
#SBATCH -e job-%J.err

# ─────────────────────────────────────────────────────────────────────────────
# G=0.25 sectoral Rossby run, seeded from an EVP eigenvector.
#
# Two stages:
#   1. evp_spherical_hd.py  — SERIAL, ~30 s. Solves the eigenvalue problem and
#      writes evp_ic.h5, the eigenvector of the EXACT operator stage 2 will
#      integrate. Both stages read the same input.json, so the parameters
#      cannot drift apart.
#   2. spherical-2d-hd-s2s4.py — MPI, ~3.5 h on 48 ranks for t=200.
#
# Files required in this directory:
#   run.sh                      (this)
#   input.json                  (parameters)
#   sw_model.py                 (shared: bases/fields/equations/BCs)
#   spherical-2d-hd-s2s4.py     (IVP; imports sw_model)
#   evp_spherical_hd.py         (EVP; imports sw_model)
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

# NOTE: `set -u` is deliberately NOT on for this block. `conda activate` first
# DEACTIVATES the Intel OneAPI intelpython env the cluster loads by default,
# which sources mpivars.deactivate.sh -> reads $SETVARS_CALL unbound -> under
# `set -u` that is fatal and kills the job before anything runs. Conda's own
# activation hooks are not -u clean either. Turn -u on afterwards, where it is
# actually useful (catching our own typos).
source ~/miniforge3/etc/profile.d/conda.sh
conda activate dedalus3
set -u

export LD_LIBRARY_PATH=~/miniforge3/envs/dedalus3/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
ulimit -l unlimited

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

INPUT=input.json
IC=evp_ic.h5
NRANKS="${SLURM_NTASKS:-48}"

# ── Preflight ────────────────────────────────────────────────────────────────
# Fail in 1 second rather than 20 minutes into a 48-core allocation.
echo "=============================================================="
echo " run.sh   $(date '+%Y-%m-%d %H:%M:%S')   dir=$(pwd)"
echo "=============================================================="
for f in "$INPUT" sw_model.py spherical-2d-hd-s2s4.py evp_spherical_hd.py; do
    if [[ ! -f "$f" ]]; then
        echo "MISSING: $f  — see the file list in the header of this script."
        exit 1
    fi
done
python3 -c "import dedalus.public, h5py, scipy, numpy" || {
    echo "conda env 'dedalus3' is missing a dependency"; exit 1; }

echo "--- input.json ---"
cat "$INPUT"
echo

# ── Stage 1: eigenvector IC (serial) ─────────────────────────────────────────
# Deleted first so a stale IC from a previous parameter set can never be reused.
# The IVP re-checks m/Nmu/G/s2/s4/nu4/nu4_h against the IC's attrs and refuses
# to start on a mismatch, but not creating the mismatch is better.
echo "=== Stage 1/2: EVP -> $IC  (serial) ==="
rm -f "$IC"
python3 evp_spherical_hd.py --input "$INPUT" --write-ic "$IC" --pick rossby
if [[ ! -f "$IC" ]]; then
    echo "EVP did not write $IC — aborting before burning the allocation."
    exit 1
fi
echo

# ── Stage 2: IVP (MPI) ───────────────────────────────────────────────────────
echo "=== Stage 2/2: IVP on $NRANKS ranks ==="
echo "  WATCH THE FIRST LOG LINES:"
echo "    - 'IC: EVP eigenvector from $IC'  (not 'balanced streamfunction')"
echo "    - <h2> must NOT jump 0 -> 4e-7. The eigenvector carries its own"
echo "      balanced h; a spike means geostrophic adjustment, i.e. the IC"
echo "      did not take. Kill the job in the first minute if you see it."
echo
mpiexec -n "$NRANKS" python3 spherical-2d-hd-s2s4.py "$INPUT"

echo
echo "=== done: $(date '+%Y-%m-%d %H:%M:%S') ==="
echo "Next, in order — do NOT skip the gate:"
echo "  1. python3 ../postrun_analyze.py --input $INPUT"
echo "       gate: CLEAN SINGLE MODE: YES, omega ~ -0.4374 (RETROGRADE)"
echo "  2. python3 ../lsq_from_data.py --input $INPUT --tmin <from step 1> --scan"