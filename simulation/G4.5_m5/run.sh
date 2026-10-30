#!/bin/bash
#SBATCH --job-name=rossby_G4.5_m5
#SBATCH -p ib1
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=128G
#SBATCH -o job-%J.out
#SBATCH -e job-%J.err

# ─────────────────────────────────────────────────────────────────────────────
# G=4.5, m=l=5 sectoral Rossby run, seeded from an EVP eigenvector.
#
# G=4.5 is the photospheric value (g=109 m/s^2, H0=167 km at ~70 km depth).
# NOTHING from the earlier G=0.25 / m=5 campaign transfers: the mode, its
# frequency, the field amplitude ratios and every recovered number change.
#
# Two stages:
#   1. evp_spherical_hd.py  — SERIAL, ~30 s. Solves the eigenvalue problem and
#      writes evp_ic.h5, the eigenvector of the EXACT operator stage 2 will
#      integrate. Both stages read the same input.json, so the parameters
#      cannot drift apart.
#   2. spherical-2d-hd-s2s4.py — MPI, ~3.5 h on 48 ranks for t=200.
#      Compute is unchanged by snap_dt (max_timestep sets the step count),
#      but snap_dt=0.1 writes 2000 snapshots ~= 1.6 GB per case, ~9.4 GB for
#      the six-case sweep. Check quota before launching them all.
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
# Record the seeded mode, so the post-run gate is checked against THIS run
# and not against a remembered number from a different G and m.
python3 -c "import h5py,math,sys; f=h5py.File(sys.argv[1]); wr=float(f.attrs['omega_real']); wi=float(f.attrs['omega_imag']); print(f'  seeded mode: omega = {wr:+.6f} {wi:+.6f}j   T = {2*math.pi/abs(wr):.4f}'); print(f'  decay over t=200: exp(Im*200) = {math.exp(wi*200):.3e}'); print('  -> postrun_analyse must reproduce THIS omega')" "$IC"
echo

# ── Stage 2: IVP (MPI) ───────────────────────────────────────────────────────
echo "=== Stage 2/2: IVP on $NRANKS ranks ==="
echo "  WATCH THE FIRST LOG LINES:"
echo "    - 'IC: EVP eigenvector from $IC'  (not 'balanced streamfunction')"
echo "    - <h2> must NOT jump from ~0 to a finite plateau in the first few"
echo "      writes. The eigenvector carries its own balanced h; a spike means"
echo "      geostrophic adjustment, i.e. the IC did not take. Kill the job in"
echo "      the first minute if you see it."
echo "      (no absolute threshold quoted: h/u scales as ~1/G, so h at G=4.5"
echo "       is ~18x smaller than the 4e-7 seen in the G=0.25 runs)"
echo

# MPI launch. On a SLURM-managed HPC (e.g. l1.gsc1 / cl_tmp), a bare
# `mpiexec -n N` often HANGS before the code starts, because it can't find
# the process manager / PMI from the batch allocation. srun is SLURM-native
# and inherits the allocation's task layout, so prefer it when available.
# Stage 1 (the serial EVP) is unaffected either way, which is why a hang here
# shows up as "EVP finished, IVP produced no output" — the exact symptom seen.
if command -v srun >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "  launcher: srun (SLURM-native)"

    # ── SLURM memory-variable collision ──────────────────────────────────────
    # srun aborts before starting anything if more than one of these is set:
    #   srun: fatal: SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and
    #         SLURM_MEM_PER_NODE are mutually exclusive.
    # `#SBATCH --mem=128G` sets SLURM_MEM_PER_NODE, but a cluster-wide
    # DefMemPerCPU in slurm.conf ALSO injects SLURM_MEM_PER_CPU into the job
    # environment, and srun then refuses to choose between them. The batch
    # allocation already holds the memory we asked for, so dropping the extra
    # variables changes nothing about what the job gets -- it only removes the
    # ambiguity srun is complaining about.
    echo "    memory vars in the job environment:"
    for _v in SLURM_MEM_PER_NODE SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU; do
        echo "      $_v=${!_v:-<unset>}"
    done
    unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU
    echo "    -> kept SLURM_MEM_PER_NODE, unset the other two"

    srun --mpi=pmi2 -n "$NRANKS" python3 spherical-2d-hd-s2s4.py "$INPUT"
else
    echo "  launcher: mpiexec (no SLURM detected — workstation mode)"
    mpiexec -n "$NRANKS" python3 spherical-2d-hd-s2s4.py "$INPUT"
fi

echo
echo "=== done: $(date '+%Y-%m-%d %H:%M:%S') ==="
echo "Next, in order — do NOT skip the gate:"
echo "  1. python3 ../postrun_analyse.py --input $INPUT"
echo "       gate: CLEAN SINGLE MODE: YES, and omega RETROGRADE (Re < 0),"
echo "       matching the seeded omega printed by Stage 1 above."
echo "       Do NOT use the old -0.4374 — that was G=0.25 with m=5."
echo "       Thin-shell estimate at m=3 is omega ~ -2/(m+1) = -0.5, but the"
echo "       EVP is the authority."
echo "  2. python3 ../lsq_from_data.py --input $INPUT --tmin <from step 1> --cheb 30 --svd 0"
echo "       (clean data wants cheb=30/svd=0; the script defaults to the"
echo "        degraded-data values cheb=20/svd=2 and will warn you)"