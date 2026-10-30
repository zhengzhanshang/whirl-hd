# whirl-hd

Code for **WHIRL** (Whole-field Hidden-parameter Inversion of Rossby waves by
Least squares), the closed-form solver described in the paper

> **Shang, Z.-Z., et al. (in review).**
> *WHIRL: Closed-form recovery of solar differential rotation from equatorial Rossby waves.*
> Journal, volume, page. DOI: `[placeholder]` · arXiv: `[placeholder]`

The paper is the full description of the method, the tests and the results. This
repository holds the code that produced every number and figure in it, so that
anyone can check them or apply the method to their own data.

## What the method does

Rossby waves on a differentially rotating sphere carry the rotation profile in
their flow pattern. We write the surface rotation as
Ω(μ) = 1 − s₂μ² − s₄μ⁴ (μ = cos θ) and evolve the waves with the 2D spherical
shallow-water equations in [Dedalus](https://dedalus-project.org).

The key property is that these equations are **exactly linear in s₂ and s₄**.
The residual of the equations, evaluated on the wave fields, therefore splits as
R = c₀ + s₂ c₂ + s₄ c₄, and recovering (s₂, s₄) is a linear least-squares solve:
no optimisation, no training, seconds on a laptop. The paper tests this on clean
simulation output, under observational degradations (polar masking,
centre-to-limb noise, irregular sampling, cadence), through a synthetic
line-of-sight (Dopplergram) inversion, and at two values of the gravity
parameter G (0.25 and 4.5). It also derives an *h-free* form that uses the two
velocity components only, in which G cancels exactly.

## Repository layout

```
simulation/
  G0.25_m5/                 forward model at G = 0.25, m = 5 (Sects. 3-5)
    spherical-2d-hd-s2s4.py   initial-value problem (time integration)
    evp_spherical_hd.py       eigenvalue problem; writes the eigenvector initial condition
    sw_model.py               the equations, shared by both solvers
    cases/s2=..._s4=.../      input.json and run.sh of the six (s2, s4) runs
  G4.5_m5/                  forward model at G = 4.5, m = 5 (Sect. 6)
    (same three solver files) run.sh, input_template.json, pairs.txt, make_cases.sh
analysis/
  lsq_from_data.py          closed-form recovery from u_theta, u_phi and h (all equations)
  lsq_hfree.py              closed-form recovery from u_theta and u_phi only (h-free)
  los_inversion.py          synthetic Dopplergrams and their inversion to the three fields
postprocess/
  postrun_analyse.py        post-run check: single clean mode, its frequency and decay rate
  fit_growth.py             growth / decay rate from the kinetic-energy log
figures/
  make_figures.py           the three figures of the paper
reproduce/
  G0.25/  G4.5/  both/      batch scripts that produced the paper's numbers (see below)
tools/
  check_versions.sh         finds every copy of these scripts on a machine and checks versions
MANIFEST.txt, checksums.md5 where every file came from, with its MD5 checksum
```

## Installation

```bash
conda env create -f environment.yml
conda activate closed-form-hd
CC=mpicc pip install --no-cache --no-build-isolation dedalus==3.0.4
```

The paper used Python 3.12, Dedalus 3.0.4, NumPy 2.2.6, SciPy 1.16.2, h5py 3.14,
mpi4py 4.1, OpenMPI 5.0.8, Matplotlib 3.10.8 and SymPy 1.14 (SymPy only for the
self-tests).

## Quick check (no simulation data needed)

Two self-tests check the core algebra in a few seconds:

```bash
python3 analysis/lsq_hfree.py --selftest       # re-derives the h-free coefficients with SymPy
python3 analysis/los_inversion.py --selftest   # projection geometry and inversion round trip
```

Both should end with `SELFTEST PASSED`. The first one also confirms that the
h-free relation is linear in (s₂, s₄) and that G cancels.

## Running a simulation

Each run needs its own folder holding `input.json`, `run.sh` and the three solver
files. For one G = 0.25 case:

```bash
cd simulation/G0.25_m5/cases/s2=0.19_s4=0.23
cp ../../{spherical-2d-hd-s2s4.py,evp_spherical_hd.py,sw_model.py} .
sbatch run.sh          # on a SLURM cluster; or: bash run.sh on a workstation with MPI
```

`run.sh` first solves the eigenvalue problem (serial, about 30 s) to build an
initial condition that excites one Rossby mode only, then integrates the
initial-value problem on 48 MPI ranks (about 3.5 h to t = 200). For G = 4.5,
`bash simulation/G4.5_m5/make_cases.sh` builds all seven case folders from the
template.

Then check the run and recover (s₂, s₄):

```bash
python3 ../../../../postprocess/postrun_analyse.py snapshots/*.h5 --input input.json
python3 ../../../../analysis/lsq_from_data.py --input input.json --tmin 10 --cheb 30 --svd 0
```

`postrun_analyse.py` must report `CLEAN SINGLE MODE: YES` with a retrograde
frequency matching the one the eigenvalue solver printed. For clean data use
`--cheb 30 --svd 0`; the script's defaults (`--cheb 20 --svd 2`) are meant for
noisy data and it warns if they are used on a clean run.

## Reproducing the paper

All drivers run from the folder that holds the case folders (`s2=..._s4=.../`),
and take the location of the analysis scripts from environment variables:

```bash
cd simulation/G0.25_m5/cases
export LSQ=$PWD/../../../analysis/lsq_from_data.py \
       LOS=$PWD/../../../analysis/los_inversion.py \
       HFREE=$PWD/../../../analysis/lsq_hfree.py
bash ../../../reproduce/G0.25/run_table1_G025.sh
```

Each driver writes one log with the exact command, settings and exit code of
every run, and prints a summary at the end.

| Paper item | G | Driver (in `reproduce/`) |
|---|---|---|
| Table 1 (six pairs, clean data) | 0.25 | `G0.25/run_table1_G025.sh` |
| Sect. 2.5 (K = 20 on clean data); Sects. 4.2-4.4, reference pair; Fig. 1c | 0.25 | `G0.25/run_verify_numbers_G025.sh` |
| Sect. 4.1, Fig. 1a (polar masking, three pairs) | 0.25 | `G0.25/run_pole_mask_G025.sh` |
| Sect. 4.4, Fig. 1b (cadence, two stress pairs) | 0.25 | `G0.25/run_sec44_cadence.sh` |
| Sects. 4.5 and 5, Fig. 2, Table 2 (LOS chain, reference pair) | 0.25 | `G4.5/run_gdep_G4.5_m5.sh` with `STAGES=C` |
| Table 2 (LOS chain, two stress pairs) | 0.25 | `G0.25/run_sec5_los.sh` |
| Sect. 5.5, Table 4 (h-free recovery on the LOS fields) | both | `both/run_hfree_los.sh` |
| Sect. 6, Table 3 and polar mask at G = 4.5 | 4.5 | `G4.5/run_gdep_G4.5_m5.sh` (stages A, B) |
| Sect. 6, LOS chain at G = 4.5 | 4.5 | `G4.5/run_gdep_G4.5_m5.sh` (stage C) |
| The (0.15, 0.15) case at G = 4.5 | 4.5 | `G4.5/make_case_s015_s015.sh` |
| Fig. 3 (geometry and conditioning) | - | `figures/make_figures.py` (analytic) |
| Appendix B (h-free coefficients) | - | `analysis/lsq_hfree.py --selftest` |

`run_gdep_G4.5_m5.sh` is general: run from the G = 0.25 case folder with
`STAGES=C REF=s2=0.19_s4=0.23`, it reproduces the G = 0.25 line-of-sight
results, even though its banner says G = 4.5.

`make_figures.py` draws Figs. 1 and 2 from numbers transcribed from the logs
above (Fig. 3 is computed from the geometry). Run it from the `figures/` folder;
it writes the PDFs into `figures/figures/`.

## Provenance and checksums

Every script is a byte-identical copy of the file that produced the paper's
numbers, collected from the two machines that ran them with
`tools/check_versions.sh`. `MANIFEST.txt` lists each file's MD5 checksum and
original location. To check your copy:

```bash
md5sum -c checksums.md5
```

Because the files are kept unchanged, some contain settings of our own systems:
the SLURM partition and conda path in the `run.sh` files, and absolute paths in
comments. Adjust these for your machine; they do not affect the results.

Four files differ from that rule, and `MANIFEST.txt` marks them:

* `reproduce/G0.25/run_sec44_cadence.sh` and `run_sec5_los.sh` were recovered
  from the record of the session in which they were written (21 Sept 2026),
  because no copy survived on disk. They first ran with an earlier version of
  `lsq_from_data.py`. On 30 Oct 2026 `reproduce/G0.25/rerun_stress_pairs_pinned.sh`
  re-ran both with the archived code and reproduced the paper's values.
* `simulation/G4.5_m5/make_cases.sh` and `rerun_stress_pairs_pinned.sh` were
  added for this repository.

Two small naming notes: the `run.sh` files point to `postrun_analyze.py`, which
is `postprocess/postrun_analyse.py` here, and the G = 0.25 (0.15, 0.17) run
writes snapshots every 0.1 instead of 0.5, so `run_table1_G025.sh` subsamples
it to the same time grid as the other five pairs.

## Data

The simulation snapshots (up to about 1.6 GB per run)
are too large for GitHub. They will be archived on Zenodo: `[placeholder DOI]`.
Until then, they are available from the author on request.

## License

The code in this repository is released under the [MIT License](LICENSE).
It uses [Dedalus](https://dedalus-project.org), which is distributed separately
under its own licence (GPL-3.0).

## Citation

If you use this code, please cite the paper above. `[BibTeX placeholder, to be
added on publication]`

## Contact

Zheng-Zhan Shang, University of Graz ·
zhengzhan.shang@uni-graz.at
