# spherical-2d-hydro.py
# 2D Spherical Shallow-Water Rossby Wave Simulation — HYDRODYNAMIC (h formulation)
# Physical-domain IVP.  Induction equations omitted; b-fields removed.
# Prognostic height h' (reduced gravity G); no pressure/constraint.
#
# Coordinates: phi in [0, 2pi)  x  mu = cos(theta) in [-1, 1]
# Basis:       RealFourier(phi) x Chebyshev(mu)
#
# Variables (primitive, real): u_theta, u_phi, h   (all prognostic)
#
# Governing equations (primes dropped, normalised, t' = Omega_0 t):
#   (1)  dt u_θ + Ω_d ∂_φ u_θ − 2(1+Ω_d)μ u_φ = G sinθ ∂_μ h
#   (2)  dt u_φ + Ω_d ∂_φ u_φ + C_φ u_θ        = −(G/sinθ) ∂_φ h
#   (5)  dt h  + Ω_d ∂_φ h + (1/sinθ)[ μ u_θ − D u_θ + ∂_φ u_φ ] = 0
#
# with D = (1−μ²) ∂_μ ,  Ω_d = −s₂μ² − s₄μ⁴ ,
#      C_φ = 2(1+Ω_d)μ + (1−μ²)(2s₂μ + 4s₄μ³),  G = ε = gH₀/(Ω₀²R₀²).
#
# Tau placement:
#   eq (5) carries lift(tau_uth1,-1)+lift(tau_uth2,-2) → enforces u_theta(±1)=0.
#   h is prognostic and free at the poles (no BC); no pressure gauge needed.
#
# Usage:
#   python3 spherical-2d-hydro.py input.json [--restart]
#   mpiexec -n 4 python3 spherical-2d-hydro.py input.json

import json
import sys
import numpy as np
import dedalus.public as d3
import logging
from mpi4py import MPI
from scipy.special import lpmv

comm = MPI.COMM_WORLD

# ── Logging ───────────────────────────────────────────────────────────────────
if comm.rank == 0:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s :: %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("logger_info.txt", mode="w"),
        ],
        force=True,
    )
    logger = logging.getLogger(__name__)
    for _log_name in ("solvers", "subsystems"):
        _lg = logging.getLogger(_log_name)
        _lg.setLevel(logging.INFO)
        _lg.propagate = False
        for _h in logging.getLogger().handlers:
            if _h not in _lg.handlers:
                _lg.addHandler(_h)
else:
    logging.basicConfig(level=logging.WARNING, force=True)
    logger = logging.getLogger(__name__)

# ── Arguments ─────────────────────────────────────────────────────────────────
if len(sys.argv) < 2:
    print("Usage: python3 spherical-2d-hydro.py input.json [--restart]")
    sys.exit(1)
restart = "--restart" in sys.argv
params  = json.load(open(sys.argv[1]))

# ── Parameters ────────────────────────────────────────────────────────────────
m_ic      = params["m"]          # dominant azimuthal wavenumber for IC
s2        = params["s2"]         # differential rotation: Omega_d = -s2*mu^2 - s4*mu^4
s4        = params["s4"]
G         = params["G"]          # reduced gravity  G = eps = g H0 / (Omega0^2 R0^2)
gamma     = params["gamma"]      # magnetic param (UNUSED in hydrodynamic case)
alpha_r   = params["alpha_r"]    # normalised alpha (UNUSED in hydrodynamic case)
nu4       = params["nu4"]        # mu-biharmonic hyperviscosity: +nu4 * d^4/dmu^4 on momentum
nu4_h     = params.get("nu4_h", 0.0)  # OPTIONAL mu-biharmonic on the h/continuity eqn (0=off).
                                      # Needed for large-G runs where the instability lives in h.
Nmu       = params["Nmu"]
Nphi      = params["Nphi"]
IC_amp    = params["IC_amp"]
mu_center = params["mu_center"]
mu_width  = params["mu_width"]
stop_sim_time = params["stop_sim_time"]
max_timestep  = params["max_timestep"]
snap_dt       = params["snap_dt"]
max_writes    = params["max_writes"]
run_name      = params["run_name"]
ic_name       = params["ic_name"]

# CFL controls (optional in JSON; safe defaults otherwise)
cfl_safety    = params.get("cfl_safety", 0.3)   # fraction of the CFL limit to use
cfl_cadence   = params.get("cfl_cadence", 10)   # recompute dt every N steps
min_timestep  = params.get("min_timestep", 1e-8)

snapshots_dir = "snapshots"
state_dir     = run_name + "_state"

dealias     = 3/2
dtype       = np.float64
eps_reg     = 0.05           # pole regularisation for 1/sin(theta)
timestepper = d3.SBDF2

# ── Bases ─────────────────────────────────────────────────────────────────────
coords    = d3.CartesianCoordinates("phi", "mu")
dist      = d3.Distributor(coords, dtype=dtype)
phi_basis = d3.RealFourier(coords["phi"], size=Nphi, bounds=(0, 2*np.pi), dealias=dealias)
mu_basis  = d3.Chebyshev(coords["mu"],   size=Nmu,   bounds=(-1, 1),       dealias=dealias)
basis2d   = (phi_basis, mu_basis)

phi_grid, mu_grid = dist.local_grids(phi_basis, mu_basis)

# ── Prognostic fields ─────────────────────────────────────────────────────────
u_theta = dist.Field(name="u_theta", bases=basis2d)
u_phi   = dist.Field(name="u_phi",   bases=basis2d)
h       = dist.Field(name="h",       bases=basis2d)

# ── Auxiliary fields for mu-biharmonic (s = d^2/dmu^2 of each velocity) ────────
s_theta = dist.Field(name="s_theta", bases=basis2d)
s_phi   = dist.Field(name="s_phi",   bases=basis2d)

# Optional h-biharmonic (large-G): aux s_h = d^2 h/dmu^2 + 4 taus, only if nu4_h>0
use_h_hyper = (nu4_h > 0.0)
if use_h_hyper:
    s_h   = dist.Field(name="s_h", bases=basis2d)
    tau_h1 = dist.Field(name="tau_h1", bases=phi_basis)
    tau_h2 = dist.Field(name="tau_h2", bases=phi_basis)
    tau_h3 = dist.Field(name="tau_h3", bases=phi_basis)
    tau_h4 = dist.Field(name="tau_h4", bases=phi_basis)

# ── Tau fields ────────────────────────────────────────────────────────────────
# 4 per velocity component: clamped BCs  u = 0  and  d u/dmu = 0  at mu = ±1.
tau_ut1 = dist.Field(name="tau_ut1", bases=phi_basis)
tau_ut2 = dist.Field(name="tau_ut2", bases=phi_basis)
tau_ut3 = dist.Field(name="tau_ut3", bases=phi_basis)
tau_ut4 = dist.Field(name="tau_ut4", bases=phi_basis)
tau_up1 = dist.Field(name="tau_up1", bases=phi_basis)
tau_up2 = dist.Field(name="tau_up2", bases=phi_basis)
tau_up3 = dist.Field(name="tau_up3", bases=phi_basis)
tau_up4 = dist.Field(name="tau_up4", bases=phi_basis)

# ── Metric / parameter fields (mu-dependent only) ─────────────────────────────
mu_f = dist.Field(name="mu_f", bases=mu_basis)
mu_f["g"] = mu_grid

one_minus_mu2 = dist.Field(name="one_minus_mu2", bases=mu_basis)
one_minus_mu2["g"] = 1.0 - mu_grid**2

st_f = dist.Field(name="st_f", bases=mu_basis)        # sin theta = sqrt(1 - mu^2)
st_f["g"] = np.sqrt(np.maximum(1.0 - mu_grid**2, 0.0))

ist_f = dist.Field(name="ist_f", bases=mu_basis)      # 1/sin theta (regularised)
ist_f["g"] = 1.0 / (st_f["g"] + eps_reg)

Omega_d_f = dist.Field(name="Omega_d_f", bases=mu_basis)
Omega_d_f["g"] = -s2 * mu_grid**2 - s4 * mu_grid**4

# Coriolis coefficient in u_theta eq: 2(1+Omega_d)mu
corf_th_f = dist.Field(name="corf_th_f", bases=mu_basis)
corf_th_f["g"] = 2.0 * mu_grid * (1.0 + Omega_d_f["g"])

# Coriolis coefficient in u_phi eq: 2(1+Omega_d)mu + (1-mu^2)(2s2*mu + 4s4*mu^3)
corf_phi_f = dist.Field(name="corf_phi_f", bases=mu_basis)
corf_phi_f["g"] = (
    2.0 * mu_grid * (1.0 + Omega_d_f["g"])
    + (1.0 - mu_grid**2) * (2.0*s2*mu_grid + 4.0*s4*mu_grid**3)
)

# ── Operators ─────────────────────────────────────────────────────────────────
dphi = lambda f: d3.Differentiate(f, coords["phi"])
dmu  = lambda f: d3.Differentiate(f, coords["mu"])
D_op = lambda f: one_minus_mu2 * dmu(f)    # D = (1 - mu^2) d/dmu

lift_basis = mu_basis.derivative_basis(2)   # 2nd-order equations → lift at -1, -2
lift = lambda A, n: d3.Lift(A, lift_basis, n)

# ── Problem ───────────────────────────────────────────────────────────────────
_prob_vars = [u_theta, u_phi, h, s_theta, s_phi,
              tau_ut1, tau_ut2, tau_ut3, tau_ut4,
              tau_up1, tau_up2, tau_up3, tau_up4]
if use_h_hyper:
    _prob_vars += [s_h, tau_h1, tau_h2, tau_h3, tau_h4]
problem = d3.IVP(_prob_vars, namespace=locals())

# ── Eq (1): u_theta momentum (+ mu-biharmonic nu4 * d^4 u_theta/dmu^4) ─────────
# dt u_θ + Ω_d ∂_φ u_θ − 2(1+Ω_d)μ u_φ = G sinθ ∂_μ h − ν4 ∂⁴_μ u_θ
problem.add_equation(
    "dt(u_theta) + Omega_d_f*dphi(u_theta) - corf_th_f*u_phi "
    "- G*st_f*dmu(h) + nu4*dmu(dmu(s_theta)) "
    "+ lift(tau_ut1,-1) + lift(tau_ut2,-2) = 0"
)

# ── Eq (2): u_phi momentum (+ mu-biharmonic nu4 * d^4 u_phi/dmu^4) ─────────────
# dt u_φ + Ω_d ∂_φ u_φ + C_φ u_θ = −(G/sinθ) ∂_φ h − ν4 ∂⁴_μ u_φ
problem.add_equation(
    "dt(u_phi) + Omega_d_f*dphi(u_phi) + corf_phi_f*u_theta "
    "+ G*ist_f*dphi(h) + nu4*dmu(dmu(s_phi)) "
    "+ lift(tau_up1,-1) + lift(tau_up2,-2) = 0"
)

# ── Auxiliary definitions: s = d^2/dmu^2 of each velocity ──────────────────────
problem.add_equation(
    "s_theta - dmu(dmu(u_theta)) + lift(tau_ut3,-1) + lift(tau_ut4,-2) = 0"
)
problem.add_equation(
    "s_phi   - dmu(dmu(u_phi))   + lift(tau_up3,-1) + lift(tau_up4,-2) = 0"
)

# ── Eq (5): continuity (h prognostic) ─────────────────────────────────────────
# dt h + Ω_d ∂_φ h + (1/sinθ)[ μ u_θ − D u_θ + ∂_φ u_φ ] = 0
# Optional mu-biharmonic filter on h (+nu4_h ∂⁴_μ h) for large-G stabilisation,
# via aux s_h = ∂²_μ h. h has NO physical Dirichlet BC at the poles, so the
# filter uses NEUMANN/free conditions (∂μ h = 0 and ∂³μ h = ∂μ s_h = 0), which
# damp small scales without pinning h's pole value or distorting the physics.
if use_h_hyper:
    problem.add_equation(
        "dt(h) + Omega_d_f*dphi(h) "
        "+ ist_f*(mu_f*u_theta - D_op(u_theta) + dphi(u_phi)) "
        "+ nu4_h*dmu(dmu(s_h)) + lift(tau_h1,-1) + lift(tau_h2,-2) = 0"
    )
    problem.add_equation(
        "s_h - dmu(dmu(h)) + lift(tau_h3,-1) + lift(tau_h4,-2) = 0"
    )
else:
    problem.add_equation(
        "dt(h) + Omega_d_f*dphi(h) "
        "+ ist_f*(mu_f*u_theta - D_op(u_theta) + dphi(u_phi)) = 0"
    )

# ── Boundary conditions: clamped, u = 0 and du/dmu = 0 at mu = ±1 ─────────────
problem.add_equation("u_theta(mu=-1) = 0")
problem.add_equation("u_theta(mu=1)  = 0")
problem.add_equation("dmu(u_theta)(mu=-1) = 0")
problem.add_equation("dmu(u_theta)(mu=1)  = 0")
problem.add_equation("u_phi(mu=-1) = 0")
problem.add_equation("u_phi(mu=1)  = 0")
problem.add_equation("dmu(u_phi)(mu=-1) = 0")
problem.add_equation("dmu(u_phi)(mu=1)  = 0")

# h-filter BCs (Neumann/free): only when the h-biharmonic is active
if use_h_hyper:
    problem.add_equation("dmu(h)(mu=-1) = 0")
    problem.add_equation("dmu(h)(mu=1)  = 0")
    problem.add_equation("dmu(s_h)(mu=-1) = 0")
    problem.add_equation("dmu(s_h)(mu=1)  = 0")

# ── Solver ────────────────────────────────────────────────────────────────────
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# ── Initial conditions ────────────────────────────────────────────────────────
if ic_name == "none":
    # Balanced, non-divergent velocity seed from a streamfunction Psi.
    #   Psi(mu, phi) = profile(mu) * cos(m phi)
    # Velocity is u = k x grad(Psi)  (non-divergent by construction):
    #   u_theta = -(1/sin theta) d Psi/d phi
    #   u_phi   = - sin theta      d Psi/d mu
    # Height starts at h = 0.  Profile is normalised so max|Psi| = IC_amp.
    #
    # Two profiles are provided — comment out one, keep the other.

    # ---- (A) Associated Legendre P_l^m, sectoral l = m  (ACTIVE) ─────────────
    #   P_l^m vanishes at mu = ±1 for l >= m >= 1, so u vanishes at the poles.
    l_deg = int(params.get("l_deg", m_ic))     # sectoral mode: l = m
    prof = lpmv(m_ic, l_deg, mu_grid)          # (1, Nmu)
    prof_label = f"Legendre P_l^m (l={l_deg}, m={m_ic})"

    # ---- (B) Boundary-vanishing Gaussian  (comment (A), uncomment this) ──────
    #   raw bump minus the endpoint line, so profile(±1) = 0 exactly.
    #   (mu_width is the e-folding scale on (mu-mu_center)^2, = sqrt(2)*sigma.)
    # Graw = np.exp(-((mu_grid - mu_center)**2) / mu_width**2)
    # G_lo = np.exp(-((-1.0 - mu_center)**2) / mu_width**2)      # at mu = -1
    # G_hi = np.exp(-(( 1.0 - mu_center)**2) / mu_width**2)      # at mu = +1
    # baseline = G_lo + (G_hi - G_lo) * (mu_grid + 1.0) / 2.0    # endpoint line
    # prof = Graw - baseline                                      # zero at mu = ±1
    # prof_label = f"Gaussian (mu0={mu_center}, width={mu_width})"

    # normalise profile so that max|Psi| = IC_amp  (global max, MPI-safe)
    local_max = float(np.max(np.abs(prof))) if prof.size else 0.0
    prof_max  = comm.allreduce(local_max, op=MPI.MAX)
    if prof_max == 0.0:
        if comm.rank == 0:
            logger.error("IC profile is identically zero (check l>=m or mu_width).")
        sys.exit(1)
    prof = IC_amp * prof / prof_max

    # streamfunction field
    Psi = dist.Field(name="Psi", bases=basis2d)
    Psi["g"] = prof * np.cos(m_ic * phi_grid)

    # velocity from Psi — derivatives via Dedalus operators (MPI-safe);
    # metric factors are local, pointwise, and exact (sin theta > 0 on the grid).
    dPsi_dphi = dphi(Psi).evaluate(); dPsi_dphi.change_scales(1)
    dPsi_dmu  = dmu(Psi).evaluate();  dPsi_dmu.change_scales(1)
    st_grid = np.sqrt(np.maximum(1.0 - mu_grid**2, 0.0))

    u_theta["g"] = -(1.0 / st_grid) * dPsi_dphi["g"]
    u_phi["g"]   = -st_grid * dPsi_dmu["g"]
    h["g"]       = 0.0

    if comm.rank == 0:
        logger.info("IC: balanced streamfunction  %s  amp=%.4e", prof_label, IC_amp)
        logger.info("IC: u = k x grad(Psi)  (non-divergent);  h = 0")
elif ic_name == "evp":
    # ── EVP eigenvector IC ───────────────────────────────────────────────────
    # Seed exactly one eigenmode of the SAME discrete operator this IVP will
    # integrate, so nothing else is excited and h starts already balanced
    # (no geostrophic adjustment transient).
    #
    # Convention (sw_model.py): fields ~ exp(i m phi) exp(-i omega t), so the
    # real-space field is Re[ u_hat(mu) exp(i m phi) ].
    import h5py
    ic_path = params.get("evp_ic_file", "evp_ic.h5")
    with h5py.File(ic_path, "r") as f:
        mu_e   = np.asarray(f["mu"][:]).ravel()
        uth_e  = np.asarray(f["u_theta_hat"][:]).ravel()
        uph_e  = np.asarray(f["u_phi_hat"][:]).ravel()
        h_e    = np.asarray(f["h_hat"][:]).ravel()
        a      = dict(f.attrs)

    # Refuse a stale IC rather than silently integrating the wrong mode.
    for key, want in (("m", m_ic), ("Nmu", Nmu), ("G", G),
                      ("s2", s2), ("s4", s4), ("nu4", nu4)):
        if key not in a:
            continue
        got = a[key]
        if abs(float(got) - float(want)) > 1e-12 * max(1.0, abs(float(want))):
            if comm.rank == 0:
                logger.error("IC %s mismatch: %s has %s=%s, input.json has %s",
                             ic_path, ic_path, key, got, want)
            sys.exit(1)

    # The EVP and IVP share the same Chebyshev(Nmu) grid on (-1,1), but this
    # rank may hold only part of it. Interpolate onto the LOCAL mu so the
    # branch is correct however Dedalus distributes the axis.
    mu_l  = np.asarray(mu_grid).ravel()
    order = np.argsort(mu_e)
    def _local(z):
        re = np.interp(mu_l, mu_e[order], np.real(z)[order])
        im = np.interp(mu_l, mu_e[order], np.imag(z)[order])
        return re, im

    uth_re, uth_im = _local(uth_e)
    uph_re, uph_im = _local(uph_e)
    h_re,   h_im   = _local(h_e)

    cos_mp = np.cos(m_ic * phi_grid)          # (Nphi_local, 1)
    sin_mp = np.sin(m_ic * phi_grid)

    u_theta["g"] = uth_re[None, :] * cos_mp - uth_im[None, :] * sin_mp
    u_phi["g"]   = uph_re[None, :] * cos_mp - uph_im[None, :] * sin_mp
    h["g"]       = h_re[None, :]   * cos_mp - h_im[None, :]   * sin_mp

    if comm.rank == 0:
        logger.info("IC: EVP eigenvector from %s", ic_path)
        logger.info("IC: omega = %+.6f %+.6fj   (retrograde => Re < 0)",
                    a.get("omega_real", float("nan")),
                    a.get("omega_imag", float("nan")))
        logger.info("IC: max|u_theta|=%.4e  max|u_phi|=%.4e  max|h|=%.4e",
                    np.max(np.abs(u_theta["g"])) if u_theta["g"].size else 0.0,
                    np.max(np.abs(u_phi["g"]))   if u_phi["g"].size else 0.0,
                    np.max(np.abs(h["g"]))       if h["g"].size else 0.0)
        logger.info("IC: h is NOT zero here — the eigenvector carries its own "
                    "balanced height. A <h2> spike in the first writes means "
                    "this branch did not run.")

else:
    ic_file = ic_name + "_state/" + ic_name + "_state_s1.h5"
    write, initial_time = solver.load_state(ic_file)
    if restart:
        solver.sim_time      = initial_time
        solver.stop_sim_time = initial_time + stop_sim_time
    if comm.rank == 0:
        logger.info("Loaded state from %s  t=%.4f", ic_file, solver.sim_time)

# ── Analysis ──────────────────────────────────────────────────────────────────
snapshots = solver.evaluator.add_file_handler(
    snapshots_dir, sim_dt=snap_dt, max_writes=max_writes, mode="overwrite"
)
for fld, name in ((u_theta, "u_theta"), (u_phi, "u_phi"), (h, "h")):
    snapshots.add_task(fld, name=name)

# Relative vorticity: zeta = ist * dphi(u_theta) + dmu(st * u_phi)
snapshots.add_task(ist_f*dphi(u_theta) + dmu(st_f*u_phi), name="vorticity")
# Absolute vorticity: planetary (2mu) + relative
snapshots.add_task(2*mu_f + ist_f*dphi(u_theta) + dmu(st_f*u_phi), name="abs_vorticity")
# Kinetic energy density
snapshots.add_task(u_theta**2 + u_phi**2, name="KE_density")

# Final-state checkpoint (for restart)
final_state = solver.evaluator.add_file_handler(
    state_dir, sim_dt=stop_sim_time - max_timestep, mode="overwrite"
)
final_state.add_tasks(solver.state)

# ── Flow diagnostics ──────────────────────────────────────────────────────────
flow = d3.GlobalFlowProperty(solver, cadence=10)
flow.add_property(u_theta**2 + u_phi**2, name="KE")
flow.add_property(h**2, name="h2")

# ── CFL time-step control ─────────────────────────────────────────────────────
# The fastest signal is the shallow-water GRAVITY wave, c_g = sqrt(G) (normalised);
# the flow is a small linear perturbation, so advection is sub-dominant but is
# included for generality. The velocity components (u_theta, u_phi) are PHYSICAL
# spherical components, not conjugate to the (phi, mu) grid, so the metric enters
# explicitly:
#   phi-cell crossing rate: ist * (c_g + |u_phi|)   / dphi     (ist = 1/(sinθ+eps))
#   mu -cell crossing rate:  st * (c_g + |u_theta|) / dmu       (st  = sinθ)
# dt_CFL = cfl_safety / max_over_grid(sum of the two rates), capped by max_timestep.
cg   = float(np.sqrt(max(G, 1e-30)))                     # gravity-wave speed sqrt(G)
dphi = 2.0 * np.pi / Nphi
# analytic Gauss-Chebyshev (roots) mu grid -> identical on every rank (MPI-safe)
_mu_cheb = np.cos(np.pi * (2*np.arange(Nmu) + 1) / (2*Nmu))
_order   = np.argsort(_mu_cheb)
_mu_srt  = _mu_cheb[_order]
_dmu     = np.gradient(_mu_srt)                          # per-point mu spacing
_st_cheb = np.sqrt(np.maximum(1.0 - _mu_srt**2, 0.0))    # sin(theta)
_ist_cheb = 1.0 / (_st_cheb + eps_reg)                   # regularised 1/sin(theta)

def compute_timestep():
    """Adaptive CFL time step (gravity wave + advection), MPI-safe, capped."""
    # global max physical flow speeds (0 contribution if a rank holds no grid)
    umax_phi = comm.allreduce(
        float(np.max(np.abs(u_phi["g"]))) if u_phi["g"].size else 0.0, op=MPI.MAX)
    umax_th  = comm.allreduce(
        float(np.max(np.abs(u_theta["g"]))) if u_theta["g"].size else 0.0, op=MPI.MAX)
    inv_dt = (_ist_cheb * (cg + umax_phi) / dphi
              + _st_cheb * (cg + umax_th) / np.abs(_dmu))
    inv_dt_max = float(np.max(inv_dt))
    dt = cfl_safety / inv_dt_max if inv_dt_max > 0 else max_timestep
    return float(np.clip(dt, min_timestep, max_timestep))

# ── Main loop ─────────────────────────────────────────────────────────────────
try:
    # Initial CFL dt — compute_timestep() is an MPI collective (allreduce),
    # so ALL ranks must call it here, not just rank 0.
    timestep = compute_timestep()
    if comm.rank == 0:
        logger.info("Starting main loop (HYDRODYNAMIC, shallow-water h)")
        logger.info(
            "Params: m=%d  s2=%.4f  s4=%.4f  G=%.4g  nu4=%.2e  nu4_h=%.2e%s  "
            "(gamma, alpha_r ignored: %.1e, %.1e)",
            m_ic, s2, s4, G, nu4, nu4_h,
            "" if use_h_hyper else " (h-filter OFF)", gamma, alpha_r,
        )
        logger.info(
            "CFL: c_g=sqrt(G)=%.4g  safety=%.2f  cadence=%d  ->  initial dt=%.3e  (cap=%.3e)",
            cg, cfl_safety, cfl_cadence, timestep, max_timestep,
        )
    while solver.proceed:
        # Adaptive gravity-wave CFL; recomputed every cfl_cadence steps
        # (called by all ranks together — the collective stays synchronised).
        if solver.iteration % cfl_cadence == 0:
            timestep = compute_timestep()
        solver.step(timestep)
        if solver.iteration % 1000 == 0:
            ke  = flow.volume_integral("KE")   # MPI collective — all ranks call
            h2  = flow.volume_integral("h2")
            if comm.rank == 0:
                logger.info(
                    "Iter=%i  t=%.4f  dt=%.2e  KE=%.4e  <h²>=%.4e",
                    solver.iteration, solver.sim_time, timestep, ke, h2,
                )
except Exception:
    logger.error("Exception raised, triggering end of main loop.")
    raise
finally:
    solver.log_stats()

if comm.rank == 0:
    logger.info("Simulation completed successfully.")