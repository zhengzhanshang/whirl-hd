#!/usr/bin/env python3
"""
make_figures.py -- figures for the closed-form Rossby-wave recovery paper
(method paper, HD case, G = 0.25, m = 5).

Three figures. LaTeX numbers them in order of appearance in main.tex, which
is NOT the order of the file names:

    paper Fig. 1  <- fig3_degradations.pdf          (Sect. 4, measured)
    paper Fig. 2  <- fig2_noise_eqsubset.pdf        (Sect. 4.5, measured)
    paper Fig. 3  <- fig1_geometry_conditioning.pdf (Sect. 5, analytic)

ANALYTIC (fig1_geometry_conditioning): computed from the projection geometry
every run; nothing is hard-coded.

MEASURED (fig2, fig3): numbers transcribed from the run logs into the block
below. Before submission, replace that block with a parser that reads the
.log files, so the figures cannot drift from the tables in the text.

A&A figure widths: 88 mm single column, 180 mm double column.

    python3 make_figures.py            # writes the three PDFs into ./figures
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MM = 1/25.4
COL1, COL2 = 88*MM, 180*MM
OUT = "figures"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.6, "lines.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

# =====================================================================
# MEASURED RESULTS -- transcribed from run logs. REPLACE WITH A PARSER.
# =====================================================================
# End-to-end LOS recovery, (s2,s4)=(0.19,0.23), B0=7.25 deg, mask 60 deg,
# K=20, r=2. Errors in per cent. None where not run.
LOS_NOISE = [0.001, 0.01, 0.02, 0.05, 0.10]
LOS = {
    "eq1+2+3": {"s2": [2.30, 5.74, 14.53, 72.30, 273.0],
                "s4": [3.29, 13.43, 38.45, 200.40, 760.0]},
    "eq1+2":   {"s2": [0.13, 0.06, 0.07, 0.66, 3.05],
                "s4": [0.22, 0.03, 0.28, 1.44, 5.59]},
    "eq2":     {"s2": [None, 0.29, None, 0.29, 0.13],
                "s4": [None, 0.43, None, 0.41, 0.29]},
}
# 5-seed scatter at noise 0.05 (std dev, per cent)
LOS_ERRBAR = {"eq1+2": {"s2": (0.05, 0.28), "s4": (0.05, 1.12)},
              "eq2":   {"s2": (0.05, 0.12), "s4": (0.05, 0.33)}}

# ---------------------------------------------------------------------
# THREE-PAIR EXTENSION (runs of 2026-09-21).
# (0.19,0.23) is the reference pair and its numbers are UNCHANGED from
# above; (0.29,0.03) and (0.03,0.33) are the two stress pairs that bracket
# the 11:1 ratio range. Errors in per cent throughout.
#
# NOTE ON THE STRESS PAIRS: in each one the LARGE relative errors sit on
# whichever parameter equals 0.03, i.e. they are small-denominator effects,
# not a loss of accuracy. Absolute errors are given alongside so the figures
# can be redrawn either way -- see ABS_* and the ERR_MODE switch.
# ---------------------------------------------------------------------
PAIRS = ["(0.19,0.23)", "(0.29,0.03)", "(0.03,0.33)"]
PAIR_TRUE = {"(0.19,0.23)": (0.19, 0.23),
             "(0.29,0.03)": (0.29, 0.03),
             "(0.03,0.33)": (0.03, 0.33)}
# marker / linestyle per pair; colour still encodes the PARAMETER
# (crimson = s2, cornflowerblue = s4) exactly as in the original figures.
PAIR_STYLE = {"(0.19,0.23)": ("o", "-"),
              "(0.29,0.03)": ("^", "--"),
              "(0.03,0.33)": ("s", ":")}
C_S2, C_S4 = "crimson", "cornflowerblue"

# --- Sect. 4.3  irregular sampling, reference pair only -------------
# noise 0.1 + limb weighting + 60 deg mask, cheb=20, svd=2, single seed.
IRREG_LABEL = ["baseline", "20% dropout", "40% dropout", "gap\n(~3 periods)"]
IRREG_S2 = [0.32, 0.00, 0.00, 0.31]
IRREG_S4 = [0.49, 0.25, 0.11, 0.46]

# --- Sect. 4.1  polar masking, clean data, FIXED K=30 r=0 -------------
# Source: pole_mask_rerun_G025.log (run_pole_mask_G025.sh, 2026-10-19,
# lsq md5 176d096a...). Relative errors in per cent computed from the
# recovered values in the log (|s_rec - s_true|/s_true), not from the
# rounded per-cent column, so no entry is zero. cutoff 90 means no mask.
POLE = {
    "(0.19,0.23)": {"cut": [90, 60, 40, 20],
                    "s2": [0.0237, 0.0263, 0.0279, 0.0537],
                    "s4": [0.0209, 0.0278, 0.0326, 0.1896],
                    "cond": [7.18, 7.24, 8.99, 18.97]},
    "(0.29,0.03)": {"cut": [90, 60, 40, 20],
                    "s2": [0.0255, 0.0234, 0.0252, 0.0462],
                    "s4": [0.3733, 0.3067, 0.3600, 1.9700],
                    "cond": [7.20, 7.26, 8.97, 19.29]},
    "(0.03,0.33)": {"cut": [90, 60, 40, 20],
                    "s2": [0.0667, 0.0767, 0.1000, 0.1900],
                    "s4": [0.0045, 0.0070, 0.0148, 0.0824],
                    "cond": [6.61, 6.72, 8.63, 18.44]},
}

# --- Sect. 4.4  cadence -----------------------------------------------
# noise 0.1 + limb weighting + 60 deg mask, cheb=20, svd=2.
# The SAMPLES PER PERIOD differ between pairs because the wave period does:
# T = 13.71, 13.02, 15.71 for (0.19,0.23), (0.29,0.03), (0.03,0.33).
# The stress pairs are 5-seed means with standard deviations; the reference
# pair is a single seed, so it carries no error bar.
CAD = {
    "(0.19,0.23)": {"pts": [27.4, 13.7, 6.85],
                    "s2": [0.32, 0.00, 5.79], "s2e": None,
                    "s4": [0.49, 0.44, 5.99], "s4e": None},
    "(0.29,0.03)": {"pts": [26.04, 13.02, 6.51],
                    "s2": [0.19, 0.49, 5.59], "s2e": [0.09, 0.16, 0.34],
                    "s4": [7.07, 10.10, 76.29], "s4e": [2.62, 5.04, 11.90]},
    "(0.03,0.33)": {"pts": [31.42, 15.71, 7.86],
                    "s2": [1.51, 1.18, 18.93], "s2e": [0.70, 1.10, 2.05],
                    "s4": [0.43, 0.35, 1.97], "s4e": [0.23, 0.23, 0.67]},
}
# One reference-pair value rounds to 0.00 per cent in the log. On a log axis
# it cannot be drawn, so it is plotted at this floor with a HOLLOW marker.
CAD_FLOOR = 0.02

# Plot relative error (as the published figures do) or absolute error.
# "rel" keeps the degradation figure consistent with the noise figure and
# with the tables; "abs" removes the small-denominator excursions of the
# stress pairs. Switch here rather than editing the plotting code.
ERR_MODE = "rel"               # "rel" or "abs"


def _scale(pair, par, vals):
    """Convert per-cent relative errors to absolute if ERR_MODE == 'abs'."""
    if ERR_MODE == "rel" or vals is None:
        return vals
    true = PAIR_TRUE[pair][0 if par == "s2" else 1]
    return [None if v is None else v/100.0*true for v in vals]


def _ylabel():
    return ("relative error (per cent)" if ERR_MODE == "rel"
            else "absolute error")
# =====================================================================


# Fig. 3 (fig1_geometry_conditioning) colours: one colour per latitude, shared
# by panels (a) and (b). Red is reserved for the B0 = 0 line in panel (a).
LAT_A = (15, 45, 60)
LAT_COLOR = {15: "cornflowerblue", 30: "darksalmon", 45: "mediumseagreen",
             60: "mediumpurple"}
C_B0ZERO = "crimson"


def proj(mu, phi, b0):
    """(r_hat.e, theta_hat.e, phi_hat.e) with observer e = (cos b0, 0, sin b0)."""
    st = np.sqrt(np.maximum(1 - mu**2, 0))
    cb, sb = np.cos(b0), np.sin(b0)
    return (st*np.cos(phi)*cb + mu*sb,
            mu*np.cos(phi)*cb - st*sb,
            -np.sin(phi)*cb)


def design(mu, phi, m, b0):
    pw, pt, pp = proj(mu, phi, b0)
    c, s = np.cos(m*phi), np.sin(m*phi)
    return np.column_stack([2*pw*c, -2*pw*s, 2*pt*c, -2*pt*s, 2*pp*c, -2*pp*s])


# ---------------------------------------------------------------------
# Fig. 1  ANALYTIC -- the B0 = 0 degeneracy and what it does to conditioning
# ---------------------------------------------------------------------
def fig1():
    fig, ax = plt.subplots(1, 2, figsize=(1.2*COL2, 0.58*COL2))

    # (a) ratio of the two projections, NORMALISED to its central-meridian
    #     value. At B0 = 0 this is identically 1 at every longitude and every
    #     latitude -- the two columns are proportional, so the system is
    #     singular. Any nonzero B0 breaks it.
    phi = np.linspace(-np.pi/2 + 2e-2, np.pi/2 - 2e-2, 400)
    # colour = latitude (same colours as panel b); line style = tilt.
    # B0 = 0 stays one bold red line: the three latitudes coincide there.
    for b0deg, ls in ((0.0, "-"), (1.0, ":"), (7.25, "--")):
        b0 = np.deg2rad(b0deg)
        for lat0 in LAT_A:
            mu0 = np.sin(np.deg2rad(lat0))
            pw, pt, _ = proj(mu0, phi, b0)
            r = (pt/pw) / (pt/pw)[len(phi)//2]
            if b0deg == 0:
                ax[0].plot(np.rad2deg(phi), r, ls, color=C_B0ZERO, lw=2.2, zorder=5)
            else:
                ax[0].plot(np.rad2deg(phi), r, ls, color=LAT_COLOR[lat0], lw=1.2,
                           zorder=2)
    ax[0].set_xlabel("longitude from central meridian (deg)", fontsize=9)
    ax[0].set_ylabel(r"$(\hat\theta\cdot\hat e)/(\hat r\cdot\hat e)$,"
                     "\n normalised to central meridian", fontsize=9)
    ax[0].set_ylim(0.05, 1.18)
    ax[0].set_title("(a) three latitudes: 15, 45, 60 deg", loc="left")
    tilt_handles = [plt.Line2D([], [], color=C_B0ZERO, ls="-", lw=2.2, label=r"$B_0=0^\circ$"),
                    plt.Line2D([], [], color="0.3", ls=":", lw=1.2, label=r"$B_0=1^\circ$"),
                    plt.Line2D([], [], color="0.3", ls="--", lw=1.2, label=r"$B_0=7.25^\circ$")]
    lat_handles = [plt.Line2D([], [], color=LAT_COLOR[l], lw=1.6, label=rf"${l}^\circ$")
                   for l in LAT_A]
    # both keys sit side by side in the empty space inside the innermost
    # curve (15 deg at B0 = 7.25 deg), centred on the central meridian
    leg1 = ax[0].legend(handles=tilt_handles, frameon=False, loc="center right",
                        bbox_to_anchor=(0.53, 0.22), fontsize=8, title="tilt",
                        title_fontsize=8, handlelength=2.0, borderpad=0.2,
                        handletextpad=0.4, labelspacing=0.3)
    ax[0].add_artist(leg1)
    ax[0].legend(handles=lat_handles, frameon=False, loc="center left",
                 bbox_to_anchor=(0.53, 0.22), fontsize=8, title="latitude",
                 title_fontsize=8, handlelength=1.8, borderpad=0.2,
                 handletextpad=0.4, labelspacing=0.3)
    ax[0].text(0.5, 0.93,
               r"$B_0=0$: identically flat $\Rightarrow$ columns proportional",
               transform=ax[0].transAxes, ha="center", fontsize=8.7,
               style="italic", color=C_B0ZERO)

    # (b) condition number vs B0
    m, Nphi = 5, 256
    phi_all = 2*np.pi*np.arange(Nphi)/Nphi
    b0s = np.logspace(-3, np.log10(7.25), 60)
    for lat in (15, 30, 45, 60):
        c = LAT_COLOR[lat]
        mu = np.sin(np.deg2rad(lat))
        conds = []
        for b0deg in b0s:
            b0 = np.deg2rad(b0deg)
            pw, _, _ = proj(mu, phi_all, b0)
            conds.append(np.linalg.cond(design(mu, phi_all[pw > 0], m, b0)))
        ax[1].loglog(b0s, conds, color=c, lw=1.4, alpha=0.8,
                     label=rf"${lat}^\circ$")
    ax[1].loglog(b0s, 100/b0s, "k:", lw=0.9, label=r"$100/B_0$")
    ax[1].axvspan(7.25, 10, color="0.9", zorder=0)
    ax[1].text(8.6, 8e1, "beyond\nsolar $|B_0|$", fontsize=8.7, color="gray",
               ha="center", va="center")
    ax[1].set_xlim(1e-3, 10)
    ax[1].set_xlabel(r"$B_0$ (deg)", fontsize=9)
    ax[1].set_ylabel("condition number", fontsize=9)
    ax[1].set_title("(b) conditioning of the inversion", loc="left")
    ax[1].legend(frameon=False, title="latitude (curves coincide)", fontsize=8.7,
                 title_fontsize=8.7, loc="lower left")

    # stacked panels: hspace (not wspace) keeps the (b) title clear of the
    # (a) x-axis label
    fig.subplots_adjust(hspace=0.32)
    fig.savefig(f"{OUT}/fig1_geometry_conditioning.pdf")
    plt.close(fig)
    print("fig1_geometry_conditioning.pdf   [analytic]")


# ---------------------------------------------------------------------
# Fig. 2  MEASURED -- end-to-end error vs LOS noise, by equation subset
# ---------------------------------------------------------------------
def fig2():
    fig, ax = plt.subplots(1, 1, figsize=( 0.55*COL2, 0.45*COL2))

    # Using CSS color names: (s2_color, s4_color, marker)
    style = {
        "eq1+2+3": ("firebrick", "lightcoral", "o"),
        "eq1+2": ("royalblue", "skyblue", "s"),
        "eq2": ("seagreen", "lightgreen", "D"),
    }

    param_config = {
        "s2": ("-", 0, r"$s_2$"),
        "s4": ("--", 1, r"$s_4$"),
    }

    for par, (ls, col_idx, par_label) in param_config.items():
        for name, (*colors, mk) in style.items():
            c = colors[col_idx]  # Retrieves CSS color name string
            x = [n for n, v in zip(LOS_NOISE, LOS[name][par]) if v is not None]
            y = [v for v in LOS[name][par] if v is not None]

            ax.loglog(
                x,
                y,
                linestyle=ls,
                color=c,
                marker=mk,
                ms=3.5,
                label=f"{name} ({par_label})",
            )

            eb = LOS_ERRBAR.get(name, {}).get(par)
            if eb:
                xi, si = eb
                yi = LOS[name][par][LOS_NOISE.index(xi)]
                ax.errorbar([xi], [yi], yerr=[si], color=c, capsize=5, lw=1)

    ax.axhline(1.0, color="0.7", lw=0.6, ls=":")

    ax.set_xlabel(
        "line-of-sight noise (fraction of RMS)", fontsize=10
    )
    ax.set_ylabel("relative error (per cent)", fontsize=10)
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    '''ax.text(
        0.03,
        0.60,
        "1 per cent",
        transform=ax.transAxes,
        fontsize=8.7,
        color="0.5",
    )'''

    fig.tight_layout()
    fig.savefig(f"{OUT}/fig2_noise_eqsubset.pdf")
    plt.close(fig)
    print("fig2_noise_eqsubset.pdf   [MEASURED -- replace with log parser]")


# ---------------------------------------------------------------------
# Fig. 3  MEASURED -- the three degradations that are not noise
#
# EXTENDED to three (s2,s4) pairs in panels (a) and (b). Encoding, chosen
# to stay consistent with the rest of the paper:
#     COLOUR    = parameter   (crimson s2, cornflowerblue s4), as before
#     MARKER/LS = pair        (circle/solid = reference, triangle/dashed and
#                              square/dotted = the two stress pairs)
# Panel (c) is unchanged: irregular sampling was only run on the reference
# pair, and the panel says so.
#
# Panels (a) and (b) use LOG ordinates: the stress pairs span 0.01 to 76
# per cent. The reference-pair cadence value that rounds to 0.00 per cent
# is drawn at CAD_FLOOR with a hollow marker.
# ---------------------------------------------------------------------
def fig3():
    fig, ax = plt.subplots(1, 3, figsize=(1.1*COL2, 0.34*COL2))

    # ---- (a) polar masking -------------------------------------------
    for pair in PAIRS:
        d = POLE[pair]
        mk, ls = PAIR_STYLE[pair]
        ax[0].plot(d["cut"], _scale(pair, "s2", d["s2"]), marker=mk, ls=ls,
                   color=C_S2, ms=3.2, lw=1.0)
        ax[0].plot(d["cut"], _scale(pair, "s4", d["s4"]), marker=mk, ls=ls,
                   color=C_S4, ms=3.2, lw=1.0)
    ax[0].set_yscale("log")
    ax[0].invert_xaxis()
    ax[0].set_xlabel("mask cutoff latitude (deg)", fontsize=8.7)
    ax[0].set_ylabel(_ylabel())
    ax[0].set_title("(a) polar masking\n", loc="left", fontsize=8.7,
                    linespacing=1.25)
    ax[0].margins(y=0.35)

    # condition number on the right axis, all three pairs (they coincide)
    a0 = ax[0].twinx()
    for pair in PAIRS:
        a0.plot(POLE[pair]["cut"], POLE[pair]["cond"], ":", color="darkgray",
                lw=0.9, alpha=0.85)
    a0.set_ylabel("condition number", color="darkgray", fontsize=8.7)
    a0.tick_params(colors="0.55", labelsize=8)
    a0.set_ylim(0, 22)

    # pair legend lives here, drawn in neutral black so it reads as a
    # linestyle key rather than a parameter key
    pair_handles = [plt.Line2D([], [], color="0.25", marker=PAIR_STYLE[p][0],
                               ls=PAIR_STYLE[p][1], ms=3.2, lw=1.0,
                               label=rf"$(s_2,s_4)=${p}")
                    for p in PAIRS]
    ax[0].legend(handles=pair_handles, frameon=False, loc="upper left",
                 fontsize=6.2)

    # ---- (b) cadence --------------------------------------------------
    for pair in PAIRS:
        d = CAD[pair]
        mk, ls = PAIR_STYLE[pair]
        for par, col in (("s2", C_S2), ("s4", C_S4)):
            y = list(_scale(pair, par, d[par]))
            e = _scale(pair, par, d[par + "e"]) if d[par + "e"] else None
            # a measured 0.00 cannot be drawn on a log axis: floor it and
            # mark it hollow so the reader knows it is a limit, not a value
            floored = [i for i, v in enumerate(y) if v == 0.0]
            for i in floored:
                y[i] = CAD_FLOOR if ERR_MODE == "rel" else CAD_FLOOR/100.0
            ax[1].errorbar(d["pts"], y, yerr=e, marker=mk, ls=ls, color=col,
                           ms=3.2, lw=1.0, capsize=2, elinewidth=0.7)
            for i in floored:
                ax[1].plot(d["pts"][i], y[i], marker=mk, ms=4.5, mfc="white",
                           mec=col, mew=0.9, ls="none", zorder=6)
    ax[1].set_yscale("log")
    ax[1].axvline(10, color="0.6", lw=0.8, ls="--")
    # keep this label off the curves: park it just inside the right edge
    # (the x axis is inverted, so "right" is the low-samples end)
    ax[1].text(0.97, 0.03, "suggested\nfloor", transform=ax[1].transAxes,
               fontsize=6.2, color="0.45", ha="right", va="bottom")
    ax[1].invert_xaxis()
    ax[1].set_xlabel("samples per wave period", fontsize=8.7)
    ax[1].set_ylabel(_ylabel(), fontsize=8.7)
    ax[1].set_title("(b) cadence\n", loc="left", fontsize=8.7,
                    linespacing=1.25)
    ax[1].margins(y=0.30)

    # parameter legend lives here
    par_handles = [plt.Line2D([], [], color=C_S2, lw=1.4, label=r"$s_2$"),
                   plt.Line2D([], [], color=C_S4, lw=1.4, label=r"$s_4$"),
                   plt.Line2D([], [], color="0.4", marker="o", ls="none",
                              ms=4.5, mfc="white", mew=0.9,
                              label="rounds to 0.00")]
    ax[1].legend(handles=par_handles, frameon=False, loc="upper left",
                 fontsize=6.2)

    # ---- (c) irregular sampling (reference pair only) -----------------
    xs = np.arange(len(IRREG_LABEL))
    ax[2].bar(xs-0.18, IRREG_S2, 0.36, color=C_S2, label=r"$s_2$")
    ax[2].bar(xs+0.18, IRREG_S4, 0.36, color=C_S4, label=r"$s_4$")
    ax[2].set_xticks(xs)
    ax[2].set_xticklabels(["base", "20%\ndrop", "40%\ndrop", "gap"],
                          fontsize=8.7)
    ax[2].set_ylim(0, 0.72)
    ax[2].set_ylabel("relative error (per cent)", fontsize=8.7)
    # wrapped onto two lines: the one-line version runs past the panel edge
    ax[2].set_title("(c) irregular sampling,\n" + r"$(0.19,0.23)$ only",
                    loc="left", fontsize=8.7, linespacing=1.25)
    ax[2].legend(frameon=False, loc="upper right", fontsize=8)

    fig.subplots_adjust(wspace=0.62)
    fig.savefig(f"{OUT}/fig3_degradations.pdf")
    plt.close(fig)
    print("fig3_degradations.pdf   [MEASURED, 3 pairs, log axes]")


if __name__ == "__main__":
    fig1(); fig2(); fig3()
    print(f"\nwrote 3 figures to {OUT}/  (ERR_MODE = {ERR_MODE!r})")
    print("fig1_geometry_conditioning is analytic and reproducible as-is.")
    print("fig2 and fig3 use transcribed numbers -- swap in a log parser")
    print("before submission so they cannot drift from the tables.")