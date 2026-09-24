"""Render the report figures from the experiment results.

Reads data/results/<id>.csv (written by the executables in experiments/) and writes
docs/figures/<id>.svg. SVG output is deterministic (fixed hash salt, no timestamp), so an
unchanged result gives an unchanged file.

Usage: python3 tools/make_figures.py [--results data/results] [--out docs/figures] [id ...]
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

# Chart palette: categorical slots in fixed order, recessive chrome, text in ink colors.
SERIES = ["#2a78d6", "#eb6834"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

plt.rcParams.update({
    "svg.hashsalt": "riskengine",
    "svg.fonttype": "none",
    "font.size": 10,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK_SECONDARY,
    "axes.titlecolor": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_SECONDARY,
    "ytick.labelcolor": INK_SECONDARY,
    "legend.frameon": False,
    "legend.labelcolor": INK_SECONDARY,
})


def read_csv(path):
    """Numeric CSV as a dict of columns."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {key: [float(row[key]) for row in rows] for key in rows[0]}


def fd_vcurve(data, meta, out):
    eps = meta["parameters"]["machine_epsilon"]
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    h = data["h_rel"]
    series = [("delta", data["delta_abs_error"], "1/3"), ("gamma", data["gamma_abs_error"], "1/4")]
    for (name, err, power), color in zip(series, SERIES):
        num, den = map(int, power.split("/"))
        h_opt = eps ** (num / den)
        ax.loglog(h, err, color=color, linewidth=2, label=rf"{name}: optimum $h \approx \varepsilon^{{{power}}}$")
        ax.axvline(h_opt, color=color, linewidth=1, linestyle=(0, (4, 3)), alpha=0.6)
        ax.annotate(name, (h[-1], err[-1]), xytext=(6, 0), textcoords="offset points",
                    va="center", color=INK_SECONDARY)
    ax.text(0.02, 0.04, "← rounding: error ∝ ε/h (delta), ε/h² (gamma)", transform=ax.transAxes, color=MUTED)
    ax.text(0.98, 0.04, "truncation: error ∝ h² →", transform=ax.transAxes, color=MUTED, ha="right")
    ax.set_xlabel("relative bump h")
    ax.set_ylabel("absolute error vs closed form")
    ax.set_title("Finite-difference Greeks: error is V-shaped in the bump size", loc="left")
    ax.legend(loc="upper right", frameon=True, facecolor=SURFACE, edgecolor="none", framealpha=1.0)
    fig.tight_layout()
    fig.savefig(out, metadata={"Date": None})
    plt.close(fig)


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def iv_roundtrip(path, meta, out):
    rows = [r for r in read_rows(path) if r["status"] == "ok"]
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    lo, hi = 1e-17, 1.0
    ax.loglog([lo, hi], [lo, hi], color=MUTED, linewidth=1, linestyle=(0, (4, 3)))
    ax.annotate("error = noise bound", (1e-9, 1e-9), xytext=(-8, 10), textcoords="offset points",
                color=MUTED, ha="right")
    for (name, itm), color in zip([("out-of-the-money quote", "0"), ("in-the-money quote", "1")], SERIES):
        pts = [(float(r["noise_bound"]), float(r["rel_vol_error"])) for r in rows if r["in_the_money"] == itm]
        shown = [(x, y) for x, y in pts if y > 0.0]
        exact = len(pts) - len(shown)
        ax.scatter([x for x, _ in shown], [y for _, y in shown], s=30, color=color, edgecolors=SURFACE,
                   linewidths=0.8, label=f"{name} ({len(pts)} solved, {exact} exact to the bit, not shown)")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(r"noise bound $\varepsilon\,(1 + (1+d^2)(a+b)/(\sigma\,\mathrm{vega}))$")
    ax.set_ylabel("relative error of the recovered vol")
    ax.set_title("Implied vol round trip: the error never exceeds the quote's own noise", loc="left")
    ax.legend(loc="upper left", frameon=True, facecolor=SURFACE, edgecolor="none", framealpha=1.0)
    fig.tight_layout()
    fig.savefig(out, metadata={"Date": None})
    plt.close(fig)


def mc_convergence(path, meta, out):
    data = read_csv(path)
    n, se, err = data["paths"], data["std_error"], data["abs_error"]
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.loglog(n, se, color=SERIES[0], linewidth=2, marker="o", markersize=5, label="standard error")
    ax.loglog(n, err, color=SERIES[1], linewidth=0, marker="o", markersize=6, markeredgecolor=SURFACE,
              label="|estimate − Black-Scholes|")
    # Reference slope, drawn a factor 3 below the data so it does not hide the SE line.
    ref = [se[0] / 3 * (x / n[0]) ** -0.5 for x in (n[0], n[-1])]
    ax.loglog([n[0], n[-1]], ref, color=MUTED, linewidth=1, linestyle=(0, (4, 3)))
    ax.annotate(r"reference slope $-1/2$", (n[-1], ref[-1]), xytext=(0, -16), textcoords="offset points",
                color=MUTED, ha="right")
    ax.set_xlabel("paths N")
    ax.set_ylabel("price error")
    ax.set_title(r"Monte Carlo on an ATM call: error and SE fall as $N^{-1/2}$", loc="left")
    ax.legend(loc="upper right", frameon=True, facecolor=SURFACE, edgecolor="none", framealpha=1.0)
    fig.tight_layout()
    fig.savefig(out, metadata={"Date": None})
    plt.close(fig)


def mc_coverage(path, meta, out):
    import math

    rows = read_rows(path)
    panels = [("atm_call", "ATM call"), ("otm_digital_k130", "OTM digital, K = 130")]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.4), sharey=True)
    grid = [-4 + 0.05 * i for i in range(161)]
    bins = [-4.5 + 0.5 * i for i in range(19)]
    for ax, (key, title) in zip(axes, panels):
        z = [float(r["z_score"]) for r in rows if r["payoff"] == key]
        covered = sum(int(r["covered_95"]) for r in rows if r["payoff"] == key)
        ax.hist(z, bins=bins, density=True, color=SERIES[0], edgecolor=SURFACE, linewidth=2)
        ax.plot(grid, [math.exp(-x * x / 2) / math.sqrt(2 * math.pi) for x in grid], color=INK_SECONDARY,
                linewidth=1.5, linestyle=(0, (4, 3)), label="N(0, 1)")
        ax.set_title(f"{title}: {covered / len(z):.1%} of 95 % CIs cover", loc="left", fontsize=10)
        ax.set_xlabel("z = (estimate − exact) / SE")
        ax.set_xlim(-4.5, 4.5)
    axes[0].set_ylabel("density")
    axes[0].legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(out, metadata={"Date": None})
    plt.close(fig)


def log_slope(xs, ys):
    """Least-squares slope of log(ys) against log(xs)."""
    import math

    lx, ly = [math.log(x) for x in xs], [math.log(y) for y in ys]
    mx, my = sum(lx) / len(lx), sum(ly) / len(ly)
    return sum((a - mx) * (b - my) for a, b in zip(lx, ly)) / sum((a - mx) ** 2 for a in lx)


def qmc_convergence(path, meta, out):
    rows = read_rows(path)
    panels = [("atm_call", "ATM call (1-D, kink)"), ("digital_k130", "Digital, K = 130 (1-D, jump)"),
              ("asian_12", "Arithmetic Asian (12-D, kink)"),
              ("asian_digital_12", "Digital on the average (12-D, jump)")]
    methods = [("pseudo_random", "pseudo-random"), ("rqmc", "RQMC"), ("rqmc_bridge", "RQMC + bridge")]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 6.0))
    for ax, (problem, title) in zip(axes.flat, panels):
        for (method, label), color in zip(methods, SERIES + ["#1baf7a"]):
            pts = [(float(r["paths"]), float(r["sd_single_estimate"])) for r in rows
                   if r["problem"] == problem and r["method"] == method]
            if not pts:
                continue
            n, sd = zip(*pts)
            slope = log_slope(n, sd)
            ax.loglog(n, sd, color=color, linewidth=2, marker="o", markersize=3.5, label=label)
            ax.annotate(f"{slope:+.2f}", (n[-1], sd[-1]), xytext=(4, 0), textcoords="offset points",
                        va="center", color=INK_SECONDARY, fontsize=9)
        ax.set_title(title, loc="left", fontsize=10)
        ax.set_xlim(right=2 ** 16 * 2.2)
    for ax in axes[1]:
        ax.set_xlabel("points per estimate N")
    for ax in axes[:, 0]:
        ax.set_ylabel("sd of one estimate")
    fig.legend(*axes[1, 0].get_legend_handles_labels(), loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.0))
    fig.suptitle("Labels: fitted slope of the error in N (pseudo-random: −0.5)", y=0.035, fontsize=9,
                 color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(out, metadata={"Date": None})
    plt.close(fig)


# Each renderer takes (csv path, metadata dict, output path). Experiments without an entry (tables)
# are reported directly from their CSV.
FIGURES = {
    "fd_vcurve": lambda path, meta, out: fd_vcurve(read_csv(path), meta, out),
    "iv_roundtrip": iv_roundtrip,
    "mc_convergence": mc_convergence,
    "mc_coverage": mc_coverage,
    "qmc_convergence": qmc_convergence,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, default=Path("data/results"))
    parser.add_argument("--out", type=Path, default=Path("docs/figures"))
    parser.add_argument("ids", nargs="*", help=f"figures to render (default: all of {sorted(FIGURES)})")
    args = parser.parse_args()

    unknown = set(args.ids) - FIGURES.keys()
    if unknown:
        sys.exit(f"unknown figure(s): {sorted(unknown)}")
    args.out.mkdir(parents=True, exist_ok=True)
    for fig_id in args.ids or sorted(FIGURES):
        meta = json.loads((args.results / f"{fig_id}.meta.json").read_text())
        out = args.out / f"{fig_id}.svg"
        FIGURES[fig_id](args.results / f"{fig_id}.csv", meta, out)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
