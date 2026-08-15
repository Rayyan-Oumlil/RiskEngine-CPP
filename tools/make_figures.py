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


# Each renderer takes (csv path, metadata dict, output path).
FIGURES = {
    "fd_vcurve": lambda path, meta, out: fd_vcurve(read_csv(path), meta, out),
    "iv_roundtrip": iv_roundtrip,
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
