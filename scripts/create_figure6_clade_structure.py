#!/usr/bin/env python3
"""
Create Figure 6 for the influenza B/Victoria HA-NA co-occurrence analysis.

The figure summarizes clade composition, clade separation of high HA-epitope
burden, and the unadjusted association between HA-epitope burden and carriage
of at least one focal NA substitution.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import fisher_exact


BLUE = "#3B74B4"
GOLD = "#DCA135"
GREEN = "#579E7A"
GREY = "#BDBDBD"
EDGE = "#333333"
TEXT = "#222222"


def odds_ratio_ci_woolf(a: int, b: int, c: int, d: int, alpha: float = 0.05) -> tuple[float, float]:
    """
    Woolf log-odds confidence interval for a 2 x 2 table.

    Table convention:
        [[a, b],
         [c, d]]
    where a = high burden with focal NA substitution, b = high burden without,
    c = low burden with focal NA substitution, d = low burden without.
    """
    if min(a, b, c, d) == 0:
        # Haldane-Anscombe correction only if required.
        a, b, c, d = (a + 0.5, b + 0.5, c + 0.5, d + 0.5)

    odds = (a * d) / (b * c)
    se = math.sqrt((1 / a) + (1 / b) + (1 / c) + (1 / d))
    z = 1.959963984540054
    lo = math.exp(math.log(odds) - z * se)
    hi = math.exp(math.log(odds) + z * se)
    return lo, hi


def read_inputs(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    clade = pd.read_csv(data_dir / "figure6_clade_counts.csv")
    high = pd.read_csv(data_dir / "figure6_high_burden_by_clade.csv")
    burden = pd.read_csv(data_dir / "figure6_na_by_burden.csv")
    return clade, high, burden


def save_figure(fig: plt.Figure, outdir: Path, basename: str, formats: Iterable[str], dpi: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fmt = fmt.lower().lstrip(".")
        path = outdir / f"{basename}.{fmt}"
        if fmt in {"png", "tif", "tiff"}:
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
        else:
            fig.savefig(path, bbox_inches="tight")


def make_figure(clade: pd.DataFrame, high: pd.DataFrame, burden: pd.DataFrame) -> tuple[plt.Figure, dict[str, float]]:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#555555",
        "xtick.color": "#555555",
        "ytick.color": "#555555",
        "axes.labelcolor": TEXT,
        "text.color": TEXT,
    })

    fig = plt.figure(figsize=(10.5, 6.4), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.45], hspace=0.45, wspace=0.55)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    clade_colors = [BLUE if c == "V1A.3a.2" else GOLD if c == "V1A" else "#D9D9D9" for c in clade["clade"]]

    # Panel A
    ax_a.bar(clade["clade"], clade["n"], color=clade_colors, edgecolor=EDGE, linewidth=0.8)
    ax_a.set_title("Clade composition of paired HA-NA isolates", pad=12)
    ax_a.set_ylabel("Matched isolates (n)")
    ax_a.set_ylim(0, 430)
    for i, row in clade.iterrows():
        ax_a.text(i, row["n"] + 10, f"{int(row['n'])}\n({row['percent']:.1f}%)", ha="center", va="bottom", color="#555555", linespacing=0.9)

    # Panel B
    high["percent"] = high["high_burden"] / high["total"] * 100
    ax_b.bar(high["clade"], high["percent"], color=clade_colors, edgecolor=EDGE, linewidth=0.8)
    ax_b.set_title("High HA-epitope burden is clade-separated", pad=12)
    ax_b.set_ylabel("High HA-epitope burden (>5), %")
    ax_b.set_ylim(0, 110)
    for i, row in high.iterrows():
        y = max(row["percent"], 1.5)
        ax_b.text(i, y + 2.5, f"{int(row['high_burden'])}/{int(row['total'])}", ha="center", va="bottom", color="#555555")

    # Panel C
    low = burden.loc[burden["burden_class"].str.startswith("Low")].iloc[0]
    high_b = burden.loc[burden["burden_class"].str.startswith("High")].iloc[0]

    classes = ["Low burden\n<=5 (n=375)", "High burden\n>5 (n=82)"]
    no_counts = [low["no_focal_NA"], high_b["no_focal_NA"]]
    yes_counts = [low["at_least_one_focal_NA"], high_b["at_least_one_focal_NA"]]
    totals = [low["total"], high_b["total"]]
    no_pct = [100 * no_counts[i] / totals[i] for i in range(2)]
    yes_pct = [100 * yes_counts[i] / totals[i] for i in range(2)]

    x = range(2)
    width = 0.48
    ax_c.bar(x, no_pct, width=width, color=GREY, edgecolor=EDGE, linewidth=0.8, label="No focal NA substitution")
    ax_c.bar(x, yes_pct, bottom=no_pct, width=width, color=GREEN, edgecolor=EDGE, linewidth=0.8, label="At least one focal NA substitution")
    ax_c.set_title("Unadjusted HA-epitope burden and focal NA-substitution carriage", pad=14)
    ax_c.set_ylabel("Isolates within burden class (%)")
    ax_c.set_xticks(list(x), classes)
    ax_c.set_ylim(0, 100)
    ax_c.legend(loc="upper left", frameon=True, framealpha=1.0, edgecolor="#CCCCCC")

    for i in x:
        ax_c.text(i, no_pct[i] / 2, f"{int(no_counts[i])}", ha="center", va="center", color=TEXT, fontsize=10)
        ax_c.text(i, no_pct[i] + yes_pct[i] / 2, f"{int(yes_counts[i])}", ha="center", va="center", color="white", fontsize=11, fontweight="bold")

    # Fisher exact test with high-burden row first, low-burden row second.
    table = [[int(high_b["at_least_one_focal_NA"]), int(high_b["no_focal_NA"])],
             [int(low["at_least_one_focal_NA"]), int(low["no_focal_NA"])]]
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    ci_low, ci_high = odds_ratio_ci_woolf(table[0][0], table[0][1], table[1][0], table[1][1])

    ax_c.text(
        0.5,
        0.63,
        f"Unadjusted Fisher test\nOR = {odds_ratio:.2f} (95% CI {ci_low:.2f}-{ci_high:.2f})\np = {p_value:.4f}; clade-separated",
        transform=ax_c.transAxes,
        ha="center",
        va="center",
        fontsize=9.5,
        bbox=dict(facecolor="white", edgecolor="#CCCCCC", boxstyle="round,pad=0.35"),
    )

    # Panel labels aligned in figure coordinates to avoid overlap with axis labels.
    fig.text(0.015, 0.965, "A", fontsize=14, fontweight="bold", va="top", ha="left")
    fig.text(0.585, 0.965, "B", fontsize=14, fontweight="bold", va="top", ha="left")
    fig.text(0.015, 0.535, "C", fontsize=14, fontweight="bold", va="top", ha="left")

    for ax in [ax_a, ax_b, ax_c]:
        ax.spines["left"].set_linewidth(0.9)
        ax.spines["bottom"].set_linewidth(0.9)

    fig.subplots_adjust(left=0.07, right=0.98, top=0.93, bottom=0.11)

    stats = {
        "low_total": int(low["total"]),
        "low_focal_NA_yes": int(low["at_least_one_focal_NA"]),
        "low_focal_NA_no": int(low["no_focal_NA"]),
        "high_total": int(high_b["total"]),
        "high_focal_NA_yes": int(high_b["at_least_one_focal_NA"]),
        "high_focal_NA_no": int(high_b["no_focal_NA"]),
        "odds_ratio": odds_ratio,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "fisher_p": p_value,
    }
    return fig, stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create Figure 6 for HA-NA clade-structure analysis.")
    parser.add_argument("--data-dir", type=Path, default=Path("data"), help="Directory containing input CSV files.")
    parser.add_argument("--outdir", type=Path, default=Path("figures"), help="Output directory for figure files.")
    parser.add_argument("--results-dir", type=Path, default=Path("results"), help="Output directory for statistics summary.")
    parser.add_argument("--basename", default="Figure6_clade_structure", help="Base filename for figure outputs.")
    parser.add_argument("--formats", nargs="+", default=["png", "pdf", "svg"], help="Output formats.")
    parser.add_argument("--dpi", type=int, default=600, help="Raster image resolution.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    clade, high, burden = read_inputs(args.data_dir)
    fig, stats = make_figure(clade, high, burden)
    save_figure(fig, args.outdir, args.basename, args.formats, args.dpi)

    args.results_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([stats]).to_csv(args.results_dir / f"{args.basename}_statistics.csv", index=False)
    plt.close(fig)


if __name__ == "__main__":
    main()
