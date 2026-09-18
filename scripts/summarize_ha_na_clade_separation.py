#!/usr/bin/env python3
"""
Summarize clade separation and the unadjusted HA-burden/NA-substitution association.

The script uses aggregate count tables supplied in the data directory. It reports:
  1. high HA-epitope burden by clade,
  2. the unadjusted Fisher exact test for focal NA-substitution carriage by
     HA-burden class,
  3. a short interpretation note for the binary high-burden contrast.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


def woolf_ci(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    odds = (a * d) / (b * c)
    se = math.sqrt((1 / a) + (1 / b) + (1 / c) + (1 / d))
    z = 1.959963984540054
    return math.exp(math.log(odds) - z * se), math.exp(math.log(odds) + z * se)


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize HA-NA clade separation and crude Fisher association.")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--outdir", type=Path, default=Path("results"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    high_by_clade = pd.read_csv(args.data_dir / "figure6_high_burden_by_clade.csv")
    burden = pd.read_csv(args.data_dir / "figure6_na_by_burden.csv")

    high_by_clade["percent_high_burden"] = 100 * high_by_clade["high_burden"] / high_by_clade["total"]
    high_by_clade.to_csv(args.outdir / "high_burden_by_clade_summary.csv", index=False)

    low = burden.loc[burden["burden_class"].str.startswith("Low")].iloc[0]
    high = burden.loc[burden["burden_class"].str.startswith("High")].iloc[0]

    table = [[int(high["at_least_one_focal_NA"]), int(high["no_focal_NA"])],
             [int(low["at_least_one_focal_NA"]), int(low["no_focal_NA"])]]
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    ci_low, ci_high = woolf_ci(table[0][0], table[0][1], table[1][0], table[1][1])

    fisher_summary = pd.DataFrame([
        {
            "comparison": "high_vs_low_HA_epitope_burden",
            "high_burden_focal_NA_yes": table[0][0],
            "high_burden_focal_NA_no": table[0][1],
            "low_burden_focal_NA_yes": table[1][0],
            "low_burden_focal_NA_no": table[1][1],
            "odds_ratio": odds_ratio,
            "ci_low_woolf": ci_low,
            "ci_high_woolf": ci_high,
            "fisher_exact_p": p_value,
            "interpretation": "descriptive_unadjusted_association_clade_separated",
        }
    ])
    fisher_summary.to_csv(args.outdir / "ha_burden_na_fisher_summary.csv", index=False)

    notes = [
        "HA-NA clade-separation summary",
        "================================",
        "",
        "High HA-epitope burden is separated by clade in the UAE matched HA-NA cohort.",
        "The crude high-burden versus low-burden comparison is therefore interpreted as descriptive.",
        "Adjusted inference for this binary contrast should use exact or stratified analyses;",
        "Firth-penalized logistic regression can be used as a sensitivity method where estimable.",
        "",
        f"Unadjusted Fisher OR = {odds_ratio:.2f}",
        f"95% CI = {ci_low:.2f}-{ci_high:.2f}",
        f"Fisher exact p = {p_value:.4g}",
    ]
    (args.outdir / "ha_na_clade_separation_summary.txt").write_text("\n".join(notes) + "\n")


if __name__ == "__main__":
    main()
