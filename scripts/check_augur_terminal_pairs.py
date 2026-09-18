#!/usr/bin/env python3
"""
Check the unit of analysis when matching HA and NA mutation profiles.

Augur mutation JSON files contain phylogenetic records, including reconstructed
internal nodes and terminal virus tips. This script counts matched HA-NA records
at the all-node level and at the terminal-tip level so that burden summaries can
be reported with the correct unit of analysis.
"""

from __future__ import annotations

import argparse
import json
import re
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd


def normalize_name(name: str) -> str:
    return re.sub(r"_(HA|NA)$", "", name)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_json_from_zip(zip_path: Path, inner_path: str) -> dict[str, Any]:
    with zipfile.ZipFile(zip_path, "r") as archive:
        with archive.open(inner_path) as handle:
            return json.load(handle)


def terminal_names(tree: dict[str, Any]) -> set[str]:
    tips: set[str] = set()
    stack = [tree]
    while stack:
        node = stack.pop()
        children = node.get("children", [])
        if children:
            stack.extend(children)
        else:
            tips.add(normalize_name(node.get("name", "")))
    return tips


def node_names_from_aa_muts(aa_muts: dict[str, Any]) -> set[str]:
    nodes = aa_muts.get("nodes", {})
    return {normalize_name(str(name)) for name in nodes.keys()}


def main() -> None:
    parser = argparse.ArgumentParser(description="Count matched HA-NA all-node and terminal-tip records.")
    parser.add_argument("--ha-aa-muts", type=Path, required=True, help="HA aa_muts JSON file.")
    parser.add_argument("--ha-auspice", type=Path, required=True, help="HA Auspice JSON file.")
    parser.add_argument("--na-aa-muts", type=Path, help="NA aa_muts JSON file.")
    parser.add_argument("--na-auspice", type=Path, help="NA Auspice JSON file.")
    parser.add_argument("--na-zip", type=Path, help="ZIP file containing NA JSON files.")
    parser.add_argument("--na-aa-muts-inner", help="Path to NA aa_muts JSON inside ZIP.")
    parser.add_argument("--na-auspice-inner", help="Path to NA Auspice JSON inside ZIP.")
    parser.add_argument("--outdir", type=Path, default=Path("results/unit_of_analysis"))
    args = parser.parse_args()

    if args.na_zip:
        if not args.na_aa_muts_inner or not args.na_auspice_inner:
            raise SystemExit("When --na-zip is used, provide --na-aa-muts-inner and --na-auspice-inner.")
        na_aa = load_json_from_zip(args.na_zip, args.na_aa_muts_inner)
        na_auspice = load_json_from_zip(args.na_zip, args.na_auspice_inner)
    else:
        if not args.na_aa_muts or not args.na_auspice:
            raise SystemExit("Provide either --na-zip with inner paths, or --na-aa-muts and --na-auspice.")
        na_aa = load_json(args.na_aa_muts)
        na_auspice = load_json(args.na_auspice)

    ha_aa = load_json(args.ha_aa_muts)
    ha_auspice = load_json(args.ha_auspice)

    ha_nodes = node_names_from_aa_muts(ha_aa)
    na_nodes = node_names_from_aa_muts(na_aa)
    ha_tips = terminal_names(ha_auspice["tree"])
    na_tips = terminal_names(na_auspice["tree"])

    all_node_matches = ha_nodes & na_nodes
    internal_node_matches = {x for x in all_node_matches if x.startswith("NODE_")}
    non_node_matches = all_node_matches - internal_node_matches
    terminal_tip_matches = ha_tips & na_tips

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame([
        {"metric": "HA aa_muts node records", "n": len(ha_nodes)},
        {"metric": "NA aa_muts node records", "n": len(na_nodes)},
        {"metric": "Matched HA-NA aa_muts records", "n": len(all_node_matches)},
        {"metric": "Matched internal NODE records", "n": len(internal_node_matches)},
        {"metric": "Matched non-NODE records", "n": len(non_node_matches)},
        {"metric": "Matched terminal Auspice tips", "n": len(terminal_tip_matches)},
    ])
    summary.to_csv(args.outdir / "ha_na_unit_of_analysis_summary.csv", index=False)

    pd.Series(sorted(terminal_tip_matches), name="strain").to_csv(args.outdir / "matched_terminal_tips.csv", index=False)
    pd.Series(sorted(internal_node_matches), name="node").to_csv(args.outdir / "matched_internal_nodes.csv", index=False)

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
