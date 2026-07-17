#!/usr/bin/env python3
"""
patristic_distance.py

Compute patristic (tree-path) distance from every leaf in an Auspice v2
phylogeny (Nextstrain / Augur ``*.json``) to a chosen reference tip, and bin
the results into divergence categories.

Patristic distance between two tips A and B on a tree is the sum of branch
lengths along the unique path connecting them. Using the cumulative
root-to-tip divergence (``node_attrs.div``, in substitutions/site) that Auspice
stores on every node, this reduces to:

    d(A, B) = div(A) + div(B) - 2 * div(LCA(A, B))

where LCA is the lowest common ancestor of A and B. The ``- 2 * div(LCA)`` term
removes the shared ``root -> LCA`` path that would otherwise be double-counted.

This is the standard tree-path distance used in Nextstrain, TreeTime, BEAST,
ete3 and dendropy. It is NOT raw Hamming/alignment distance; it inherits the
branch lengths already estimated in the input tree.

Example
-------
    python patristic_distance.py flu.json \
        --reference EPI_ISL_983345 \
        --group-countries "UAE" "United Arab Emirates" \
        --out uae_patristic_distances.csv

Author: (add your name)
License: MIT (or your choice)
"""

import argparse
import csv
import json
import sys
from statistics import mean, median


# ---------------------------------------------------------------------------
# Auspice node accessors
# ---------------------------------------------------------------------------
def name(node):
    return node.get("name", "")


def div(node):
    """Cumulative root-to-tip divergence (substitutions/site)."""
    return node.get("node_attrs", {}).get("div")


def country(node):
    return node.get("node_attrs", {}).get("country", {}).get("value")


def date_str(node):
    return node.get("node_attrs", {}).get("date", {}).get("value")


def is_leaf(node):
    return not node.get("children")


# ---------------------------------------------------------------------------
# Tree indexing
# ---------------------------------------------------------------------------
def build_index(root):
    """Single pass over the tree.

    Returns
    -------
    all_nodes : dict[int, node]   id(node) -> node
    parent    : dict[int, node]   id(node) -> parent node (None for root)
    """
    all_nodes = {}
    parent = {}

    stack = [(root, None)]
    while stack:
        node, par = stack.pop()
        all_nodes[id(node)] = node
        parent[id(node)] = par
        for child in node.get("children", []):
            stack.append((child, node))

    return all_nodes, parent


def find_reference(all_nodes, ref_id):
    for node in all_nodes.values():
        if name(node) == ref_id:
            return node
    return None


def ancestor_set(node, parent):
    """Set of id()s for every node on the path from `node` up to the root."""
    ancestors = set()
    cur = node
    while cur is not None:
        ancestors.add(id(cur))
        cur = parent[id(cur)]
    return ancestors


def find_lca(node, parent, ref_ancestors):
    """First ancestor of `node` that is also an ancestor of the reference."""
    cur = node
    while cur is not None:
        if id(cur) in ref_ancestors:
            return cur
        cur = parent[id(cur)]
    return None


# ---------------------------------------------------------------------------
# Distance categories (final manuscript thresholds, subs/site)
# ---------------------------------------------------------------------------
CATEGORIES = [
    ("Very close",   lambda d: d < 0.005),
    ("Close",        lambda d: 0.005 <= d < 0.007),
    ("Moderate",     lambda d: 0.007 <= d < 0.010),
    ("Distant",      lambda d: 0.010 <= d < 0.020),
    ("Very distant", lambda d: d >= 0.020),
]


def categorize(d):
    for label, test in CATEGORIES:
        if test(d):
            return label
    return "Unclassified"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Patristic distance from tree leaves to a reference tip "
                    "in an Auspice v2 JSON."
    )
    ap.add_argument("json_file", help="Auspice v2 phylogeny JSON (e.g. flu.json)")
    ap.add_argument("--reference", required=True,
                    help="Reference tip name/id (e.g. EPI_ISL_983345)")
    ap.add_argument("--group-countries", nargs="*", default=None,
                    metavar="COUNTRY",
                    help="If given, only leaves whose country matches one of "
                         "these values are reported (e.g. UAE \"United Arab Emirates\"). "
                         "Omit to score every leaf in the tree.")
    ap.add_argument("--out", default="patristic_distances.csv",
                    help="Output CSV path (default: patristic_distances.csv)")
    args = ap.parse_args(argv)

    with open(args.json_file) as fh:
        data = json.load(fh)

    tree = data.get("tree")
    if tree is None:
        sys.exit("Error: no 'tree' key found — is this an Auspice v2 JSON?")

    all_nodes, parent = build_index(tree)

    reference = find_reference(all_nodes, args.reference)
    if reference is None:
        sys.exit(f"Error: reference '{args.reference}' not found in tree.")
    if div(reference) is None:
        sys.exit(f"Error: reference '{args.reference}' has no node_attrs.div.")

    ref_div = div(reference)
    ref_ancestors = ancestor_set(reference, parent)

    print(f"Reference : {name(reference)}")
    print(f"  country : {country(reference)}")
    print(f"  date    : {date_str(reference)}")
    print(f"  div     : {ref_div:.5f} subs/site\n")

    # Select leaves
    group = set(args.group_countries) if args.group_countries else None
    leaves = [
        n for n in all_nodes.values()
        if is_leaf(n) and div(n) is not None
        and (group is None or country(n) in group)
    ]
    print(f"Scoring {len(leaves)} leaves against the reference...\n")

    rows = []
    for leaf in leaves:
        lca = find_lca(leaf, parent, ref_ancestors)
        if lca is None:
            continue
        d = div(leaf) + ref_div - 2 * div(lca)
        # Numerical guard: identical/degenerate paths can yield tiny negatives
        if d < 0:
            d = 0.0
        rows.append({
            "name": name(leaf),
            "country": country(leaf),
            "date": date_str(leaf),
            "div": round(div(leaf), 6),
            "patristic_distance": round(d, 6),
            "category": categorize(d),
        })

    rows.sort(key=lambda r: r["patristic_distance"])

    # Write CSV
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["name", "country", "date", "div",
                        "patristic_distance", "category"],
        )
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    dists = [r["patristic_distance"] for r in rows]
    if dists:
        print("Distance summary (subs/site):")
        print(f"  n      : {len(dists)}")
        print(f"  min    : {min(dists):.5f}")
        print(f"  median : {median(dists):.5f}")
        print(f"  mean   : {mean(dists):.5f}")
        print(f"  max    : {max(dists):.5f}\n")

        print("Category distribution:")
        n = len(dists)
        for label, _ in CATEGORIES:
            c = sum(1 for r in rows if r["category"] == label)
            print(f"  {label:<13}: {c:>4}  ({100 * c / n:5.1f}%)")

    print(f"\nWrote {len(rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()