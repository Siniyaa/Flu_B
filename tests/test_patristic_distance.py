"""Tests for patristic distances using a synthetic Auspice tree."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "patristic_distance.py"
EXAMPLE = ROOT / "examples" / "patristic_example.json"
EXPECTED = ROOT / "examples" / "patristic_expected.csv"

spec = importlib.util.spec_from_file_location("patristic_distance", SCRIPT)
assert spec is not None and spec.loader is not None
pd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pd)


class PatristicDistanceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tree = json.loads(EXAMPLE.read_text())["tree"]
        self.nodes, self.parents = pd.build_index(self.tree)
        self.reference = pd.find_reference(self.nodes, "example_reference")
        self.assertIsNotNone(self.reference)

    def run_cli(self, *extra: str, reference: str = "example_reference"):
        """Run the command and return the process result and parsed CSV rows."""
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "distances.csv"
            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(EXAMPLE),
                 "--reference", reference, "--out", str(output), *extra],
                capture_output=True, text=True, check=False,
            )
            rows = []
            if output.exists():
                with output.open(newline="", encoding="utf-8") as handle:
                    rows = list(csv.DictReader(handle))
            return result, rows

    def test_shared_ancestor_path_is_subtracted(self) -> None:
        tip = pd.find_reference(self.nodes, "example_uae")
        lca = pd.find_lca(tip, self.parents,
                          pd.ancestor_set(self.reference, self.parents))
        self.assertEqual(pd.name(lca), "example_shared_ancestor")
        distance = pd.div(tip) + pd.div(self.reference) - 2 * pd.div(lca)
        self.assertAlmostEqual(distance, 0.005)
        self.assertNotAlmostEqual(distance, 0.013)

    def test_reference_distance_is_zero(self) -> None:
        result, rows = self.run_cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        reference_row = next(row for row in rows if row["name"] == "example_reference")
        self.assertEqual(float(reference_row["patristic_distance"]), 0.0)

    def test_expected_csv(self) -> None:
        result, rows = self.run_cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        with EXPECTED.open(newline="", encoding="utf-8") as handle:
            expected = list(csv.DictReader(handle))
        self.assertEqual(rows, expected)

    def test_internal_nodes_are_not_output(self) -> None:
        result, rows = self.run_cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        names = {row["name"] for row in rows}
        self.assertNotIn("example_root", names)
        self.assertNotIn("example_shared_ancestor", names)
        self.assertEqual(len(rows), 4)

    def test_country_aliases_are_selected_explicitly(self) -> None:
        result, rows = self.run_cli("--group-countries", "UAE", "United Arab Emirates")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual({row["name"] for row in rows}, {"example_uae", "example_uae_alias"})
        self.assertEqual(len(rows), 2)

    def test_country_matching_is_case_sensitive(self) -> None:
        result, rows = self.run_cli("--group-countries", "uae")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows, [])

    def test_tip_without_divergence_is_omitted(self) -> None:
        result, rows = self.run_cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("example_missing_div", {row["name"] for row in rows})

    def test_category_boundaries(self) -> None:
        values = [
            (0.0, "Very close"), (0.004999, "Very close"),
            (0.005, "Close"), (0.006999, "Close"),
            (0.007, "Moderate"), (0.009999, "Moderate"),
            (0.010, "Distant"), (0.019999, "Distant"),
            (0.020, "Very distant"), (0.1, "Very distant"),
        ]
        for distance, category in values:
            with self.subTest(distance=distance):
                self.assertEqual(pd.categorize(distance), category)

    def test_missing_reference_exits_with_error(self) -> None:
        result, rows = self.run_cli(reference="absent_reference")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not found in tree", result.stderr)
        self.assertEqual(rows, [])

    def test_missing_reference_divergence_exits_with_error(self) -> None:
        result, rows = self.run_cli(reference="example_missing_div")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("has no node_attrs.div", result.stderr)
        self.assertEqual(rows, [])

    def test_help_command(self) -> None:
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--group-countries", result.stdout)
        self.assertIn("--reference", result.stdout)


if __name__ == "__main__":
    unittest.main()
