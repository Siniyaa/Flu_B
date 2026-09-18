"""Newick export and annotation tests using an artificial tree."""
from pathlib import Path
import copy
import csv
import json
import sys
import tempfile
import unittest
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from export_itol import export, validate_tree, newick_name, category_colour


def node(name, div, year, **kwargs):
    return dict(name=name, node_attrs={'div':div,'num_date':{'value':year},
                                     'country':{'value':'UAE'},'clade':{'value':'A'}}, **kwargs)


class ItolExportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)
        self.root = node('ROOT',0.001,2020.0,children=[
            node('INTERNAL',0.004,2021.0,children=[node('ref',0.006,2022.0),node('tip:A',0.007,2023.0)]),
            node('tipB',0.008,2024.0)])
        self.json = self.path / 'tree.json'
        self.out = self.path / 'itol'

    def run_export(self, **kwargs):
        self.json.write_text(json.dumps({'version':'v2','tree':self.root}))
        return export(self.json,self.out,'HA',['country','clade'],**kwargs)

    def test_divergence_tree_preserves_path_distance(self):
        self.run_export()
        tree=Phylo.read(self.out/'HA_divergence.nwk','newick')
        self.assertAlmostEqual(tree.distance('ref','tip:A'),0.005)
        self.assertAlmostEqual(tree.distance('ref','tipB'),0.012)

    def test_time_tree_has_year_durations(self):
        summary=self.run_export()
        tree=Phylo.read(self.out/'HA_time_years.nwk','newick')
        self.assertAlmostEqual(tree.distance('ref','tip:A'),3.0)
        self.assertTrue(summary['time_tree_written'])

    def test_only_terminal_annotations(self):
        self.run_export()
        with (self.out/'HA_tip_metadata.tsv').open() as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual({r['strain'] for r in rows},{'ref','tip:A','tipB'})
        self.assertNotIn('INTERNAL\t', (self.out/'HA_country.txt').read_text())

    def test_observed_country_overrides_inferred_country(self):
        meta=self.path/'metadata.tsv'
        meta.write_text('strain\tcountry\tclade\nref\tAustria\tA\ntip:A\t\tA\ntipB\tUAE\tB\n')
        self.run_export(metadata=meta)
        text=(self.out/'HA_country.txt').read_text()
        self.assertIn('Austria',text)
        self.assertNotIn('tip:A\t',text)

    def test_duplicate_node_names_rejected(self):
        self.root['children'][1]['name']='ref'
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            self.run_export()

    def test_invalid_divergence_rejected(self):
        self.root['children'][0]['node_attrs']['div']=0.0
        with self.assertRaisesRegex(ValueError,'decreases'):
            self.run_export()

    def test_nonfinite_divergence_rejected(self):
        self.root['node_attrs']['div']=float('nan')
        with self.assertRaisesRegex(ValueError,'nonfinite'):
            self.run_export()

    def test_nonmonotonic_dates_skip_time_export_without_changing_source(self):
        self.root['children'][0]['node_attrs']['num_date']['value']=2019.0
        before=copy.deepcopy(self.root)
        summary=self.run_export()
        self.assertTrue((self.out/'HA_divergence.nwk').exists())
        self.assertFalse((self.out/'HA_time_years.nwk').exists())
        self.assertIn('decreases',summary['time_export_warning'])
        self.assertEqual(self.root,before)

    def test_missing_dates_skip_time_export(self):
        del self.root['node_attrs']['num_date']
        summary=self.run_export()
        self.assertFalse(summary['time_tree_written'])
        self.assertIsNotNone(summary['time_export_warning'])

    def test_stale_time_file_removed(self):
        self.run_export()
        self.root['children'][0]['node_attrs']['num_date']['value']=2019.0
        self.run_export()
        self.assertFalse((self.out/'HA_time_years.nwk').exists())

    def test_zero_length_branch_preserved(self):
        self.root['children'][0]['children'][0]['node_attrs']['div']=0.004
        self.run_export()
        tree=Phylo.read(self.out/'HA_divergence.nwk','newick')
        self.assertEqual(next(n for n in tree.get_terminals() if n.name=='ref').branch_length,0.0)

    def test_palette_is_deterministic(self):
        self.assertEqual(category_colour('country:UAE'),category_colour('country:UAE'))
        self.assertRegex(category_colour('country:UAE'),r'^#[a-f0-9]{6}$')

    def test_explicit_no_time_export(self):
        summary=self.run_export(time_tree=False)
        self.assertFalse(summary['time_tree_requested'])
        self.assertIsNone(summary['time_export_warning'])

    def test_newick_reserved_characters_are_quoted(self):
        self.assertEqual(newick_name('tip:A'),"'tip:A'")
        self.assertEqual(newick_name('plain_tip'),'plain_tip')


if __name__=='__main__':
    unittest.main()
