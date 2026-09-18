"""Input validation checks using small artificial sequences, not study data."""
from pathlib import Path
import csv
import json
import sys
import tempfile
import unittest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from prepare_phylogeny_inputs import prepare, fasta_records, read_metadata


class InputValidationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)
        self.fasta = self.path / 'sequences.fasta'
        self.fasta.write_text('>ref\nATGAAATAA\n>tip1\nATGAGATAA\n>tip2\nATGACATAA\n>tip3\nATGANNTAA\n')
        self.meta = self.path / 'metadata.tsv'
        self.meta.write_text('strain\tdate\tcountry\tregion\tclade\tage\n'
                             'ref\t2021-01-09\tAustria\tEurope\tA\t30\n'
                             'tip1\t2022-01-01\tUAE\tAsia\tA\t8\n'
                             'tip2\t2022-XX-XX\t\t\tB\t11\n'
                             'tip3\t2023-03-01\tUAE\tAsia\tB\t15\n')
        self.gb = self.path / 'reference.gb'
        rec = SeqRecord(Seq('ATGAAATAA'), id='ref', name='ref', description='Software test reference')
        rec.annotations['molecule_type'] = 'DNA'
        rec.features = [SeqFeature(FeatureLocation(0, 9), type='CDS', qualifiers={'gene':['HA']})]
        SeqIO.write(rec, self.gb, 'genbank')
        self.out = self.path / 'out'

    def run_prepare(self, **kwargs):
        return prepare(self.fasta, self.meta, self.gb, 'ref', 'HA', self.out, **kwargs)

    def test_valid_inputs(self):
        result = self.run_prepare()
        self.assertEqual(result['retained_sequences'], 4)
        self.assertEqual((self.out / 'include_reference.txt').read_text(), 'ref\n')

    def test_duplicate_fasta_identifier_rejected(self):
        self.fasta.write_text(self.fasta.read_text() + '>tip1\nATGAAATAA\n')
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            self.run_prepare()

    def test_duplicate_metadata_identifier_rejected(self):
        self.meta.write_text(self.meta.read_text() + 'tip1\t2022-01-01\tUAE\tAsia\tA\t8\n')
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            self.run_prepare()

    def test_whitespace_header_rejected(self):
        self.fasta.write_text(self.fasta.read_text().replace('>tip1', '>tip1 extra'))
        with self.assertRaisesRegex(ValueError, 'single identifiers'):
            self.run_prepare()

    def test_internal_style_identifier_rejected(self):
        self.fasta.write_text(self.fasta.read_text().replace('>tip1', '>NODE_001'))
        with self.assertRaisesRegex(ValueError, 'invalid terminal'):
            self.run_prepare()

    def test_reference_sequence_mismatch_rejected(self):
        self.fasta.write_text(self.fasta.read_text().replace('>ref\nATGAAATAA', '>ref\nATGACATAA'))
        with self.assertRaisesRegex(ValueError, 'differ'):
            self.run_prepare()

    def test_missing_metadata_rejected(self):
        self.meta.write_text('\n'.join(line for line in self.meta.read_text().splitlines() if not line.startswith('tip1\t'))+'\n')
        with self.assertRaisesRegex(ValueError, 'lack metadata'):
            self.run_prepare()

    def test_unknown_calendar_date_not_invented(self):
        self.run_prepare()
        fields, rows = read_metadata(self.out / 'metadata_observed.tsv')
        self.assertEqual(rows['tip2']['date'], '2022-XX-XX')

    def test_invalid_complete_date_rejected(self):
        self.meta.write_text(self.meta.read_text().replace('2022-01-01', '2022-02-31'))
        with self.assertRaises(ValueError):
            self.run_prepare()

    def test_ambiguity_threshold_recorded(self):
        result = self.run_prepare(max_ambiguous_fraction=0.1)
        self.assertEqual(result['retained_sequences'], 3)
        records = fasta_records(self.out / 'sequences.fasta')
        self.assertNotIn('tip3', records)
        self.assertIn('ambiguity_above_threshold', (self.out / 'sequence_qc.tsv').read_text())

    def test_browser_export_excludes_clinical_fields(self):
        self.run_prepare()
        fields, rows = read_metadata(self.out / 'metadata_export.tsv')
        self.assertNotIn('age', fields)
        self.assertIn('clade', fields)

    def test_observed_metadata_distinct_from_inference_copy(self):
        self.run_prepare()
        _, observed = read_metadata(self.out / 'metadata_observed.tsv')
        _, analysis = read_metadata(self.out / 'metadata.tsv')
        self.assertEqual(observed['tip2']['country'], '')
        self.assertEqual(analysis['tip2']['country'], '?')

    def test_identity_does_not_merge_distinct_tips(self):
        self.fasta.write_text(self.fasta.read_text().replace('ATGAGATAA', 'ATGAAATAA'))
        result = self.run_prepare()
        self.assertEqual(result['retained_sequences'], 4)


if __name__ == '__main__':
    unittest.main()
