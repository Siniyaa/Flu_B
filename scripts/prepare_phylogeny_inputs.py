#!/usr/bin/env python3
"""Validate and prepare consensus FASTA and metadata for an Augur segment build."""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import json
import re
from pathlib import Path
from typing import Any
from Bio import SeqIO

NUCLEOTIDES = set('ACGTRYSWKMBDHVN-')


def fasta_records(path: Path) -> dict[str, Any]:
    records = {}
    with path.open(encoding='utf-8') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.description != record.id:
                raise ValueError(f'{path}: FASTA headers must be single identifiers without spaces: {record.description!r}')
            if record.id in records:
                raise ValueError(f'{path}: duplicate sequence identifier {record.id!r}; resolve duplicate records before running.')
            if not record.id or record.id.startswith('NODE_'):
                raise ValueError(f'{path}: invalid terminal identifier {record.id!r}.')
            record.seq = record.seq.upper()
            if not len(record):
                raise ValueError(f'{path}: empty sequence {record.id!r}.')
            invalid = set(str(record.seq)) - NUCLEOTIDES
            if invalid:
                raise ValueError(f'{path}: unsupported nucleotide symbols in {record.id}: {sorted(invalid)}')
            records[record.id] = record
    if not records:
        raise ValueError(f'{path}: no FASTA records.')
    return records


def read_metadata(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError(f'{path}: duplicate metadata column names.')
        if not {'strain', 'date'}.issubset(fields):
            raise ValueError(f'{path}: tab-separated metadata requires strain and date columns.')
        rows = {}
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'{path}: malformed tab-separated row.')
            identifier = row['strain'].strip()
            if not identifier or re.search(r'\s', identifier):
                raise ValueError(f'{path}: empty or whitespace-containing strain identifier.')
            if identifier in rows:
                raise ValueError(f'{path}: duplicate metadata identifier {identifier!r}.')
            row['strain'] = identifier
            date = row['date'].strip()
            if not date:
                raise ValueError(f'{path}: missing date for {identifier}.')
            # Augur accepts incomplete dates; record them without inventing precision.
            if not re.fullmatch(r'[0-9X]{4}(?:-[0-9X]{2}(?:-[0-9X]{2})?)?', date):
                raise ValueError(f'{path}: invalid date {date!r} for {identifier}.')
            if re.fullmatch(r'\d{4}-\d{2}-\d{2}', date):
                dt.date.fromisoformat(date)
            row['date'] = date
            rows[identifier] = row
    return fields, rows


def validate_reference(path: Path, records: dict[str, Any], reference: str, gene: str) -> None:
    with path.open(encoding='utf-8') as handle:
        refs = list(SeqIO.parse(handle, 'genbank'))
    if len(refs) != 1:
        raise ValueError('Reference GenBank must contain exactly one record.')
    ref = refs[0]
    if reference not in records:
        raise ValueError(f'Reference tip {reference!r} is absent from the input FASTA.')
    if str(ref.seq).upper() != str(records[reference].seq).upper():
        raise ValueError('Reference GenBank sequence and selected reference-tip FASTA sequence differ. Check coordinate systems and identifiers.')
    features = [f for f in ref.features if f.type == 'CDS' and gene in f.qualifiers.get('gene', [])]
    if len(features) != 1:
        raise ValueError(f'Reference GenBank must contain one CDS with /gene="{gene}".')
    feature = features[0]
    if feature.location is None or int(feature.location.start) < 0 or int(feature.location.end) > len(ref.seq):
        raise ValueError('Reference CDS is outside the supplied sequence.')
    if feature.qualifiers.get('codon_start', ['1'])[0] != '1':
        raise ValueError('The CDS requires a non-default codon_start; verify an Augur-compatible annotation before running.')
    if len(feature.extract(ref.seq)) % 3:
        raise ValueError('Reference CDS length is not divisible by three. Verify the annotation; do not infer coordinates from a mutation label.')


def prepare(sequences: Path, metadata: Path, reference_genbank: Path, reference_tip: str,
            gene: str, outdir: Path, max_ambiguous_fraction: float | None = None,
            require_exact_ids: bool = True, trait_columns: tuple[str, ...] = ('country', 'region')) -> dict:
    records = fasta_records(sequences)
    fields, rows = read_metadata(metadata)
    validate_reference(reference_genbank, records, reference_tip, gene)
    missing_metadata = sorted(set(records) - set(rows))
    metadata_only = sorted(set(rows) - set(records))
    if require_exact_ids and missing_metadata:
        raise ValueError(f'{len(missing_metadata)} sequences lack metadata, including {missing_metadata[:3]}.')
    if reference_tip not in rows:
        raise ValueError('Reference metadata are required; do not invent a collection date.')
    if not set(trait_columns).issubset(fields):
        raise ValueError(f'Metadata are missing trait columns: {sorted(set(trait_columns) - set(fields))}.')
    if max_ambiguous_fraction is not None and not 0 <= max_ambiguous_fraction <= 1:
        raise ValueError('max_ambiguous_fraction must be null or between zero and one.')
    selected, qc_rows = [], []
    for identifier, record in records.items():
        seq = str(record.seq)
        called = sum(seq.count(base) for base in 'ACGT')
        fraction = (len(seq) - called) / len(seq)
        reason = ''
        if identifier not in rows:
            reason = 'missing_metadata'
        elif max_ambiguous_fraction is not None and fraction > max_ambiguous_fraction:
            reason = 'ambiguity_above_threshold'
        if identifier == reference_tip and reason:
            raise ValueError(f'Reference fails input QC: {reason}.')
        qc_rows.append({'strain': identifier, 'length': len(seq), 'called_bases': called,
                        'ambiguous_fraction': fraction, 'status': reason or 'retained'})
        if not reason:
            selected.append(record)
    if len(selected) < 3:
        raise ValueError('Fewer than three eligible sequences remain.')
    outdir.mkdir(parents=True, exist_ok=True)
    SeqIO.write(selected, str(outdir / 'sequences.fasta'), 'fasta')
    with (outdir / 'metadata.tsv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for rec in selected:
            row = dict(rows[rec.id])
            # '?' is Augur's missing-state marker. Only the analysis copy changes.
            for field in trait_columns:
                if not row[field].strip():
                    row[field] = '?'
            writer.writerow(row)
    # Restrict the browser export to sequence, sampling and lineage fields.
    export_fields = ['strain', 'date'] + [f for f in ('country', 'region', 'clade') if f in fields]
    with (outdir / 'metadata_export.tsv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, export_fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows[rec.id] for rec in selected)
    # Keep observed metadata for iTOL, separate from reconstructed trait values.
    with (outdir / 'metadata_observed.tsv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows[rec.id] for rec in selected)
    with (outdir / 'sequence_qc.tsv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, list(qc_rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(qc_rows)
    summary = {'input_sequences': len(records), 'input_metadata_rows': len(rows),
               'retained_sequences': len(selected), 'metadata_without_sequence': metadata_only,
               'sequences_without_metadata': missing_metadata,
               'max_ambiguous_fraction': max_ambiguous_fraction,
               'reference_tip': reference_tip, 'gene': gene}
    (outdir / 'input_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (outdir / 'include_reference.txt').write_text(reference_tip + '\n')
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sequences', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--reference-genbank', type=Path, required=True)
    parser.add_argument('--reference-tip', required=True)
    parser.add_argument('--gene', required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--max-ambiguous-fraction', type=float)
    parser.add_argument('--allow-unmatched-sequences', action='store_true')
    parser.add_argument('--trait-columns', nargs='*', default=['country', 'region'])
    args = parser.parse_args()
    try:
        result = prepare(args.sequences, args.metadata, args.reference_genbank, args.reference_tip,
                         args.gene, args.outdir, args.max_ambiguous_fraction,
                         not args.allow_unmatched_sequences, tuple(args.trait_columns))
    except (ValueError, OSError) as error:
        parser.exit(1, f'Input error: {error}\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
