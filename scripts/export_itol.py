#!/usr/bin/env python3
"""Export a Nextstrain Auspice tree as Newick and iTOL annotation datasets."""
from __future__ import annotations
import argparse
import colorsys
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any


def attr(node: dict, key: str) -> Any:
    value = node.get('node_attrs', {}).get(key)
    return value.get('value') if isinstance(value, dict) else value


def walk(root: dict):
    stack = [(root, None)]
    while stack:
        node, parent = stack.pop()
        yield node, parent
        stack.extend((child, node) for child in reversed(node.get('children', [])))


def validate_tree(root: dict, key: str) -> int:
    seen, tips = set(), 0
    for node, parent in walk(root):
        identifier = node.get('name')
        if not isinstance(identifier, str) or not identifier or any(c in identifier for c in '\t\n\r'):
            raise ValueError('Every tree node needs a nonempty identifier without tabs or line breaks.')
        if identifier in seen:
            raise ValueError(f'Duplicate tree identifier: {identifier}')
        seen.add(identifier)
        value = attr(node, key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f'Missing or nonfinite {key!r} at {identifier}.')
        if parent is not None and value < attr(parent, key) - 1e-8:
            raise ValueError(f'{key!r} decreases along the branch to {identifier}.')
        if not node.get('children'):
            tips += 1
    return tips


def newick_name(identifier: str) -> str:
    if re.fullmatch(r'[A-Za-z0-9_.|/+-]+', identifier):
        return identifier
    return "'" + identifier.replace("'", "''") + "'"


def write_newick(root: dict, path: Path, key: str) -> None:
    """Iterative writer; no recursion limit for large phylogenies."""
    validate_tree(root, key)
    stack = [('node', root, None)]
    with path.open('w', encoding='utf-8') as handle:
        while stack:
            kind, item, parent = stack.pop()
            if kind == 'text':
                handle.write(item)
                continue
            length = 0.0 if parent is None else max(0.0, attr(item, key) - attr(parent, key))
            ending = f'{newick_name(item["name"])}:{length:.12g}'
            children = item.get('children', [])
            if not children:
                handle.write(ending)
                continue
            handle.write('(')
            stack.append(('text', ')' + ending, None))
            for index in range(len(children) - 1, -1, -1):
                stack.append(('node', children[index], item))
                if index:
                    stack.append(('text', ',', None))
        handle.write(';\n')


def category_colour(label: str) -> str:
    hue = int.from_bytes(hashlib.sha256(label.encode()).digest()[:4], 'big') / 2**32
    rgb = colorsys.hsv_to_rgb(hue, 0.58, 0.72)
    return '#' + ''.join(f'{round(component * 255):02x}' for component in rgb)


def metadata_rows(path: Path | None) -> dict[str, dict]:
    if path is None:
        return {}
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if 'strain' not in (reader.fieldnames or []):
            raise ValueError('Annotation metadata require a strain column.')
        rows = {}
        for row in reader:
            if row['strain'] in rows:
                raise ValueError('Duplicate metadata identifiers are not permitted.')
            rows[row['strain']] = row
        return rows


def export(auspice: Path, outdir: Path, prefix: str, columns: list[str],
           metadata: Path | None = None, time_tree: bool = True) -> dict:
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', prefix):
        raise ValueError('Prefix must be a filename-safe identifier.')
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+', col) for col in columns):
        raise ValueError('Column names must be filename-safe.')
    data = json.loads(auspice.read_text())
    root = data.get('tree')
    if not isinstance(root, dict):
        raise ValueError('Input must contain an Auspice v2 tree object.')
    n = validate_tree(root, 'div')
    time_requested = time_tree
    time_warning = None
    if time_tree:
        try:
            validate_tree(root, 'num_date')
        except ValueError as error:
            time_tree = False
            time_warning = str(error)
    tips = [node for node, _ in walk(root) if not node.get('children')]
    observed = metadata_rows(metadata)
    outdir.mkdir(parents=True, exist_ok=True)
    # Do not leave a stale temporal file when this export cannot validate dates.
    (outdir / f'{prefix}_time_years.nwk').unlink(missing_ok=True)
    for col in columns:
        (outdir / f'{prefix}_{col}.txt').unlink(missing_ok=True)
    write_newick(root, outdir / f'{prefix}_divergence.nwk', 'div')
    if time_tree:
        write_newick(root, outdir / f'{prefix}_time_years.nwk', 'num_date')
    rows = []
    for tip in tips:
        identifier = tip['name']
        # A supplied metadata file is authoritative for observed tip annotations.
        row = {'strain': identifier}
        for col in columns:
            value = observed.get(identifier, {}).get(col) if metadata is not None else attr(tip, col)
            value = '' if value is None else str(value).strip()
            row[col] = '' if value in {'?', 'NA', 'nan', 'None'} else value
            if any(c in row[col] for c in '\t\r\n'):
                raise ValueError('Annotation values cannot contain tabs or newlines.')
        rows.append(row)
    with (outdir / f'{prefix}_tip_metadata.tsv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, ['strain'] + columns, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    for col in columns:
        labels = sorted({row[col] for row in rows if row[col]})
        if not labels:
            continue
        palette = {label: category_colour(f'{col}:{label}') for label in labels}
        lines = ['DATASET_COLORSTRIP', 'SEPARATOR TAB', f'DATASET_LABEL\t{col}',
                 'COLOR\t#555555', 'COLOR_BRANCHES\t0', f'LEGEND_TITLE\t{col}',
                 'LEGEND_SHAPES\t' + '\t'.join('1' for _ in labels),
                 'LEGEND_COLORS\t' + '\t'.join(palette[label] for label in labels),
                 'LEGEND_LABELS\t' + '\t'.join(labels), 'DATA']
        lines.extend(f'{row["strain"]}\t{palette[row[col]]}\t{row[col]}' for row in rows if row[col])
        (outdir / f'{prefix}_{col}.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    summary = {'terminal_tips': n, 'source': str(auspice),
               'annotations': 'observed metadata' if metadata else 'Auspice node attributes; values may include trait inference',
               'divergence_units': 'same units as source node_attrs.div',
               'time_units': 'years' if time_tree else None,
               'time_tree_requested': time_requested,
               'time_tree_written': time_tree,
               'time_export_warning': time_warning,
               'country_comparator_created': False,
               'upload_performed': False}
    (outdir / 'export_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (outdir / 'README.md').write_text(
        '# iTOL files\n\n'
        f'Upload `{prefix}_divergence.nwk` for the divergence-scaled tree. '
        'Its branch lengths retain the source Auspice divergence scale.\n\n'
        + (f'`{prefix}_time_years.nwk` contains branch durations in years, not substitutions. '
           'It does not automatically configure a calendar-date axis in iTOL.\n\n' if time_tree else '')
        + (f'Temporal export was omitted: {time_warning} Check the original temporal inference before exporting branch durations; no dates were changed.\n\n' if time_warning else '')
        + 'Add the generated colour-strip `.txt` files after uploading a tree. '
        'Keep terminal identifiers unchanged. These datasets do not change the phylogeny. '
        'No data have been uploaded by this script.\n\n'
        'Review identifiers, metadata and data-use permissions before external upload.\n', encoding='utf-8')
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--auspice', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--prefix', default='flu_b')
    parser.add_argument('--columns', nargs='*', default=['country', 'clade'])
    parser.add_argument('--metadata', type=Path)
    parser.add_argument('--no-time-tree', action='store_true')
    args = parser.parse_args()
    try:
        summary = export(args.auspice, args.outdir, args.prefix, args.columns,
                         args.metadata, not args.no_time_tree)
    except (ValueError, OSError, json.JSONDecodeError) as error:
        parser.exit(1, f'Export error: {error}\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
