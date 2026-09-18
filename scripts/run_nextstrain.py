#!/usr/bin/env python3
"""Run a consensus-sequence-to-Auspice/iTOL workflow for one influenza segment."""
from __future__ import annotations
import argparse
from dataclasses import dataclass
import datetime as dt
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from typing import Any

REPO = Path(__file__).resolve().parents[1]


@dataclass
class Step:
    name: str
    argv: list[str]
    outputs: list[Path]


def absolute(path: str | Path) -> Path:
    path = Path(path)
    return path if path.is_absolute() else REPO / path


def load_configuration(path: Path, segment: str) -> tuple[dict, dict]:
    config = json.loads(path.read_text())
    if config.get('schema_version') != 1:
        raise ValueError('Expected configuration schema_version 1.')
    if segment not in config.get('segments', {}):
        raise ValueError(f'Segment {segment!r} is not present in the configuration.')
    settings = dict(config.get('defaults', {}))
    settings.update(config['segments'][segment])
    for required in ['sequences', 'metadata', 'reference_genbank', 'reference_tip', 'gene']:
        if not settings.get(required):
            raise ValueError(f'Missing configuration field: {required}.')
    return config, settings


def build_steps(config: dict, s: dict, segment: str, threads: int) -> list[Step]:
    work = absolute(config.get('output_root', 'work')) / segment
    prep = work / 'prepared'
    ref = absolute(s['reference_genbank'])
    tip = s['reference_tip']
    py = sys.executable
    augur = config.get('augur_executable', 'augur')
    seq, meta = prep / 'sequences.fasta', prep / 'metadata.tsv'
    filtered, fmeta = work / 'filtered.fasta', work / 'metadata.tsv'
    aligned, rawtree, tree = work / 'aligned.fasta', work / 'tree_raw.nwk', work / 'tree.nwk'
    bl, nt, aa, tr = [work / name for name in ['branch_lengths.json','nt_muts.json','aa_muts.json','traits.json']]
    aus = work / 'auspice' / f'flu_b_{segment}.json'
    steps = []
    prepare = [py, str(REPO / 'scripts/prepare_phylogeny_inputs.py'),
               '--sequences', str(absolute(s['sequences'])), '--metadata', str(absolute(s['metadata'])),
               '--reference-genbank', str(ref), '--reference-tip', tip,
               '--gene', s['gene'], '--outdir', str(prep), '--trait-columns'] + s.get('trait_columns', [])
    if s.get('max_ambiguous_fraction') is not None:
        prepare += ['--max-ambiguous-fraction', str(s['max_ambiguous_fraction'])]
    if not s.get('require_exact_ids', True):
        prepare += ['--allow-unmatched-sequences']
    steps.append(Step('prepare', prepare, [seq, meta, prep / 'metadata_export.tsv']))
    index = work / 'sequence_index.tsv'
    steps.append(Step('index', [augur,'index','--sequences',str(seq),'--output',str(index)], [index]))
    filt = [augur,'filter','--sequences',str(seq),'--metadata',str(meta),
            '--sequence-index',str(index),'--metadata-id-columns','strain',
            '--output-sequences',str(filtered),'--output-metadata',str(fmeta),
            '--output-strains',str(work / 'filtered_strains.txt'),
            '--output-log',str(work / 'filter_log.tsv'),'--non-nucleotide']
    for key, flag in [('min_length','--min-length'), ('max_length','--max-length'),
                      ('min_date','--min-date'),('max_date','--max-date'),
                      ('exclude_ambiguous_dates_by','--exclude-ambiguous-dates-by')]:
        if s.get(key) is not None:
            filt += [flag, str(s[key])]
    include_files = []
    if s.get('retain_reference', True):
        include_files.append(str(prep / 'include_reference.txt'))
    if s.get('include_file'):
        include_files.append(str(absolute(s['include_file'])))
    if include_files:
        filt += ['--include'] + include_files
    if s.get('exclude_file'):
        filt += ['--exclude', str(absolute(s['exclude_file']))]
    steps.append(Step('filter', filt, [filtered, fmeta, work / 'filtered_strains.txt']))
    steps.append(Step('align', [augur,'align','--sequences',str(filtered),'--reference-name',tip,
                             '--method','mafft','--output',str(aligned),'--nthreads',str(threads)], [aligned]))
    model = s.get('substitution_model') or 'SET_SUBSTITUTION_MODEL'
    tree_cmd = [augur,'tree','--alignment',str(aligned),'--method','iqtree',
                '--substitution-model',model,'--nthreads',str(threads),'--output',str(rawtree)]
    extra = s.get('tree_builder_args', '')
    if s.get('tree_seed') is not None:
        extra = f'{extra} -seed {int(s["tree_seed"])}'.strip()
    if extra:
        tree_cmd += ['--tree-builder-args=' + extra]
    steps.append(Step('tree', tree_cmd, [rawtree]))
    refine = [augur,'refine','--tree',str(rawtree),'--alignment',str(aligned),
              '--metadata',str(fmeta),'--metadata-id-columns','strain','--timetree',
              '--root', s.get('root') or tip, '--output-tree',str(tree),
              '--output-node-data',str(bl),'--divergence-units','mutations-per-site',
              '--branch-length-inference',s.get('branch_length_inference','auto'),
              '--date-inference',s.get('date_inference','joint')]
    for key, flag in [('clock_filter_iqd','--clock-filter-iqd'),('clock_rate','--clock-rate'),
                      ('clock_std_dev','--clock-std-dev'),('coalescent','--coalescent'),('seed','--seed')]:
        if s.get(key) is not None:
            refine += [flag, str(s[key])]
    if s.get('date_confidence',True):
        refine += ['--date-confidence']
    if s.get('keep_polytomies',False):
        refine += ['--keep-polytomies']
    steps.append(Step('refine', refine, [tree, bl]))
    ancestral = [augur,'ancestral','--tree',str(tree),'--alignment',str(aligned),
                 '--root-sequence',str(ref),'--output-node-data',str(nt),
                 '--inference',s.get('ancestral_inference','joint')]
    if s.get('keep_ambiguous',True):
        ancestral += ['--keep-ambiguous']
    if s.get('keep_overhangs',True):
        ancestral += ['--keep-overhangs']
    if s.get('seed') is not None:
        ancestral += ['--seed',str(s['seed'])]
    steps.append(Step('ancestral', ancestral, [nt]))
    steps.append(Step('translate', [augur,'translate','--tree',str(tree),
        '--ancestral-sequences',str(nt),'--reference-sequence',str(ref),'--genes',s['gene'],
        '--output-node-data',str(aa),'--alignment-output',str(work / 'aa_%GENE.fasta')],
        [aa,work / f'aa_{s["gene"]}.fasta']))
    node_data = [bl, nt, aa]
    traits = s.get('trait_columns', [])
    if traits:
        steps.append(Step('traits', [augur,'traits','--tree',str(tree),'--metadata',str(fmeta),
            '--metadata-id-columns','strain','--columns'] + traits +
            ['--confidence','--output-node-data',str(tr)], [tr]))
        node_data.append(tr)
    steps.append(Step('export', [augur,'export','v2','--tree',str(tree),
        '--metadata',str(prep / 'metadata_export.tsv'),'--metadata-id-columns','strain',
        '--node-data'] + [str(p) for p in node_data] +
        ['--auspice-config',str(absolute(config['auspice_config'])),
         '--title',f'Influenza B/Victoria {segment}','--output',str(aus)], [aus]))
    steps.append(Step('validate', [augur,'validate','export-v2',str(aus)], []))
    steps.append(Step('itol', [py,str(REPO / 'scripts/export_itol.py'),'--auspice',str(aus),
        '--metadata',str(prep / 'metadata_observed.tsv'),'--outdir',str(work / 'itol'),
        '--prefix',segment,'--columns'] + s.get('itol_columns',['country','clade']),
        [work / 'itol' / f'{segment}_divergence.nwk',work / 'itol' / 'export_summary.json']))
    if s.get('run_patristic',False):
        steps.append(Step('patristic', [py,str(REPO / 'scripts/patristic_distance.py'),str(aus),
            '--reference',s.get('patristic_reference') or tip,'--out',str(work / 'patristic_distances.csv')],
            [work / 'patristic_distances.csv']))
    return steps


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def tool_info(command: list[str]) -> dict:
    resolved = shutil.which(command[0])
    if not resolved:
        raise ValueError(f'Required executable is unavailable: {command[0]}. Activate the workflow environment.')
    result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    if result.returncode != 0:
        raise ValueError(f'Version check failed: {shlex.join(command)}\n{result.stderr[:500]}')
    return {'path': resolved, 'version_output': (result.stdout + result.stderr).strip()}


def check_settings(s: dict) -> None:
    if s.get('parameters_confirmed') is not True:
        raise ValueError('Set parameters_confirmed to true only after entering the study input paths, QC thresholds, model and clock settings in the configuration.')
    if not s.get('substitution_model'):
        raise ValueError('Set the segment substitution_model from the original IQ-TREE command/log.')
    for key in ('seed', 'tree_seed'):
        if s.get(key) is not None and (isinstance(s[key],bool) or not isinstance(s[key],int) or s[key] < 0):
            raise ValueError(f'{key} must be null or a nonnegative integer.')
    if s.get('clock_filter_iqd') is not None and s['clock_filter_iqd'] <= 0:
        raise ValueError('clock_filter_iqd must be positive or null (disabled).')
    # Browser exports are deliberately limited to these nonsensitive annotation fields.
    allowed = {'country','region','clade'}
    if not set(s.get('trait_columns',[])).issubset(allowed):
        raise ValueError('Trait columns must be country, region or clade. Clinical metadata are not exported by this workflow.')
    if not set(s.get('itol_columns',[])).issubset(allowed):
        raise ValueError('iTOL columns must be country, region or clade.')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=REPO / 'config/analysis.json')
    parser.add_argument('--segment', choices=['HA','NA'], required=True)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--dry-run', action='store_true', help='Print commands; do not run tools or read study input files.')
    parser.add_argument('--force', action='store_true', help='Archive an existing segment output directory before a fresh build.')
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive.')
    try:
        config, settings = load_configuration(args.config, args.segment)
        steps = build_steps(config, settings, args.segment, args.threads)
        if args.dry_run:
            print('# Command plan only. Input files and tools have not been validated.')
            if not settings.get('parameters_confirmed'):
                print('# Study parameter confirmation is required before execution.')
            for step in steps:
                print(f'\n# {step.name}\n{shlex.join(step.argv)}')
            return
        check_settings(settings)
        inputs = {key: absolute(settings[key]) for key in ['sequences','metadata','reference_genbank']}
        for key in ['include_file','exclude_file']:
            if settings.get(key): inputs[key] = absolute(settings[key])
        inputs['auspice_config'] = absolute(config['auspice_config'])
        for path in inputs.values():
            if not path.is_file(): raise ValueError(f'Required input not found: {path}')
        tools = {'augur': tool_info([config.get('augur_executable','augur'),'--version']),
                 'mafft': tool_info(['mafft','--version'])}
        iq = 'iqtree2' if shutil.which('iqtree2') else 'iqtree'
        tools['iqtree'] = tool_info([iq,'--version'])
        expected = config.get('expected_augur_version')
        if expected and not re.search(rf'(?<!\d){re.escape(expected)}(?!\d)', tools['augur']['version_output']):
            raise ValueError(f'Augur version differs from configured {expected}. Update the environment or record an intentional version change in config.')
        output_root = absolute(config.get('output_root','work')).resolve()
        work = output_root / args.segment
        if output_root == REPO or REPO in [work] or output_root in [Path('/'),Path.home()]:
            raise ValueError('Choose a dedicated output_root directory.')
        for path in inputs.values():
            if work in path.resolve().parents:
                raise ValueError('Input files cannot be inside the segment output directory.')
        if work.exists() and any(work.iterdir()):
            if not args.force:
                raise ValueError(f'{work} already contains outputs. Use a new output_root or --force to archive and rebuild.')
            stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
            work.rename(work.with_name(work.name + '.backup.' + stamp))
        for directory in [work / 'logs', work / 'auspice', work / 'itol']:
            directory.mkdir(parents=True,exist_ok=True)
        manifest = {'segment':args.segment,'status':'running','started_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
                    'config':config,'effective_settings':settings,'threads':args.threads,
                    'input_files':{key:{'path':str(path),'sha256':sha256(path)} for key,path in inputs.items()},
                    'tools':tools,'python':sys.version,'steps':[],
                    'note':'A successful run records this execution; it is not evidence of an earlier manuscript analysis.'}
        versions={}
        for pkg in ['nextstrain-augur','phylo-treetime','biopython','numpy','pandas','scipy']:
            try: versions[pkg]=importlib.metadata.version(pkg)
            except importlib.metadata.PackageNotFoundError: pass
        manifest['python_packages']=versions
        manifest_path=work/'run_manifest.json'
        def save(): manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
        save()
        with (work/'commands.sh').open('w') as handle:
            handle.write('#!/usr/bin/env bash\nset -euo pipefail\n')
            handle.write('cd '+shlex.quote(str(REPO))+'\n')
            for step in steps: handle.write(shlex.join(step.argv)+'\n')
        try:
            for step in steps:
                print(f'[{args.segment}] {step.name}',flush=True)
                logfile=work/'logs'/f'{step.name}.log'
                for output in step.outputs: output.parent.mkdir(parents=True,exist_ok=True)
                item={'name':step.name,'command':step.argv,'log':str(logfile),'status':'running'}
                manifest['steps'].append(item); save()
                with logfile.open('w') as log:
                    log.write('$ '+shlex.join(step.argv)+'\n\n');log.flush()
                    result=subprocess.run(step.argv,cwd=REPO,stdout=log,stderr=subprocess.STDOUT)
                item['return_code']=result.returncode
                if result.returncode:
                    item['status']='failed'; raise RuntimeError(f'Step {step.name} failed. See {logfile}')
                for output in step.outputs:
                    if not output.is_file() or output.stat().st_size == 0:
                        item['status']='failed';raise RuntimeError(f'Expected output is missing/empty: {output}')
                item['status']='completed';save()
            manifest['status']='completed'
            (work/'build_complete.txt').write_text('Completed '+dt.datetime.now(dt.timezone.utc).isoformat()+'\n')
        except BaseException:
            manifest['status']='failed';save();raise
        manifest['finished_utc']=dt.datetime.now(dt.timezone.utc).isoformat();save()
        print(f'Completed {args.segment}: {work}')
    except (ValueError,OSError,RuntimeError,subprocess.SubprocessError) as error:
        parser.exit(1,f'Workflow error: {error}\n')


if __name__ == '__main__':
    main()
