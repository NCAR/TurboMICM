"""Run native scenario configurations, recording the supplied scorer verbatim.

Run from the local conda environment. No Python execution is included in native
timings. Completed rows are flushed immediately; reruns append, never overwrite.
"""
import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = ['tier1-30s', 'diurnal-30min', 'diurnal-10min', 'diurnal-5min',
             'polluted-10min', 'clean-10min']
COUNTERS = ['internal_steps', 'accepted', 'rejected', 'function_calls', 'jacobian_updates',
            'decompositions', 'solves', 'jacobian_reuses', 'factor_reuses',
            'positivity_rejections', 'nonfinite_rejections', 'clipped_values',
            'nonlinear_iterations', 'nonlinear_failures', 'minimum_concentration',
            'clipped_amount', 'maximum_clip']
FIELDS = ['timestamp', 'id', 'scenario', 'repeat', 'compiler', 'source_sha256',
          'binary_sha256', 'method', 'options', 'threads', 'status', 'elapsed_ms',
          *COUNTERS, 'overall_digits', 'worst_species', 'worst_digits', 'gate',
          'score_json', 'stderr', 'output_sha256', 'accuracy_scope']


def run(config, scenario, repeat, args):
    name = scenario.name if scenario.is_dir() else scenario.parent.name
    scenario = (scenario/'scenario.txt' if scenario.is_dir() else scenario).resolve()
    directory = ROOT/'runs'/'scans'/config['id']/name/str(repeat)
    directory.mkdir(parents=True, exist_ok=True)
    binary = ROOT/'bin'/config.get('binary', 'ts1_scenario')
    options = config.get('options', {})
    command = [str(binary), '--threads', str(args.threads), '--method', config['method']]
    for key, value in options.items():
        if config.get('binary') == 'micm_scenario':
            command += ['--'+key.replace('_','-'), str(value)]
        else:
            command += ['--option', f'{key}={value}']
    command += [str(scenario)]
    result = subprocess.run(command, cwd=directory, capture_output=True, text=True)
    (directory/'stdout.txt').write_text(result.stdout)
    (directory/'stderr.txt').write_text(result.stderr)
    row = dict(timestamp=datetime.now(timezone.utc).isoformat(), id=config['id'], scenario=name,
               repeat=repeat, compiler=args.compiler, source_sha256=args.source,
               binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(), method=config['method'],
               options=json.dumps(options, sort_keys=True), threads=args.threads,
               status=result.returncode, stderr=result.stderr.strip())
    row.update({k:v for k,v in re.findall(r'(\w+)=([^\s]+)',result.stdout)
                if k in COUNTERS+['elapsed_ms']})
    if result.returncode == 0:
        has_snapshots = any(line.split('#')[0].strip().startswith('snapshot_every') and
                            int(line.split('#')[0].split('=')[1]) > 0
                            for line in scenario.read_text().splitlines())
        basename = 'snapshots' if has_snapshots else 'final'
        output = directory/f'ts1_{basename}.csv'
        reference = ((args.reference_root/name) if args.reference_root else scenario.parent)/f'reference_{basename}.csv'
        scored_output = output
        row['accuracy_scope'] = 'all supplied cells'
        if args.prefix_cells:
            scored_output = directory/'score-prefix.csv'
            with output.open() as source, scored_output.open('w',newline='') as target:
                reader = csv.DictReader(source)
                writer = csv.DictWriter(target,fieldnames=reader.fieldnames,lineterminator='\n')
                writer.writeheader()
                for i,values in enumerate(reader):
                    if int(values.get('cell',i)) < args.prefix_cells:
                        writer.writerow(values)
            row['accuracy_scope'] = f'first {args.prefix_cells} cells only'
        if reference.exists():
            score = subprocess.run([sys.executable, str(ROOT.parent/'scenarios'/'score.py'),
                                    str(scored_output),str(reference),'--json'],capture_output=True,text=True)
            if not score.stdout.strip():
                raise RuntimeError(f'Scorer failed: {score.stderr}')
            summary = json.loads(score.stdout)
            row.update(overall_digits=summary['overall_digits'],
                       worst_species=summary['worst_species'][0][0],
                       worst_digits=summary['worst_species'][0][1],gate=summary['gate'],
                       score_json=score.stdout.strip())
            (directory/'score.json').write_text(score.stdout)
        else:
            row['gate'] = 'UNSCORED'
        row['output_sha256'] = hashlib.sha256(output.read_bytes()).hexdigest()
    else:
        row['gate'] = 'RUN_FAILED'
    if not args.keep:
        for file in directory.glob('ts1_*.csv'):
            file.unlink()
        (directory/'score-prefix.csv').unlink(missing_ok=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path, help='JSON list of id, method, options objects')
    parser.add_argument('--log', type=Path, required=True)
    parser.add_argument('--scenarios', nargs='+', default=SCENARIOS)
    parser.add_argument('--scenario-root', type=Path, default=ROOT.parent/'scenarios')
    parser.add_argument('--reference-root',type=Path)
    parser.add_argument('--prefix-cells',type=int,default=0)
    parser.add_argument('--repeat', type=int, default=1)
    parser.add_argument('--threads', type=int, default=10)
    parser.add_argument('--compiler', default='clang')
    parser.add_argument('--keep', action='store_true')
    args = parser.parse_args()
    digest = hashlib.sha256()
    for file in sorted([*ROOT.glob('src/*.[ch]'), *ROOT.glob('include/*.h'),
                        *ROOT.glob('tools/*.c'), *ROOT.glob('tools/*.h'), ROOT/'Makefile']):
        digest.update(file.name.encode()); digest.update(file.read_bytes())
    args.source = digest.hexdigest()
    configs = json.loads(args.config.read_text())
    args.log.parent.mkdir(parents=True,exist_ok=True)
    is_new = not args.log.exists()
    with args.log.open('a',newline='') as file:
        writer = csv.DictWriter(file,fieldnames=FIELDS,lineterminator='\n')
        if is_new:
            writer.writeheader()
        for config in configs:
            for scenario in args.scenarios:
                for repeat in range(1,args.repeat+1):
                    row = run(config,args.scenario_root/scenario,repeat,args)
                    writer.writerow(row); file.flush()
                    print(config['id'],scenario,repeat,row.get('elapsed_ms','-'),
                          row.get('overall_digits','-'),row.get('worst_digits','-'),row['gate'],flush=True)


if __name__ == '__main__':
    main()
