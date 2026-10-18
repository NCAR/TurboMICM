"""Rebuild and audit both compiler gates; run only after timing jobs finish."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command(args, output, cwd=ROOT):
    with output.open('w') as f:
        result = subprocess.run(args, cwd=cwd, stdout=f, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'{args} failed; see {output}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compilers', default='clang,gcc-15')
    args = parser.parse_args()
    (ROOT / 'scans').mkdir(exist_ok=True); (ROOT / 'build').mkdir(exist_ok=True)
    required = (ROOT.parent / 'reference/micm_cases.txt').read_text().splitlines()
    own = (ROOT / 'case-list.txt').read_text().splitlines()
    expected = [line.split() for line in required if line and not line.startswith('#')]
    assert [line.split() for line in own if line and not line.startswith('#')] == expected
    ids = {line[0] for line in expected}; assert len(ids) == 680
    reference = ROOT / 'results/tier1-reference-10000'
    reference.mkdir(parents=True, exist_ok=True)
    command(['make', 'micm-ref'], ROOT / 'build/validation-reference-build.log')
    cells = ROOT / 'results/large/tier1-30s/cells.csv'
    assert cells.exists(), 'run tools/large_scenarios.py first'
    command([str(ROOT / 'bin/micm_ts1_ref'), '--threads', '10', '--group-size', '32',
             '10000', '100', '30', str(cells)], reference / 'run.log', cwd=reference)
    report = dict(case_list_sha256=sha(ROOT.parent / 'reference/micm_cases.txt'),
                  reference_10000_sha256=sha(reference / 'ts1_final.csv'),
                  reference_input_sha256=sha(cells), compilers={})
    report_path = ROOT / 'scans/validation.json'
    for compiler in args.compilers.split(','):
        build_log = ROOT / 'build' / f'validation-build-{compiler}.log'
        command(['make', 'bench', f'CC={compiler}', '-j4'], build_log)
        case_log = ROOT / 'build' / f'validation-tier1-{compiler}.log'
        command(['make', 'test-cpu', f'CC={compiler}'], case_log)
        records = {}
        for line in case_log.read_text().splitlines():
            fields = line.split(' ', 2)
            if fields[0] in ids:
                assert fields[0] not in records
                records[fields[0]] = fields[1:]
        assert set(records) == ids
        assert all(status in ['PASS', 'N/A'] and reason for status, reason in records.values())
        counts = {status: sum(value[0] == status for value in records.values()) for status in ['PASS', 'N/A', 'FAIL']}
        assert counts == {'PASS': 675, 'N/A': 5, 'FAIL': 0}
        with (ROOT / 'scans' / f'tier1-{compiler}.csv').open('w', newline='') as f:
            writer = csv.writer(f, lineterminator='\n'); writer.writerow(['identifier', 'kind', 'status', 'reason'])
            for identifier, kind in expected:
                writer.writerow([identifier, kind, *records[identifier]])
        tier3_log = ROOT / 'build' / f'validation-tier3-{compiler}.log'
        command(['make', 'test-tier3', f'CC={compiler}'], tier3_log)
        verification = ROOT / 'prototypes' / ('verification-clang.csv' if compiler == 'clang' else 'verification-gcc15.csv')
        verification.write_bytes((ROOT / 'build/prototype-verification.csv').read_bytes())
        prototype_rows = list(csv.DictReader(verification.open()))
        assert len(prototype_rows) == 45
        gates = []
        for count, input_path, ref_path in [
            (100, ROOT.parent / 'reference/ts1_cells_100_seed1.csv', ROOT.parent / 'reference/micm_ts1_final_100_seed1.csv'),
            (10000, cells, reference / 'ts1_final.csv')]:
            out = ROOT / 'results' / f'gate-{compiler}-{count}'
            gate_log = ROOT / 'build' / f'validation-gate-{compiler}-{count}.log'
            command([sys.executable, str(ROOT / 'tools/check_benchmark.py'), str(count), '100', '30',
                     str(input_path), str(ref_path), str(out), '--compiler', compiler] + (['--threads', '1'] if count == 100 else []), gate_log)
            metadata = json.loads((out / 'metadata.json').read_text())
            assert metadata['gate_exit_code'] == 0
            gates.append(dict(cells=count, input_sha256=sha(input_path), reference_sha256=sha(ref_path),
                              result_sha256=sha(out / 'ts1_final.csv'), gate=metadata['gate'], stdout=metadata['stdout']))
        report['compilers'][compiler] = dict(
            tier1=counts, tier1_log_sha256=sha(case_log), tier3_log=tier3_log.read_text(),
            prototype_cases=len(prototype_rows), prototype_max_absolute_difference=max(float(r['max_absolute_difference']) for r in prototype_rows),
            flags=(ROOT / 'build' / compiler / 'flags').read_text().strip(),
            scenario_binary_sha256=sha(ROOT / 'bin/ts1_scenario'), gates=gates)
        report_path.write_text(json.dumps(report, indent=2)+'\n')
        print(compiler, counts, '45 prototype cases, scenario checks, 100- and 10000-cell roundoff gates PASS', flush=True)


if __name__ == '__main__':
    main()
