"""Reproduce compiler gates, clean builds, timing, and supplied scenario scores.

Run with .conda/bin/python from this checkout. All writes stay in runs/ and
validation/. The sibling reference, scenarios, and micm trees are read-only.
Benchmarks execute sequentially, without simultaneous builds or tests.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT/'runs/tier2'
REFERENCE = ROOT.parent/'reference'
SCENARIOS = ROOT.parent/'scenarios'


def run(command, cwd=ROOT, log=None):
    result = subprocess.run(list(map(str, command)),cwd=cwd,text=True,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if log:
        Path(log).write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f'command failed ({result.returncode}): {command}\n{result.stdout[-6000:]}')
    return result.stdout


def stats(text):
    return {key: float(value) if key == 'elapsed_ms' else int(value)
            for key,value in re.findall(r'(\w+)=([\d.]+)',text)}


def score(result, reference, level=None):
    command = [sys.executable,SCENARIOS/'score.py',result,reference,'--json']
    if level is not None:
        command += ['--level',str(level)]
    return json.loads(run(command))


def compare(result, reference):
    return run([sys.executable,REFERENCE/'compare_final.py',result,reference]).strip()


def supplementary():
    dest = RUNS/'supplementary'
    dest.mkdir(parents=True,exist_ok=True)
    run([sys.executable,SCENARIOS/'make_scenario.py',dest/'summer-48h','--cells','65',
         '--seed','4242','--dt','900','--steps','192','--diurnal','--lat','55',
         '--doy','172','--start','19','--nox-scale','3','--snapshot-every','24'])
    interpolation = dest/'interpolation'
    interpolation.mkdir(exist_ok=True)
    with (RUNS/'cells.csv').open() as source, (interpolation/'cells.csv').open('w') as output:
        reader, writer = csv.reader(source),csv.writer(output)
        for _ in range(34):
            writer.writerow(next(reader))
    (interpolation/'photolysis.csv').write_text('time_s,factor\n600,0.2\n1234,1\n2200,0.4\n')
    (interpolation/'scenario.txt').write_text('cells = cells.csv # relative path\ndt = 450\nsteps = 8\n'
                                             'photolysis = photolysis.csv\nsnapshot_every = 3\n')
    return [dest/'summer-48h/scenario.txt', interpolation/'scenario.txt']


def validate(compiler, repetitions, extra):
    label = 'gcc' if 'gcc' in compiler else 'clang'
    dest = RUNS/label
    dest.mkdir(parents=True,exist_ok=True)
    record = dict(compiler=run([compiler,'--version']).splitlines()[0],
                  platform=platform.platform(),optimization='-O3 -DNDEBUG',group_size=32,
                  benchmark=dict(cells=10000,seed=2026,steps=100,dt=30),timings=[],scenarios={})
    output = ROOT/'validation'/f'{label}-tier2.json'
    def save():
        output.write_text(json.dumps(record,indent=2)+'\n')
    make = ['make','-j6',f'CC={compiler}']
    print(f'{label}: clean build',flush=True)
    run(['make','clean'])
    start = time.perf_counter()
    run(make+['bench'],log=dest/'build.log')
    record['clean_bench_build_seconds'] = time.perf_counter()-start
    sources = sorted((ROOT/'build/generated').glob('ts1_specialized*.c'))
    record['generated_c_bytes'] = sum(p.stat().st_size for p in sources)
    record['generated_sources'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    save()
    print(f'{label}: CPU and generator tests',flush=True)
    cpu = run(make+['test-cpu'],log=dest/'cpu.log')
    cases = [line for line in cpu.splitlines() if line.endswith(' PASS')]
    expected = [line.split()[0] for line in (ROOT/'tests/cases.txt').read_text().splitlines()
                if line and not line.startswith('#')]
    assert len(cases) == len(expected) == 680
    assert [line.removesuffix(' PASS') for line in cases] == expected
    (ROOT/'validation'/f'{label}-cases.txt').write_text('\n'.join(cases)+'\n')
    record['cpu_cases'] = len(cases)
    specialized = run(make+['test-specialized'],log=dest/'specialized.log')
    record['specialized_tests'] = [s for s in specialized.splitlines() if s.startswith('PASS ')]
    assert len(record['specialized_tests']) == 6
    save()
    print(f'{label}: small MICM gates',flush=True)
    record['small_gates'] = {}
    for program in ('ts1_bench_cpu','ts1_bench_specialized'):
        work = dest/(program+'-small')
        work.mkdir(exist_ok=True)
        run([ROOT/'bin'/program,'--threads','1','100','100','30',REFERENCE/'ts1_cells_100_seed1.csv'],cwd=work)
        record['small_gates'][program] = compare(work/'ts1_final.csv',REFERENCE/'micm_ts1_final_100_seed1.csv')
    reference_final = RUNS/'clang/micm_ts1_ref-t18-r0/ts1_final.csv'
    if label == 'clang':
        run(make+['micm-ref','micm-scenario'],log=dest/'reference-build.log')
        programs = ('micm_ts1_ref','ts1_bench_cpu','ts1_bench_specialized')
    else:
        assert reference_final.exists(), 'validate clang before GCC to create the full reference'
        programs = ('ts1_bench_cpu','ts1_bench_specialized')
    for threads in (18,1):
        for program in programs:
            for repeat in range(repetitions if threads == 18 else 1):
                print(f'{label}: {program}, {threads} threads, sample {repeat+1}',flush=True)
                work = dest/f'{program}-t{threads}-r{repeat}'
                work.mkdir(exist_ok=True)
                text = run([ROOT/'bin'/program,'--threads',threads,'10000','100','30',RUNS/'cells.csv'],
                           cwd=work,log=work/'timing.txt')
                item = dict(program=program,threads=threads,sample=repeat,**stats(text))
                item['gate'] = compare(work/'ts1_final.csv',reference_final)
                record['timings'].append(item)
                save()
    for path in sorted(SCENARIOS.glob('*/scenario.txt'))+extra:
        name = path.parent.name
        print(f'{label}: scenario {name}',flush=True)
        work = dest/name
        work.mkdir(exist_ok=True)
        text = run([ROOT/'bin/ts1_scenario_specialized','--threads','18','--group-size','32',
                    '--rtol','1e-6','--atol','1e-15',path],cwd=work,log=work/'timing.txt')
        if path in extra:
            reference = RUNS/'clang'/(name+'-micm')
            reference.mkdir(exist_ok=True)
            if label == 'clang':
                run([ROOT/'bin/micm_scenario','--threads','18','--group-size','32',
                     '--rtol','1e-6','--atol','1e-15',path],cwd=reference,log=reference/'timing.txt')
            ref_snap, ref_final = reference/'ts1_snapshots.csv',reference/'ts1_final.csv'
        else:
            ref_snap, ref_final = path.parent/'micm_snapshots.csv',path.parent/'micm_final.csv'
        item = dict(**stats(text),micm_snapshots=score(work/'ts1_snapshots.csv',ref_snap,10),
                    micm_final=score(work/'ts1_final.csv',ref_final,10))
        if path not in extra:
            item['reference_snapshots'] = score(work/'ts1_snapshots.csv',path.parent/'reference_snapshots.csv')
            item['reference_final'] = score(work/'ts1_final.csv',path.parent/'reference_final.csv')
        record['scenarios'][name] = item
        save()
    print(f'{label}: all gates passed; {output.name}',flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',action='append',help='default: clang then gcc-15')
    parser.add_argument('--repetitions',type=int,default=3)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error('positive repetitions required')
    RUNS.mkdir(parents=True,exist_ok=True)
    run([sys.executable,REFERENCE/'make_cells_csv.py','10000','2026','--output',RUNS/'cells.csv'])
    extra = supplementary()
    for compiler in args.compiler or ['clang','gcc-15']:
        validate(compiler,args.repetitions,extra)


if __name__ == '__main__':
    main()
