"""Native runner contract: restarts, forcing interpolation and CSV snapshots."""
import csv
import subprocess
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT/'runs'/'scenario-contract'
DIRECTORY.mkdir(parents=True,exist_ok=True)


def read(path):
    with path.open() as file:
        return list(csv.DictReader(file))


def write(path,rows):
    with path.open('w',newline='') as file:
        writer = csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def run(name,rows,photo,dt,steps,every=0,options=(),threads=1,method='rosenbrock3'):
    directory = DIRECTORY/name
    directory.mkdir(exist_ok=True)
    write(directory/'cells.csv',rows)
    (directory/'photolysis.csv').write_text('time_s,factor\n'+''.join(f'{t},{f}\n' for t,f in photo))
    (directory/'scenario.txt').write_text(f'# parser contract\n\n cells = cells.csv # relative\n'
        f' photolysis = photolysis.csv\n dt = {dt}\nsteps = {steps}\nsnapshot_every = {every}\n')
    command = [str(ROOT/'bin/ts1_scenario'),'--threads',str(threads),'--method',method]
    for option in options:
        command.extend(['--option',option])
    result = subprocess.run(command+['scenario.txt'],cwd=directory,capture_output=True,text=True)
    assert result.returncode == 0,result.stderr
    assert 'elapsed_ms=' in result.stdout and 'internal_steps=' in result.stdout
    return read(directory/'ts1_final.csv'), read(directory/'ts1_snapshots.csv')


def close(actual,expected):
    assert len(actual) == len(expected)
    for a,b in zip(actual,expected):
        assert set(a) == set(b)
        for key in a:
            x,y = float(a[key]),float(b[key])
            assert abs(x-y) <= 3e-13*abs(y)+1e-30,(key,x,y)


def main():
    rows = read(ROOT.parent/'reference/ts1_cells_100_seed1.csv')[:3]
    final,snapshots = run('constant',rows,[(0,1)],30,3,2)
    assert [(int(row['step']),int(row['cell'])) for row in snapshots] == [
        (step,cell) for step in (2,3) for cell in range(3)]
    # Tier 1's warmup plus two steps must equal three scenario steps.
    benchmark = DIRECTORY/'benchmark'; benchmark.mkdir(exist_ok=True)
    result = subprocess.run([str(ROOT/'bin/ts1_bench_cpu'),'--threads','1','3','2','30',
                             str(DIRECTORY/'constant/cells.csv')],cwd=benchmark,capture_output=True,text=True)
    assert result.returncode == 0,result.stderr
    close(final,read(benchmark/'ts1_final.csv'))
    photo = [(0,0),(60,1)]
    whole,snaps = run('linear',rows,photo,30,4,3)
    state = rows
    chained = []
    for k,factor in enumerate([0,.5,1,1]):
        result,_ = run(f'piece-{k}',state,[(0,factor)],30,1)
        chained.append(result)
        # Keep the original parameter values while replacing concentrations.
        state = [old|new for old,new in zip(rows,result)]
    close(whole,chained[-1])
    for step in (3,4):
        snapshot = [{k:v for k,v in row.items() if k not in ('step','cell')}
                    for row in snaps if int(row['step']) == step]
        close(snapshot,chained[step-1])
    for method in ('rosenbrock3','rodas3','rodas31','sdirk2'):
        one,_ = run(method+'-one',rows,photo,30,4,options=['group_size=1'],method=method)
        three,_ = run(method+'-three',rows,photo,30,4,options=['group_size=1'],threads=3,method=method)
        close(one,three)
    invalid = [['--method','unknown'],['--threads','0'],['--option','bogus=1'],
               ['--option','jacobian_steps=2'], # adaptive option with legacy method
               ['--method','rodas3','--option','h_start=0'],
               ['--option','max_steps=0']]
    for options in invalid:
        result = subprocess.run([str(ROOT/'bin/ts1_scenario'),*options,
                                 str(DIRECTORY/'constant/scenario.txt')],cwd=DIRECTORY,capture_output=True)
        assert result.returncode != 0,options
    print('PASS: default/Tier 1 equivalence, piecewise photolysis and boundaries, final/nonmultiple snapshots, threads, invalid CLI')


if __name__ == '__main__':
    main()
