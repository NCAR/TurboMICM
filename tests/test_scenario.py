"""Native scenario parser, photolysis, output ordering, and compatibility checks."""
import csv
from pathlib import Path
import subprocess
import sys

from api import C, I, MP, Model, State, lib, bind

ROOT = Path(__file__).resolve().parents[1]


def main():
    work = ROOT / 'build' / 'scenario-test'; work.mkdir(exist_ok=True)
    source = ROOT.parent / 'scenarios' / 'diurnal-10min' / 'cells.csv'
    with source.open(newline='') as f:
        reader = csv.DictReader(f); names = reader.fieldnames; cells = [next(reader) for _ in range(5)]
    with (work / 'cells.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=names); writer.writeheader(); writer.writerows(cells)
    (work / 'light.csv').write_text('time_s,factor\n0,0\n60,0.5\n120,1\n')
    (work / 'scenario.txt').write_text('# small independent exercise\ncells = cells.csv\ndt=30\nsteps=5\nphotolysis=light.csv\nsnapshot_every=2\n')
    command = [str(ROOT / 'bin/ts1_scenario'), '--threads', '2', '--option', 'group=3', str(work / 'scenario.txt')]
    result = subprocess.run(command, cwd=work, capture_output=True, text=True, check=True)
    assert 'elapsed_ms=' in result.stdout and 'internal_steps=' in result.stdout
    with (work / 'ts1_snapshots.csv').open() as f:
        snapshots = list(csv.DictReader(f))
    assert [(int(r['step']), int(r['cell'])) for r in snapshots] == [(s, c) for s in [2, 4, 5] for c in range(5)]
    bind('chem_load_ts1', MP, I)
    m = Model(pointer=lib.chem_load_ts1(1))
    expected = {}
    for first, count in [(0, 3), (3, 2)]:
        state = State(m, count)
        for c in range(count):
            for name, text in cells[first+c].items():
                value = float(text)
                if name in ['temperature', 'pressure']:
                    setattr(state.data.conditions[c], name, value)
                else:
                    v = lib.chem_variable_id(m.p, name.encode())
                    if v != C.c_size_t(-1).value:
                        state.data.variables[v*count+c] = value
                    else:
                        u = lib.chem_custom_id(m.p, name.encode()); state.data.custom[u*count+c] = value
            cond = state.data.conditions[c]
            cond.air_density = cond.pressure / (8.31446261815324 * cond.temperature)
        original = list(state.data.custom[:m.data.ncustom*count])
        for step, factor in enumerate([0., .25, .5, .75, 1.], 1):
            for u in range(m.data.ncustom):
                if m.data.custom_names[u].startswith(b'j'):
                    for c in range(count):
                        state.data.custom[u*count+c] = original[u*count+c]*factor
            # Recompute all rates here; the native driver caches thermal rates.
            state.rates(); state.solve(30.)
            if step in [2, 4, 5]:
                for c in range(count):
                    expected[step, first+c] = {m.data.species[m.data.variable_species[v]].name.decode(): state.data.variables[v*count+c]
                                              for v in range(m.data.nvariables)}
        state.close()
    for row in snapshots:
        for name, reference in expected[int(row['step']), int(row['cell'])].items():
            actual = float(row[name])
            assert abs(actual-reference) <= 1e-30 + 1e-10*abs(reference), (name, actual, reference)
    # Schedules change worker assignment, never group state or arithmetic.
    original_final = (work / 'ts1_final.csv').read_bytes()
    for schedule in ['dynamic', 'guided']:
        result = subprocess.run(command[:-1] + ['--option', f'schedule={schedule}', '--option', 'schedule_chunk=1', command[-1]],
                                cwd=work, capture_output=True, text=True, check=True)
        assert (work / 'ts1_final.csv').read_bytes() == original_final
    m.close()
    # Parser and option failures must be nonzero, including ignored legacy options.
    for extra in [['--method', 'unknown'], ['--option', 'group=0'], ['--option', 'rtol=nan'], ['--option', 'h_start=1']]:
        r = subprocess.run([str(ROOT / 'bin/ts1_scenario')] + extra + [str(work / 'scenario.txt')], cwd=work,
                           capture_output=True, text=True)
        assert r.returncode != 0
    print('native scenario test: interpolation, all-rate vs cached-rate agreement, partial groups, snapshots, CLI errors PASS')


if __name__ == '__main__':
    main()
