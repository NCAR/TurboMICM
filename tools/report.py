"""Generate results tables and declared operating points from final native records."""
import csv
import argparse
import json
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['tier1-30s', 'diurnal-30min', 'diurnal-10min', 'diurnal-5min', 'polluted-10min', 'clean-10min']
POINTS = ['fast', 'balanced', 'accurate', 'recommended', 'precise']


def read(name):
    with (ROOT / 'scans' / name).open() as f:
        return list(csv.DictReader(f))


def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |'] +
                     ['| '+' | '.join(str(v) for v in row)+' |' for row in rows])


def subset(rows, config, scenario):
    return [r for r in rows if r['configuration'] == config and r['scenario'] == scenario]


def timing(rows, config, scenario):
    rs = subset(rows, config, scenario)
    return median(float(r['elapsed_ms']) for r in rs)


def verify(rows, count, repeats):
    assert len(rows) == len(POINTS)*len(NAMES)*repeats, (count, len(rows))
    for config in POINTS:
        for scenario in NAMES:
            rs = subset(rows, config, scenario)
            assert len(rs) == repeats and all(int(r['cells']) == count and r['status'] == 'OK' and r['gate'] == 'PASS' for r in rs)
            assert {int(r['repeat']) for r in rs} == set(range(repeats))
            assert len({r['output_sha256'] for r in rs}) == 1
            assert len({r['options'] for r in rs}) == 1
            for key in ['internal_steps', 'accepted', 'rejected', 'function_calls', 'jacobian_updates', 'decompositions', 'solves']:
                assert len({r[key] for r in rs}) == 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points-only', action='store_true')
    args = parser.parse_args()
    small = read('final-clang-1000.csv'); verify(small, 1000, 5)
    configs = json.loads((ROOT / 'scans/operating-points-throughput.json').read_text())
    low = {point: (min(float(r['overall_digits']) for r in small if r['configuration'] == point),
                   min(float(r['worst_digits']) for r in small if r['configuration'] == point)) for point in POINTS}
    for point in POINTS:
        level = 2 if point in ['fast', 'balanced'] else 3
        assert low[point][0] >= level and low[point][1] >= level-1
    lines = ['# Declared operating points', '',
             'These five configurations are the complete evaluation list. All use the full TS1 mechanism, '
             'double precision, fresh Jacobians, maximum cell RMS error control, and `atol=1e-18` mol m-3. '
             'Use all cores by omitting `--threads`; the measured runs used ten workers.', '',
             'Append the scenario path to each command below. All other parameters retain the defaults in `DESIGN.md`.', '']
    expected = []
    for config in configs:
        name = config['name']; options = config['options']
        command = 'bin/ts1_scenario --method '+config['method']
        for key, value in options.items():
            command += f' --option {key}={value:g}' if isinstance(value, (float, int)) else f' --option {key}={value}'
        lines += [f'## {name}'+(' (chosen for an atmospheric model)' if name == 'recommended' else ''), '',
                  '```sh', command+' SCENARIO', '```', '']
        expected.append([name, options['rtol'], f'{low[name][0]:.3f}', f'{low[name][1]:.3f}'])
    lines += [table(['Point', 'rtol', 'Lowest disclosed overall digits', 'Lowest disclosed species digits'], expected), '',
              'The recommended point supplies margin over the three-digit overall / two-digit species target, '
              'including the additional 48-hour validation. `accurate` is the cheapest declared point meeting that target '
              'on every disclosed scenario. `fast` is cheapest at level two but has little accuracy margin; `balanced` '
              'provides a more conservative level-two alternative. `precise` supports tighter sensitivity work.', '',
              'These are measured expectations on the disclosed cases, not promises about unknown inputs. '
              'The complete raw results, additional scenarios, projection measurements, and scope of the 10,000-cell '
              'reference checks are in `RESULTS.md`, `ACCURACY.md`, and `scans/README.md`. The inherited default solver '
              'and research alternatives remain callable but are not additional declared evaluation points.', '']
    (ROOT / 'OPERATING_POINTS.md').write_text('\n'.join(lines))
    if args.points_only:
        print('wrote operating points from complete full-cell validation')
        return
    large = read('final-clang-10000.csv'); verify(large, 10000, 5)
    baseline = read('micm-timings.csv'); stress = read('stress-final.csv')
    assert len(stress) == 15 and all(r['status'] == 'OK' and r['gate'] == 'PASS' for r in stress)
    assert all(float(r['overall_digits']) >= 3 and float(r['worst_digits']) >= 2
               for r in stress if r['configuration'] == 'recommended')
    validation = json.loads((ROOT / 'scans/validation.json').read_text())
    assert set(validation['compilers']) == {'clang', 'gcc-15'}
    lines = ['# Tier 3 results', '',
             'Measurements use the supplied 10-core machine, clang 21 and GCC 15.2, with '
             '`OPTFLAGS=-O3 -DNDEBUG`. Numerical methods solve all 547 reactions and 209 variable species. '
             'Only native C/C++ clocks supply elapsed times. Setup, validation, and CSV writing are excluded. '
             'No benchmark jobs or compiler jobs overlap the final timing batches.', '', '## Current validation', '']
    rows = []
    for compiler, info in validation['compilers'].items():
        assert info['tier1'] == {'PASS': 675, 'N/A': 5, 'FAIL': 0}
        rows.append([compiler, 675, 5, 0, info['prototype_cases'], f"{info['prototype_max_absolute_difference']:.3g}", 'PASS'])
    lines += [table(['Compiler', 'Tier 1 PASS', 'N/A', 'FAIL', 'Prototype comparisons', 'Maximum prototype difference', '100 / 10000-cell gate'], rows), '',
              'Every one of the 680 supplied identifiers is audited, including all eight plain programs. '
              'The five N/A cases are the retained compile-time C++ category/SFINAE cases. '
              '`scans/tier1-*.csv` records each result. All new-method prototype operation counts agree exactly; '
              'the listed concentration differences are accumulated double-precision roundoff. '
              'Scenario interpolation, photolysis caching, partial groups, snapshots, schedule independence, '
              'forcing discontinuities, ownership, invalid settings, and stiff nonnegative transfer are checked. '
              '`scans/validation.json` contains commands/flags, hashes, and the supplied all-cell/all-species roundoff gates.', '',
              '## Operating-point accuracy', '',
              'The following table uses all 1,000 disclosed cells and all reference snapshots. Each cell is '
              '`overall / worst-species` significant digits from the unmodified supplied scorer.', '']
    rows = [[name]+[f"{subset(small, point, name)[0]['overall_digits']} / {subset(small, point, name)[0]['worst_digits']}" for point in POINTS] for name in NAMES]
    lines += [table(['Scenario']+POINTS, rows), '',
              'All five points pass the nonfinite/negative gate on every run. `fast` and `balanced` meet level two '
              'on every disclosed scenario; the other three meet level three. `recommended` uses `rtol=3e-4`, '
              '16-cell groups, and dynamic scheduling. `ACCURACY.md` explains the scientific target, primary-species '
              'diagnostics, and measured small-negative projections.', '', '## Native elapsed time: 1,000 cells', '',
              'Values are milliseconds, medians of five complete scenario runs. The independent MICM3 baseline '
              'uses its exact default controls, or `rtol=3e-4, atol=1e-18` for the accuracy-controlled column. '
              'MICM uses groups of 32; the selected C points use groups of 16 with dynamic scheduling. '
              'The comparison includes both the retained C implementation improvements and numerical-control changes.', '']
    rows = []
    for name in NAMES:
        default = subset(baseline, 'MICM3-default', name)
        d = f"{timing(baseline, 'MICM3-default', name):.1f}" if all(r['status'] == 'OK' for r in default) else 'FAIL'
        controlled = timing(baseline, 'MICM3-controlled', name)
        rows.append([name]+[f'{timing(small, p, name):.1f}' for p in POINTS]+[d, f'{controlled:.1f}', f"{controlled/timing(small, 'recommended', name):.2f}×"])
    lines += [table(['Scenario']+POINTS+['MICM3 defaults', 'MICM3 controlled', 'Controlled / recommended'], rows), '',
              'A finite default MICM result is not necessarily accurate: the defaults fail on two scenarios, '
              'produce very inaccurate finite output on another, and miss level three on others. '
              'Default failures are not assigned a speedup. The controlled baseline qualifies at level three, '
              'but its achieved accuracy is generally higher than that threshold; ratios compare qualifying '
              'operating points rather than identical final errors.', '',
              '![Overall work–precision curves](plots/work_precision.png)', '',
              '![Worst-species work–precision curves](plots/worst_species_precision.png)', '',
              'Each curve uses measured median C time; bars show the observed minimum–maximum range. '
              'Dashed lines show the ranking thresholds. SVG versions are available beside the PNGs. '
              'The worst-species figure is necessary because the overall score alone is insufficient.', '',
              '## Native elapsed time: 10,000 cells', '',
              'These are independently measured complete 10,000-cell extensions, again medians of five native runs. '
              'All output values pass the gate. Relative-error scores cover the first 1,000 cells, which the '
              'generator verifies match the supplied inputs exactly; there is no claim of tight-reference '
              'accuracy for the remaining 9,000 cells.', '']
    lines += [table(['Scenario']+POINTS, [[name]+[f'{timing(large, p, name):.1f}' for p in POINTS] for name in NAMES]), '',
              'Full per-run times, counters, output hashes, and measured accuracy are in `scans/final-clang-10000.csv`. '
              '`scans/micm-large.csv` records independent single-run large MICM baseline comparisons. '
              'The experimenter supplies the authoritative idle-machine medians and hidden-scenario scores.', '',
              '## Additional 48-hour validation', '',
              'Three locally generated 64-cell cases vary seed, latitude, season, duration, step length, and NOx. '
              'Their MICM6 references at `rtol=1e-10, atol=1e-23` agree with looser independent runs by more '
              'than ten overall digits and nine worst-species digits. These are additional tests, not hidden '
              'experimenter inputs. Every declared point passes the gate.', '']
    names = sorted({r['scenario'] for r in stress})
    lines += [table(['Scenario']+POINTS, [[name]+[f"{subset(stress, p, name)[0]['overall_digits']} / {subset(stress, p, name)[0]['worst_digits']}" for p in POINTS] for name in names]), '',
              '`scans/stress-scenarios.json` and `tools/stress_references.py` reproduce these inputs and their '
              'reference-convergence checks. The recommended point also clears the stated accuracy target '
              'on these longer, different-regime tests.', '',
              '## Reproduction', '', '```sh', 'make env', 'make bench CC=clang',
              '.conda/bin/python tools/large_scenarios.py',
              '.conda/bin/python tools/scan.py scans/operating-points-throughput.json --study local-small --repeats 5 --metrics',
              '.conda/bin/python tools/scan.py scans/operating-points-throughput.json --study local-large --scenario-root results/large --repeats 5',
              'make micm-scenario',
              '.conda/bin/python tools/scan.py scans/micm-baselines.json --study local-micm --backend micm --compiler clang++ --repeats 5',
              '.conda/bin/python tools/validate.py', '```', '',
              'Use fresh study names for changed code/options. The supplied input directories remain read-only. '
              '`SCANS.md`, `REFERENCES.md`, and `scans/README.md` give the controlled comparisons, provenance, '
              'schema, and caveats. `tools/report.py` regenerates this report from the committed final study names. '
              'Historical Tier 1 measurements are retained in `RESULTS_TIER1.md`; the validation above is current.', '']
    (ROOT / 'RESULTS.md').write_text('\n'.join(lines))
    print('wrote RESULTS.md and OPERATING_POINTS.md from complete validated final records')


if __name__ == '__main__':
    main()
