"""Regenerate final tables and work–precision figures from retained scan CSVs.

Requires completed final scans. Run after benchmarks so analysis does not
compete with timed integration. All plot caches and outputs stay in this repo.
"""
import csv
import hashlib
import json
import math
import os
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ['MPLCONFIGDIR'] = str(ROOT / 'build/matplotlib')
os.environ['XDG_CACHE_HOME'] = str(ROOT / 'build/cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAMES = ['tier1-30s', 'diurnal-30min', 'diurnal-10min', 'diurnal-5min', 'polluted-10min', 'clean-10min']
POINTS = json.loads((ROOT / 'scans/operating_points.json').read_text())
IDS = [p['id'] for p in POINTS]
LABELS = dict(zip(IDS, ['Fast', 'Balanced', 'Accurate', 'Recommended', 'Strict']))
LABELS.update({'c-micm3-matched': 'C MICM3', 'micm-micm3-matched': 'C++ MICM3'})
COLORS = dict(zip(IDS, ['#4379b5', '#22558b', '#e19b39', '#26967d', '#805ca1']))
OUT = ROOT / 'results'


def read(name):
    with (ROOT / 'scans' / (name + '.csv')).open(newline='') as stream:
        return list(csv.DictReader(stream))


def grouped(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row['config'], row['scenario'])].append(row)
    return groups


def criterion(row):
    return min(float(row['overall_digits']), float(row['worst_digits']) + 1)


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def savecsv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def verify(rows, count, repeats):
    assert len(rows) == count, (len(rows), count)
    for row in rows:
        assert row['status'] == 'ok' and row['gate'] == 'PASS' and row['level_result'] == 'PASS', row
    for key, group in grouped(rows).items():
        assert len(group) == repeats, key
        assert len(set(r['output_sha256'] for r in group)) == 1, key
        assert len(set(r['binary_sha256'] for r in group)) == 1, key
        for field in ['internal_steps', 'function_calls', 'jacobian_updates', 'decompositions', 'solves']:
            assert len(set(r[field] for r in group)) == 1, (key, field)


def summary(rows, cells):
    result = []
    for (config, scenario), group in grouped(rows).items():
        values = [float(r['elapsed_ms']) for r in group]
        first = group[0]
        result.append(dict(cells=cells, config=config, scenario=scenario, repeats=len(group),
                           median_ms=statistics.median(values), min_ms=min(values), max_ms=max(values),
                           overall_digits=first['overall_digits'], worst_species=first['worst_species'],
                           worst_digits=first['worst_digits'], level=first['level'],
                           level_result=first['level_result'], gate=first['gate'],
                           internal_steps=first['internal_steps'], function_calls=first['function_calls'],
                           jacobian_updates=first['jacobian_updates'], decompositions=first['decompositions'],
                           solves=first['solves'], rate_updates=first['rate_updates'],
                           clipped_mass=first['clipped_mass'], binary_sha256=first['binary_sha256']))
    return result


def limits(rows, config):
    selected = [r for r in rows if r['config'] == config]
    return min(float(r['overall_digits']) for r in selected), min(float(r['worst_digits']) for r in selected)


def plot_curves(rows, final):
    fig, axes = plt.subplots(2, 3, figsize=(13.3, 7.8), layout='constrained')
    methods = [('rodas3', 'RODAS3', '#4379b5'), ('extrap', 'Extrapolation 4', '#26967d'), ('micm', 'C MICM3', '#555555')]
    for ax, scenario in zip(axes.flat, NAMES):
        for method, label, color in methods:
            group = [r for r in rows if r['scenario'] == scenario and r['method'] == method and r['status'] == 'ok']
            group.sort(key=lambda r: -json.loads(r['options'])['rtol'])
            ax.plot([criterion(r) for r in group], [float(r['elapsed_ms']) for r in group], 'o-',
                    label=label, color=color, linewidth=1.5, markersize=4)
        for point in IDS:
            row = next(r for r in final if r['config'] == point and r['scenario'] == scenario)
            ax.scatter(criterion(row), row['median_ms'], marker='*', s=105, color=COLORS[point], edgecolor='white', linewidth=.5, zorder=4)
        baseline = next(r for r in final if r['config'] == 'micm-micm3-matched' and r['scenario'] == scenario)
        ax.scatter(criterion(baseline), baseline['median_ms'], marker='D', s=35, color='#b95056', label='C++ MICM3, median 5')
        for level in [2, 3]: ax.axvline(level, color='#aaaaaa', linestyle=':', linewidth=1)
        ax.set_title(scenario); ax.set_yscale('log'); ax.grid(axis='y', alpha=.18)
        ax.set_xlabel('min(overall, worst species + 1), digits'); ax.set_ylabel('Integration time (ms, log scale)')
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=4, frameon=False)
    fig.suptitle('1,000 cells · Clang · 10 threads · final driver\nTolerance curves: one run per point; stars: five-run medians of declared configurations', fontsize=13)
    for ext in ['png', 'svg']: fig.savefig(OUT / ('work_precision_1000.' + ext), dpi=170)
    plt.close(fig)


def plot_family(rows):
    fig, axes = plt.subplots(2, 3, figsize=(13.3, 7.8), layout='constrained')
    families = ['rodas3', 'ros3', 'extrap3', 'extrap4', 'sdirk2', 'sdirk3']
    for ax, scenario in zip(axes.flat, NAMES):
        for family in families:
            group = [r for r in rows if r['scenario'] == scenario and r['config'].startswith(family + '-') and r['status'] == 'ok']
            group.sort(key=lambda r: -json.loads(r['options'])['rtol'])
            ax.plot([criterion(r) for r in group], [float(r['elapsed_ms']) for r in group], 'o-', label=family, markersize=3, linewidth=1.2)
        ax.set_title(scenario); ax.set_yscale('log'); ax.grid(axis='y', alpha=.18)
        ax.set_xlabel('min(overall, worst species + 1), digits'); ax.set_ylabel('Integration time (ms, log scale)')
        for level in [2, 3]: ax.axvline(level, color='#aaaaaa', linestyle=':', linewidth=1)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=6, frameon=False)
    fig.suptitle('Family screen · 64 cells · 4 threads · stable tolerance sweep\nHistorical driver: complete rate recalculation and one external step per group visit', fontsize=13)
    for ext in ['png', 'svg']: fig.savefig(OUT / ('work_precision_families.' + ext), dpi=170)
    plt.close(fig)


def plot_scaling(rows):
    fig, ax = plt.subplots(figsize=(11.6, 5.3), layout='constrained')
    width = .15
    for j, point in enumerate(IDS):
        group = [next(r for r in rows if r['config'] == point and r['scenario'] == scenario) for scenario in NAMES]
        x = [i + (j-2)*width for i in range(6)]
        med = [r['median_ms']/1000 for r in group]
        err = [[(r['median_ms']-r['min_ms'])/1000 for r in group], [(r['max_ms']-r['median_ms'])/1000 for r in group]]
        ax.bar(x, med, width, label=LABELS[point], color=COLORS[point], yerr=err, capsize=2, error_kw={'linewidth':.7})
    ax.set_xticks(range(6), NAMES); ax.set_ylabel('Integration time (seconds)'); ax.grid(axis='y', alpha=.18)
    ax.set_title('10,000 cells · Clang · 10 threads\nMedian of five; error bars span minimum to maximum')
    ax.legend(ncol=5, frameon=False, loc='outside lower center' if False else 'upper left')
    for ext in ['png', 'svg']: fig.savefig(OUT / ('timing_10000.' + ext), dpi=170)
    plt.close(fig)


def species_report(rows):
    """All-species RMS and maximum errors, using the first final 1,000-cell run."""
    output = []
    for scenario in NAMES:
        record = next(r for r in rows if r['config'] == 'recommended' and r['scenario'] == scenario and r['repeat'] == '0')
        run = ROOT / record['run_log']
        suffix = 'final.csv' if scenario == 'tier1-30s' else 'snapshots.csv'
        with (ROOT.parent / 'scenarios' / scenario / ('reference_' + suffix)).open(newline='') as ref, (run.parent / ('ts1_' + suffix)).open(newline='') as result:
            ref_reader, result_reader = csv.reader(ref), csv.reader(result)
            header = next(ref_reader); assert next(result_reader) == header
            offset = 2 if suffix == 'final.csv' else 4
            count = [0] * (len(header)-offset); sums = [0.] * len(count); maximum = [0.] * len(count)
            for a, b in zip(ref_reader, result_reader, strict=True):
                assert a[:offset] == b[:offset]
                for j, (reference, value) in enumerate(zip(a[offset:], b[offset:])):
                    reference, value = float(reference), float(value)
                    if reference >= 1e-14:
                        error = abs(value-reference)/reference
                        sums[j] += error*error; count[j] += 1; maximum[j] = max(maximum[j], error)
            total = sum(sums) / sum(count)
            assert round(-math.log10(math.sqrt(total)),3) == float(record['overall_digits'])
            for j, name in enumerate(header[offset:]):
                if count[j]:
                    rms = math.sqrt(sums[j]/count[j]); digits = min(16., -math.log10(rms)) if rms else 16.
                    output.append(dict(scenario=scenario, species=name, values=count[j], rms_relative=rms, digits=digits, max_relative=maximum[j]))
    savecsv(OUT/'species_accuracy.csv', output)
    return output


def main():
    one, ten = read('operating_clang1000'), read('operating_clang10000')
    verify(one, 210, 5); verify(ten, 150, 5)
    gcc, stress_c, stress_g = read('final_gcc1000'), read('final_stress_clang'), read('final_stress_gcc')
    verify(gcc, 30, 1); verify(stress_c, 15, 1); verify(stress_g, 15, 1)
    expected_hash = json.loads((OUT/'validation-clang/build.json').read_text())['binaries']['ts1_scenario']
    assert all(r['binary_sha256'] == expected_hash for r in one+ten if r['config'] != 'micm-micm3-matched')
    small, large = summary(one,1000), summary(ten,10000)
    savecsv(OUT/'operating_summary.csv', small+large)
    indexed = {(r['cells'],r['config'],r['scenario']):r for r in small+large}
    totals = {(cells,point):sum(indexed[cells,point,s]['median_ms'] for s in NAMES) for cells,points in [(1000,list(LABELS)),(10000,IDS)] for point in points}
    species = species_report(one)
    plot_curves(read('precision_final1000'),small); plot_family(read('tolerance_final64')); plot_scaling(large)
    public_rows=[]; stress_rows=[]
    for point in IDS:
        lo = limits(one,point); hi = limits(ten,point); st = limits(stress_c+stress_g,point)
        public_rows.append([LABELS[point], f'{lo[0]:.3f}', f'{lo[1]:.3f}', f'{hi[0]:.3f}', f'{hi[1]:.3f}', f'{totals[1000,point]/1000:.3f}', f'{totals[10000,point]/1000:.3f}'])
        stress_rows.append([LABELS[point],f'{st[0]:.3f}',f'{st[1]:.3f}', 'PASS'])
    text = '''# Accuracy and performance

Fourth-order linearly implicit Euler extrapolation at `rtol=0.002`, `atol=1e-18`
is the recommended atmospheric-model configuration. It exceeds the target of
three overall digits and two digits for every scored species on every disclosed
case and on the additional 48-hour stress cases. RODAS3 provides cheaper
configurations for the two-digit tier. All five declared points pass their
intended level and finite/nonnegative gate in the final measurements.

## Scope and measurement

Measurements use the supplied 10-core Apple M5 host, macOS 27.0.1, Apple Clang
21.0.0 and GCC 15.2.0. All C objects use `-O3 -DNDEBUG -std=c11`, the documented
warnings/includes, and OpenMP. The supplied C++ MICM runner uses clang++ with
`-O3 -DNDEBUG -std=c++20`. No architecture, floating-point, LTO or profiling
flags are added. All reported times come from the compiled executables' internal
integration timers, including forcing/rate updates and validation, excluding
setup and CSV I/O. There is no scenario warm-up.

Final Clang measurements use ten threads and five serialized repetitions for
every configuration/scenario. No other benchmark or build from this task runs
concurrently; unrelated machine load is uncontrolled. Medians are local evidence,
not a promise of the evaluator's idle-machine timing. Each repetition has
bit-identical output and work counts within its configuration/scenario. Reused
scores require identical output SHA-256 and the same reference/level.

The original six disclosed inputs contain 1,000 cells each. The 10,000-cell
inputs reproduce their recipes and have a verified identical first 1,000 cells.
**Accuracy at 10,000 cells is measured only on that reference prefix.** The
finite/nonnegative gate covers all 10,000 cells, and the C runner checks every
external step for every cell. The larger input changes the grouping of the final
eight cells of the reference prefix, so prefix scores need not be identical.
The remaining 9,000 trajectories have no independent accuracy reference here.
No experimenter's withheld scenario was available.

## Declared points

The following accuracies are minima across the six scenarios, independently for
overall and worst species. Times sum the six per-scenario medians; this equal
weight summary is not an official ranking aggregate. Full commands and expected
levels are in [OPERATING_POINTS.md](OPERATING_POINTS.md).

'''
    text += table(['Point','Overall, 1k','Species, 1k','Overall, 10k prefix','Species, 10k prefix','Sum seconds, 1k','Sum seconds, 10k'],public_rows)+'\n\n'
    for cells in [1000,10000]:
        text += f'### {cells:,} cells: per-scenario median milliseconds\n\n'
        points = list(LABELS) if cells == 1000 else IDS
        text += table(['Scenario']+[LABELS[p] for p in points], [[s]+[f'{indexed[cells,p,s]["median_ms"]:.1f}' for p in points] for s in NAMES])+'\n\n'
    text += '![10,000-cell timing medians and ranges](results/timing_10000.png)\n\n'
    text += '''## MICM comparison and work–precision

The unchanged MICM3 defaults remain available. On the disclosed long runs their
`atol=1e-3` is too loose: the default fails on the 30-minute diurnal and polluted
cases and misses level three on other cases. Those runs remain in
`scans/baselines1000.csv`; their speed is not used as an accuracy-qualified baseline.
The original C and supplied C++ MICM3 solvers agree in their measured accuracy
when given matching explicit tolerances and group size.

The repeated matched baseline uses `rtol=0.001`, `atol=1e-18`, group size 32 and
100,000 maximum attempts. It passes level three on all six disclosed scenarios.
Looser tested MICM3 tolerances 0.003 and 0.01 do not pass that across-scenario
criterion. This is the fastest qualifying MICM3 setting in the tested grid,
not a proof that no intermediate tolerance would improve it. The C baseline
uses the same rate cache and scheduling as the new solvers; the supplied C++
runner is unchanged. Thus C-to-C speed ratios more directly isolate integration
methods, while C++ ratios also include library and driver differences.

'''
    for point in ['accurate','recommended']:
        c=totals[1000,'c-micm3-matched']/totals[1000,point]
        cpp=totals[1000,'micm-micm3-matched']/totals[1000,point]
        text += f'The **{LABELS[point].lower()}** point is {c:.2f}× faster than the C MICM3 baseline and {cpp:.2f}× faster than C++ MICM3 by summed scenario medians.\n\n'
    text += '''![Final full-domain work–precision curves](results/work_precision_1000.png)

The horizontal criterion is `min(overall digits, worst-species digits + 1)`;
reaching L is equivalent to the combined accuracy requirement. Curves follow
local tolerances, so nonmonotone segments are meaningful. Curves use single
runs with the final implementation; stars show the five-run declared-point
medians. C++ MICM3 is shown at the matched qualifying tolerance.

![Method-family work–precision screen](results/work_precision_families.png)

The 64-cell plot is the stable-build repetition of the tolerance screen. It uses
the historical driver before constant-rate caching and blocking; do not compare
its absolute times with the final driver plot. It shows why SDIRK2/3 and lower
extrapolation order were not selected. [SCANS.md](SCANS.md) describes controlled
structural comparisons, unsuccessful variants and provenance limitations.

## Stress tests and accuracy limits

The additional summer/winter/tropical 48-hour cases use 128 cells and different
seeds, forcing and pollution. Tight/tighter six-stage MICM references agree by
at least 10.230 overall digits and 9.483 worst-species digits. They informed
selection: the previously considered extrapolation `rtol=0.006` failed in the
tropical case and was replaced by `0.004`. These are stress tests, not untouched
holdouts. The table reports minima across all three cases and both compilers.

'''+table(['Point','Overall digits','Worst-species digits','Intended level'],stress_rows)+'\n\n'
    max_error=max(r['max_relative'] for r in species)
    worst=max(species,key=lambda r:r['max_relative'])
    text += f'The recommended point has a largest individual scored relative error of {100*max_error:.3f}% ({worst["species"]}, {worst["scenario"]}), illustrating that an RMS criterion is not a pointwise bound. All-species RMS and maximum errors are retained in `results/species_accuracy.csv`.\n\n'
    key_species=['O3','OH','HO2','NO','NO2','H2O2']
    text+=table(['Species','Minimum digits across six cases','Maximum individual relative error'],[[s,f'{min(r["digits"] for r in species if r["species"]==s):.3f}',f'{100*max(r["max_relative"] for r in species if r["species"]==s):.3f}%'] for s in key_species if any(r['species']==s for r in species)])+'\n\n'
    text+='''Projection is measured and charged to the local error norm, but does not exactly
preserve elemental invariants. A deliberately stiff conservative test documents
that effect; [ACCURACY.md](ACCURACY.md) states the budget and its limits. Reference
values below `1e-14 mol m-3` are excluded from the relative score but remain in
the chemistry and finite/nonnegative gate. The results do not bound unrecorded
instantaneous peaks or establish accuracy for every possible mechanism.

## Compatibility and verification

Both compilers pass all 680 original tests, with exactly the supplied case list
and no N/A entries, plus all nine original TS1 agreement variants at relative
`1e-10`. Both pass 44 independent prototype comparisons with matching operation
counts, ten fixed-order convergence checks, scalar/group equivalence, stiff
and invalid-input tests, generic rate-dependency tests, and scenario I/O tests.
The latter check photolysis interpolation/snapshots against MICM, generated versus
generic kernels, and bit-identical cache/block/schedule variants. New methods
are additions; the original solver formulas and mechanism remain unchanged.

The prototype uses independent dense operations and chemistry, sharing only the
explicitly disclosed ROS coefficient tables. Maximum C/prototype absolute
solution differences are below 128 double epsilons. Transcripts, coefficient
agreement results, orders and final executable hashes are in
`results/validation-clang/` and `results/validation-gcc/`. Every selected point
also passes its intended accuracy level with GCC on all disclosed and stress
cases. Portability is checked with both compilers on this host; Linux execution
was unavailable and is not claimed.

## Reproduction and evidence

`tools/validate.py --compiler clang` (or `gcc-15`) runs and records the compiler
gates. `tools/scan.py scans/final_comparison.json --scan NEW_NAME --repeats 5`
repeats the 1,000-cell comparison. Generate larger inputs with `tools/scaling.py`
and use `--scenario-root build/scaling` with `scans/operating_points.json`.
Always use `.conda/bin/python`, run scans sequentially, and do not rebuild during
measurements. `--discard-repeats --reuse-identical-scores` saves disk space while
retaining every timing/log and the first state of each configuration/scenario.

`results/operating_summary.csv` contains medians, ranges, accuracy and work;
`scans/` retains every configuration and raw measured row, including failures.
`tools/report.py` regenerates this report, operating points, and PNG/SVG figures.
The historical Tier 1 measurements are preserved separately in
[results/tier1_results.md](results/tier1_results.md); their workload includes a
warm-up and differs from the new scenario timer.
'''
    (ROOT/'RESULTS.md').write_text(text)
    op='''# Declared operating points

Five configurations are declared. Use all available cores by default, or add
`--threads N`; final local measurements use ten threads. All use the complete
mechanism, double precision, group size 32, `atol=1e-18`, fresh Jacobians,
proportional control, maximum-cell RMS error, warm step-size starts and no
intermediate-stage clipping. Rate caching is enabled with four-step blocks,
dynamic scheduling and the OpenMP default chunk. These are runtime parameters.

**Recommended: fourth-order Euler extrapolation, rtol=0.002.** It gives measured
margin above three overall digits and two digits for every scored species,
including the additional 48-hour cases. The accurate point costs less but has
less margin; the strict point provides additional margin at greater cost.

'''
    op+=table(['Point','Method','rtol','Expected level','Minimum overall/species digits, disclosed 1k'],[[LABELS[p['id']],p['method']+(' order 4' if p['method']=='extrap' else ''),str(p['options']['rtol']),str(p['level']),'/'.join(f'{v:.3f}' for v in limits(one,p['id']))] for p in POINTS])+'\n\n'
    for point in POINTS:
        op+='## '+LABELS[point['id']]+'\n\n```sh\nbin/ts1_scenario --method '+point['method']+' \\\n'
        options=list(point['options'].items())
        for start in range(0,len(options),3):
            op+='  '+' '.join('--option '+str(k)+'='+str(v) for k,v in options[start:start+3])+' \\\n'
        op+='  SCENARIO\n```\n\n'
    op+='''`SCENARIO` is the path to a scenario file. Machine-readable configurations are
in `scans/operating_points.json`; other method defaults are in DESIGN.md.
The numerical levels are expectations supported by measured scenarios, not
proofs for unseen forcing. Every declared point passes its intended level on
all disclosed and added stress cases with both compilers. The 10,000-cell gate
covers all cells, while independent accuracy covers only its first 1,000 cells.
See RESULTS.md for per-scenario accuracy, five-run timing summaries, work–precision
plots and the matched MICM comparison, and ACCURACY.md for scientific limitations.
'''
    (ROOT/'OPERATING_POINTS.md').write_text(op)
    manifest={name:hashlib.sha256((ROOT/'scans'/(name+'.csv')).read_bytes()).hexdigest() for name in ['operating_clang1000','operating_clang10000','final_gcc1000','final_stress_clang','final_stress_gcc','precision_final1000','tolerance_final64']}
    (OUT/'report_inputs.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Final reports and figures regenerated; all declared rows pass.')


if __name__ == '__main__':
    main()
