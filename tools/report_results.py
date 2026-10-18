"""Build the final measured report; missing or failed evidence stops generation."""
import csv
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = ['tier1-30s','diurnal-30min','diurnal-10min','diurnal-5min','polluted-10min','clean-10min']


def read(name):
    with (ROOT/'scans'/name).open() as file:
        return list(csv.DictReader(file))


def select(rows,id_,scenario):
    return [r for r in rows if r['id'] == id_ and r['scenario'] == scenario]


def main():
    points = json.loads((ROOT/'scans/operating-points.json').read_text())
    small,large,stress = read('final-1000.csv'),read('final-10000.csv'),read('final-stress.csv')
    baseline = read('final-baseline-10000.csv')
    gcc = read('gcc-operating-points.csv')+read('gcc-stress.csv')
    gates = json.loads((ROOT/'validation/tier3-gates.json').read_text())
    assert set(gates) == {'clang','gcc'} and all(g['cases_passed'] == 680 for g in gates.values())
    summary = {}
    for point in points:
        id_ = point['id']
        sr = [r for r in small+stress if r['id'] == id_]
        lr = [r for r in large if r['id'] == id_]
        gr = [r for r in gcc if r['id'] == id_]
        assert len(sr) == 9 and len(lr) == 30 and len(gr) == 9,id_
        assert all(r['method'] == point['method'] and json.loads(r['options']) == point['options']
                   for r in sr+lr+gr),id_
        assert all(r['gate'] == 'PASS' and r['status'] == '0' for r in sr+lr+gr),id_
        assert all(math.isfinite(float(r['minimum_concentration'])) and
                   float(r['minimum_concentration']) >= 0 for r in sr+lr+gr),id_
        level = 2 if id_ == 'economy' else 3
        assert all(float(r['overall_digits']) >= level and float(r['worst_digits']) >= level-1 for r in sr+lr+gr),id_
        assert all(r['cells'] == '10000' and r['accuracy_scope'] == 'first 1000 cells only'
                   for r in lr),id_
        assert len({r['binary_sha256'] for r in sr+lr}) == 1,id_
        times = {}
        for name in SCENARIOS:
            samples = select(large,id_,name)
            assert len(samples) == 5 and {r['repeat'] for r in samples} == {'1','2','3','4','5'}
            times[name] = statistics.median(float(r['elapsed_ms']) for r in samples)
        summary[id_] = dict(overall=min(float(r['overall_digits']) for r in sr),
            worst=min(float(r['worst_digits']) for r in sr),times=times,
            geometric_ms=math.exp(sum(math.log(t) for t in times.values())/len(times)),
            maximum_clip=max(float(r['maximum_clip']) for r in sr),
            largest_mean_projection=max(float(r['clipped_amount'])/int(r['cells']) for r in sr),
            level=level)
    lines = ['# Results','',
        'Tier 3 measurements on the supplied 10-core Apple M5, macOS 27.0.1. '
        'Native C uses clang 21 or GCC 15.2, ISO C11/OpenMP, and '
        '`OPTFLAGS="-O3 -DNDEBUG"`. MICM is compiled with clang++ and identical '
        'optimization flags. No target, fast-math, LTO, or profile flags are used.','',
        '## Declared operating points','',
        'All four points complete every disclosed scenario and all three independent '
        '48-hour cases with nonnegative states and pass their indicated accuracy level. '
        '**Balanced is recommended** for its accuracy margin. Full commands are in '
        '[OPERATING_POINTS.md](OPERATING_POINTS.md).','',
        '| Point | Target level | Lowest overall digits, 9 cases | Lowest species digits | Geometric mean of six 10,000-cell median times (ms) |',
        '| --- | ---: | ---: | ---: | ---: |']
    for point in points:
        s = summary[point['id']]
        lines.append(f"| {point['id']} | {s['level']} | {s['overall']:.3f} | {s['worst']:.3f} | {s['geometric_ms']:,.1f} |")
    lines += ['', 'The geometric mean is a descriptive summary, not a claim about the '
        'experimenter’s undisclosed ranking aggregation. Accuracy uses full disclosed '
        '1,000-cell cases and full 96-cell stress cases. The 10,000-cell runs are scored '
        'on the disclosed 1,000-cell prefix only; the C operating points’ full trajectories are checked '
        'for finite/nonnegative states by the native runner.','',
        '## Disclosed accuracy','',
        'Each entry is **overall / worst-species digits** from the supplied scorer.','',
        '| Scenario | '+ ' | '.join(p['id'] for p in points)+' |',
        '| --- | '+ ' | '.join('---:' for _ in points)+' |']
    for name in SCENARIOS:
        values = []
        for p in points:
            r, = select(small,p['id'],name)
            values.append(f"{float(r['overall_digits']):.3f} / {float(r['worst_digits']):.3f}")
        lines.append('| '+name+' | '+' | '.join(values)+' |')
    lines += ['', '## Native timing at 10,000 cells','',
        'Milliseconds, **median of five sequential runs**, ten requested threads. '
        'No build, native validation suite, or other solver benchmark from this '
        'repository was run concurrently. The C timer includes '
        'all steps and photolysis/rate updates, and excludes input, setup, diagnostics, '
        'snapshot copying and output. MICM uses its unchanged supplied timer, which '
        'also includes in-memory snapshot copying. These are local measurements; machine-wide '
        'idleness is not independently certified.','',
        '| Scenario | MICM default 3-stage | '+' | '.join(p['id'] for p in points)+' |',
        '| --- | ---: | '+' | '.join('---:' for _ in points)+' |']
    baseline_times = {}
    for name in SCENARIOS:
        br = select(baseline,'micm-default',name)
        assert len(br) == 5
        failures = sum(r['status'] != '0' or r['gate'] != 'PASS' for r in br)
        if failures:
            value = f'failed {failures}/5'; baseline_times[name] = None
        else:
            baseline_times[name] = statistics.median(float(r['elapsed_ms']) for r in br)
            value = f'{baseline_times[name]:,.1f}'
        lines.append('| '+name+' | '+value+' | '+' | '.join(f"{summary[p['id']]['times'][name]:,.1f}" for p in points)+' |')
    lines += ['', 'C uses groups of 32 cells with dynamic group scheduling. The MICM '
        'baseline uses its default one group per requested thread: 1,000 cells '
        'per group here, and 100 in the smaller comparison below. The implementations, '
        'grouping and adaptive policies therefore differ; this is a comparison '
        'against the supplied baseline, not an attribution of all time differences '
        'to the integration formula.','',
        'The default baseline fails on some long scenarios and does not meet '
        'the two-digit target on the 5-minute diurnal case. Consequently, its timings '
        'are not an equal-accuracy speed comparison, and no speedup is assigned to '
        'failed runs. For a robust three-stage baseline, scan 05 uses rtol `1e-5`, '
        'atol `1e-18`, and 100,000 maximum steps.','',
        '| 1,000-cell scenario | MICM default digits / worst | Tight MICM 3-stage ms | Tight overall / worst |',
        '| --- | ---: | ---: | ---: |']
    measured_baseline = read('05-baselines.csv')
    for name in SCENARIOS:
        normal, = select(measured_baseline,'micm-default',name)
        tight, = select(measured_baseline,'micm-tight3',name)
        value = f"{normal['overall_digits']} / {normal['worst_digits']}" if normal['status']=='0' else 'failed'
        lines.append(f"| {name} | {value} | {float(tight['elapsed_ms']):,.1f} | {tight['overall_digits']} / {tight['worst_digits']} |")
    lines += ['', 'These tight-baseline times are single native observations and obtain '
        'more accuracy than the declared targets. Their purpose is to show a working '
        'MICM three-stage comparison alongside the default failures.','',
        '## Independent 48-hour cases','',
        'The new cases use 96 cells each: winter night at 60° N; summer at 65° N '
        'with NOx multiplied by 20; and tropical conditions at 10° S with NOx '
        'multiplied by 0.03. Seeds, forcing and reference hashes are in '
        '[validation/stress-references.json](validation/stress-references.json). '
        'Two tight six-stage references agree to 10.193–10.412 overall digits '
        'and at least 9.494 digits for every species.','',
        '| Scenario | '+' | '.join(p['id'] for p in points)+' |',
        '| --- | '+' | '.join('---:' for _ in points)+' |']
    for name in ['winter-night-48h','summer-polluted-48h','tropical-clean-48h']:
        values = []
        for p in points:
            r, = select(stress,p['id'],name)
            values.append(f"{r['overall_digits']} / {r['worst_digits']}")
        lines.append('| '+name+' | '+' | '.join(values)+' |')
    lines += ['', '## Projection and implementation checks','',
        '| Point | Largest individual projection (mol m-3) | Largest trajectory-total projection per cell, averaged across cells (mol m-3) |',
        '| --- | ---: | ---: |']
    for p in points:
        s = summary[p['id']]
        lines.append(f"| {p['id']} | {s['maximum_clip']:.3e} | {s['largest_mean_projection']:.3e} |")
    lines += ['', 'Projection totals sum concentration increments over species and '
        'accepted steps; they are not elemental mass balances. All indirect effects '
        'remain included in the scored trajectories.','',
        'Both compilers pass **all 680 named Tier 1 cases**, the native extension '
        'tests, 147 C/prototype comparisons with identical operation counts, seven '
        'stiff analytical cases, linear-chain order checks, failure/rejection tests, the '
        'projected SDIRK error regression, and '
        'the scenario I/O/forcing/thread contract. Every declared point also passes '
        'all nine accuracy cases with GCC.','']
    for compiler in ['clang','gcc']:
        for name in ['small','large']:
            lines.append(f"- {compiler}, {name} Tier 1 gate: {gates[compiler]['gates'][name]['comparison']}")
    lines += ['', 'Authoritative records: [compiler gates](validation/tier3-gates.json), '
        '[clang cases](validation/tier3-clang-cases.txt), '
        '[GCC cases](validation/tier3-gcc-cases.txt), and '
        '[scans/](scans/). The legacy implementation, configuration and all named '
        'cases are retained. Historical Tier 1 measurements are archived in '
        '[validation/TIER1_RESULTS.md](validation/TIER1_RESULTS.md). Linux was not '
        'available for execution; portability is supported by source inspection '
        'and the two supplied compilers on this host.','',
        '## Work–precision plots','',
        'The x coordinate includes the species constraint: '
        '`min(overall digits, worst-species digits + 1)`. Vertical lines mark '
        'levels 2 and 3. Exploratory scans show individual observations; final '
        'points and the MICM default baseline show five-run medians and full timing '
        'ranges. Failed baseline cases are labeled without assigning a time. The '
        'six exploratory observations known to overlap reference compilation and '
        'the six superseded SDIRK observations are omitted from the plot and remain '
        'in the raw scan logs. Scan 11 supplies the corrected SDIRK points.','',
        '![Exploratory work–precision](figures/work-precision-scans.png)','',
        '![Operating-point work–precision](figures/work-precision-operating-points.png)','',
        'Vector exports: [scan SVG](figures/work-precision-scans.svg), '
        '[operating-point SVG](figures/work-precision-operating-points.svg). '
        'Regenerate with `.conda/bin/python tools/plot_scans.py`.','']
    (ROOT/'RESULTS.md').write_text('\n'.join(lines))
    (ROOT/'validation/operating-summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__ == '__main__':
    main()
