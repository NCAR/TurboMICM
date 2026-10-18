#!/usr/bin/env python3
"""Build the results report from current run records; reject incomplete evidence."""
import csv
import itertools
import json
from pathlib import Path
import re
from statistics import median
from scan import ROOT,NAMES

def read(name):
    return list(csv.DictReader((ROOT/'scans'/name).open()))

def matched(config,point):
    c=json.loads(config)
    return c=={'method':point['method'],**point['options']}

def main():
    points=json.loads((ROOT/'operating_points.json').read_text())
    ids=[p['id'] for p in points]
    rows=read('operating-points.csv')
    assert len(rows)==150
    assert {(r['point'],r['scenario'],int(r['repeat'])) for r in rows}==set(itertools.product(ids,NAMES,range(5)))
    assert all(r['returncode']=='0' and r['gate']=='PASS' and r['level_result']=='PASS' for r in rows)
    gcc=read('operating-points-gcc.csv')
    assert len(gcc)==30 and all(r['gate']=='PASS' and r['level_result']=='PASS' for r in gcc)
    assert {(r['point'],r['scenario']) for r in gcc}==set(itertools.product(ids,NAMES))
    stress=[r for r in read('stress.csv') if r['repeat']=='1']
    assert len(stress)==15 and all(r['gate']=='PASS' and r['level_result']=='PASS' for r in stress)
    scale=read('scale-points.csv')
    assert len(scale)==5 and all(r['gate']=='PASS' and r['level_result']=='PASS' for r in scale)
    for cc in ['clang','gcc']:
        text=(ROOT/f'gate-{cc}-final.log').read_text()
        for expected in ['cases=680 PASS=675 FAIL=0 N/A=5','48 C/prototype comparisons PASS',
                         '64 forcing-cache checks PASS','Analytic order and invariant checks PASS',
                         'Independent MICM scenario driver parity PASS']:
            assert expected in text,(cc,expected)
    assert 'cases=680 PASS=675 FAIL=0 N/A=5' in (ROOT/'gate-clang-restored.log').read_text()
    gates={}
    for cc,size,path in [('clang',100,'tier1-clang-final.log'),('clang',10000,'tier1-clang-10000.log'),
                         ('gcc',100,'tier1-gcc-100.log'),('gcc',10000,'tier1-gcc-10000.log')]:
        text=(ROOT/path).read_text()
        assert f'PASS: {size*211} values in {size} cells' in text
        gates[cc,size]=re.search(r'largest relative difference ([^ ]+)',text).group(1)
    sanitizer=(ROOT/'sanitizer-new.log').read_text()
    assert 'Analytic order and invariant checks PASS' in sanitizer
    assert 'invalid options PASS' in sanitizer and 'ERROR: AddressSanitizer' not in sanitizer and 'runtime error:' not in sanitizer
    linked=(ROOT/'linked-libraries.log').read_text()
    assert 'libomp' in linked and all(s not in linked.lower() for s in ['python','accelerate','dispatch','asan','ubsan'])
    assert '680 IDs' in (ROOT/'case-manifest.log').read_text()
    baseline=read('micm-tuned.csv')
    assert len(baseline)==18 and all(r['level_result']=='PASS' for r in baseline)
    baseline2=read('micm-level2.csv')
    assert len(baseline2)==18 and all(r['gate']=='PASS' and float(r['overall_digits'])>2.001 and float(r['worst_digits'])>1.001 for r in baseline2)
    default=read('micm-default.csv')
    lines=['# Tier 3 results','',
      'Five operating points pass their declared accuracy levels on every disclosed',
      'scenario with both compilers, and on three independent 48-hour cases. The',
      'recommendation uses RODAS3.1 with rtol=3e-4 and atol=1e-18 mol/m3. It retains',
      'a current Jacobian at every accepted step and solves the complete mechanism.',
      'See OPERATING_POINTS.md for exact commands, ACCURACY.md for the error budget,',
      'and SCANS.md / scans/ for controlled experiments and raw records.','',
      '## Accuracy','',
      'Digits are from the unmodified supplied score.py. Each entry below is the',
      'minimum over scenarios, not the score of an average state. All state values',
      'pass the finite/nonnegative gate. The last columns use independent references',
      'for high-latitude summer, winter and equatorial 48-hour conditions.','',
      '| Point | Declared level | Disclosed overall | Disclosed worst species | 48-hour overall | 48-hour worst species |',
      '|---|---:|---:|---:|---:|---:|']
    for p in points:
        rs=[r for r in rows if r['point']==p['id']]
        ss=[r for r in stress if matched(r['configuration'],p)]
        assert len(ss)==3
        assert {r['scenario'] for r in ss}=={'stress-summer','stress-winter','stress-equator'}
        lines.append(f"| {p['id']} | {p['level']} | {min(float(r['overall_digits']) for r in rs):.3f} | {min(float(r['worst_digits']) for r in rs):.3f} | {min(float(r['overall_digits']) for r in ss):.3f} | {min(float(r['worst_digits']) for r in ss):.3f} |")
    lines+=['','The additional references agree under tolerance refinement by 10.076–10.349',
      'digits overall and at least 8.992 digits per species. Their settings, work',
      'and convergence are recorded in scans/stress-reference-*.csv and',
      'scans/stress_conditions.json. This tests additional conditions and horizon;',
      'it does not prove a uniform error bound on undisclosed inputs.','',
      '## C timings on disclosed 1000-cell scenarios','',
      'Apple clang 21, -O3 -DNDEBUG, ten threads, groups of 32, guided scheduling',
      'with minimum chunk 1. Numbers are median milliseconds from five sequential',
      'runs in a deterministic shuffled order. No other CPU-heavy task processes',
      'ran during these measurements. Input, setup, snapshot copying and CSV output',
      'are excluded; rate updates, all chemistry steps and returned-state checks',
      'are included. Python provides no timing measurements.','',
      '| Scenario | L2 | L2_margin | L3 | recommended | L4 | MICM3 tuned for level 3 |',
      '|---|---:|---:|---:|---:|---:|---:|']
    times={}
    for name in NAMES:
        vals=[]
        for p in points:
            rs=[r for r in rows if r['point']==p['id'] and r['scenario']==name]
            assert len(rs)==5
            for field in ['overall_digits','worst_digits','internal_steps','function_calls','jacobian_updates','decompositions','solves']:
                assert len({r[field] for r in rs})==1,(name,p['id'],field)
            t=median(float(r['elapsed_ms']) for r in rs);times[p['id'],name]=t;vals.append(f'{t:.1f}')
        b=median(float(r['elapsed_ms']) for r in baseline if r['scenario']==name)
        times['baseline3',name]=b
        lines.append('| '+name+' | '+' | '.join(vals+[f'{b:.1f}'])+' |')
    ratios3=[times['baseline3',n]/times['L3',n] for n in NAMES]
    ratiosrec=[times['baseline3',n]/times['recommended',n] for n in NAMES]
    ratios2=[]
    for name in NAMES:
        b=median(float(r['elapsed_ms']) for r in baseline2 if r['scenario']==name)
        ratios2.append(b/times['L2',name])
    lines+=['',
      'The comparison baseline is the independently compiled, unmodified supplied',
      'C++ MICM three-stage scenario runner, with the same OPTFLAGS and groups of',
      '32. Its level-3 setting is rtol=1e-3, atol=1e-18; its level-2 setting is',
      'rtol=1e-2, atol=1e-18. Each baseline time is a median of three runs. Both',
      'settings meet the corresponding level on all six scenarios.','',
      f"Against these tested baselines, L2 is {min(ratios2):.2f}–{max(ratios2):.2f} times faster, L3 is {min(ratios3):.2f}–{max(ratios3):.2f} times faster,",
      f"and the recommendation is {min(ratiosrec):.2f}–{max(ratiosrec):.2f} times faster. These ratios combine the retained Tier 1",
      'C implementation improvements with the new numerical choices; they are not',
      'pure method-only speedups. SCANS.md and the C controls separate those effects.',
      'The supplied C++ runner also retains its own scheduling and snapshot buffering.',
      'These are local development measurements, not the evaluator’s 10000-cell',
      'median-of-five ranking.','',
      '## Original default baseline','',
      'The required default MICM3 method/tolerances remain available. Its atol=1e-3',
      'mol/m3 does not provide acceptable long-run accuracy. With groups of 32 the',
      'fresh C++ baseline gives:','',
      '| Scenario | Overall digits | Worst-species digits | C++ elapsed ms / status |',
      '|---|---:|---:|---|']
    for r in default:
        lines.append(f"| {r['scenario']} | {r['overall_digits'] or '—'} | {r['worst_digits'] or '—'} | {r['elapsed_ms'] if r['returncode']=='0' else 'solver failure'} |")
    lines+=['','Different group norms explain differences from the supplied README’s',
      '100-cell groups. Failed/inaccurate default runs are not used as an',
      'accuracy-qualified speed baseline. The tighter rtol=1e-4 baseline is also',
      'recorded in scans/micm-tight.csv.','',
      '## Scaling','',
      'Each point additionally runs all 10000 cells of the diurnal-10min scenario',
      'without nonfinite or negative output. The accuracy score covers only the',
      'disclosed first 1000 cells, whose input prefix is verified exactly; it is',
      '**not** a full 10000-cell reference comparison. The gate covers all 10000',
      'cells, and the C driver checks every returned external-step state.','',
      '| Point | C elapsed ms, one run | First 1000 overall digits | First 1000 worst species |',
      '|---|---:|---:|---:|']
    for p in points:
        r=next(r for r in scale if r['point']==p['id'])
        lines.append(f"| {p['id']} | {float(r['elapsed_ms']):.1f} | {r['overall_digits']} | {r['worst_digits']} |")
    lines+=['','Grouping/scheduling development trials are in scans/scale.csv; full point',
      'runs are in scans/scale-points.csv. All disclosed timings and operation',
      'counts are in scans/operating-points.csv. Counts sum over groups, not cells.',
      'Every one of the five repeats has identical scores and operation counts',
      'for the same point/scenario.','',
      '## Work–precision plots','',
      'Curves show C development scans; stars are the declared-point median-of-five',
      'timings. Failed runs and negative digit values remain in CSV records but are',
      'omitted from these plots. The species plot prevents pooled RMS from hiding',
      'a weak component.','',
      '![Overall work–precision](plots/work_precision.png)','',
      '![Worst-species work–precision](plots/work_precision_species.png)','',
      'Vector exports are plots/work_precision.pdf and plots/work_precision_species.pdf.',
      'Recreate both with `.conda/bin/python tools/plot_scans.py`.','',
      '## Retained Tier 1 gate and portability','',
      '| Check | clang | GCC 15.2 |','|---|---:|---:|',
      '| Exact case manifest | 680 IDs | 680 IDs |','| Runtime/numerical PASS | 675 | 675 |',
      '| FAIL | 0 | 0 |','| Existing C++-only N/A | 5 | 5 |',
      f"| TS1, 100 cells, maximum relative difference | {gates['clang',100]} | {gates['gcc',100]} |",
      f"| TS1, 10000 cells, maximum relative difference | {gates['clang',10000]} | {gates['gcc',10000]} |",
      '| Eight-method strict C/prototype comparisons | 48 PASS | 48 PASS |',
      '| Exact forcing-cache comparisons | 64 PASS | 64 PASS |',
      '| Linear/nonlinear analytic order and invariants | PASS | PASS |',
      '| All declared points, all six disclosed cases | PASS | PASS |','',
      'The TS1 comparisons inspect every cell/species at 1e-30 + 1e-10 relative;',
      'the 10000-cell reference is a fresh independent MICM run of the same full',
      'input, not a prefix or pooled-error check. GCC point records are in',
      'scans/operating-points-gcc.csv. The five exclusions and historical Tier 1',
      'timing/group-size studies are retained in TIER1_RESULTS.md.','',
      'Address/undefined-behavior sanitizer checks pass the new analytic methods,',
      'all eight method probes with reuse, and the scenario/partial-group driver',
      'tests. A source audit finds no forbidden intrinsics, target extensions or',
      'long double. The release C executable links only the system C/math runtime',
      'and OpenMP; Python and reference C++ are absent. No optimization flags are',
      'added beyond caller OPTFLAGS.','',
      'Aggressive H211b b=1 has documented roundoff sensitivity in independent',
      'prototype comparisons; none of the declared points uses it or approximate',
      'Jacobian reuse. No uniform error guarantee on unobserved conditions or',
      'performance claim for untested integrator families is made.','']
    historical=ROOT/'TIER1_RESULTS.md'
    if not historical.exists():
        historical.write_text('# Historical Tier 1 results\n\nRetained from the initial repository. Fresh Tier 3 gates and measurements are\nin RESULTS.md.\n\n'+(ROOT/'RESULTS.md').read_text())
    (ROOT/'RESULTS.md').write_text('\n'.join(lines))
    print('RESULTS.md generated from complete current validation and timing evidence')
if __name__=='__main__': main()
