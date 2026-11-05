#!/usr/bin/env python3
"""Fill the paper's timing results from the session's timings.csv files: the
timing subsection and table, the abstract, the key points, the summary, the
discussion, the conclusions, the scans section, and the title."""
import csv, re, statistics, sys
from pathlib import Path

OUT = Path.home() / "NSF-NCAR/tier3-paper-work/out"
P = Path.home() / "DavidFillmore/ODE/paper_efficient_ode"
ORDER = [("disclosed", s) for s in ("tier1-30s", "diurnal-30min", "diurnal-10min", "diurnal-5min", "polluted-10min", "clean-10min")] + \
        [("heldout", s) for s in ("summer-15min", "tropics-20min", "winter60-10min", "polluted-48h-30min", "clean-48h-5min")]
med = {}
for set_, s in ORDER:
    byp = {}
    for r in csv.DictReader(open(OUT / set_ / s / "timings.csv")):
        assert r["gate"] == "PASS", (s, r["program"], r["gate"])
        byp.setdefault(r["program"], []).append(float(r["elapsed_ms"]) / 1000.0)
    for p, v in byp.items():
        assert len(v) == 3, (s, p, len(v))
        med[(s, p)] = statistics.median(v)
tot = lambda p: sum(med[(s, p)] for _, s in ORDER)
ratio = lambda a, b: tot(a) / tot(b)
per = lambda a, b: [med[(s, a)] / med[(s, b)] for _, s in ORDER]
v = dict(
    micmC2=tot("micm/C2"), micmC3=tot("micm/C3"), loopC2=tot("loop/C2"), loopC3=tot("loop/C3"), retC2=tot("retaining/C2"), retC3=tot("retaining/C3"),
    rC2=ratio("loop/C2", "retaining/C2"), rC3=ratio("loop/C3", "retaining/C3"), mC2=ratio("micm/C2", "retaining/C2"), mC3=ratio("micm/C3", "retaining/C3"),
    kC2=ratio("micm/C2", "loop/C2"), kC3=ratio("micm/C3", "loop/C3"),
    minC3=min(per("loop/C3", "retaining/C3")), maxC3=max(per("loop/C3", "retaining/C3")),
    minC2=min(per("loop/C2", "retaining/C2")), maxC2=max(per("loop/C2", "retaining/C2")))
for k, x in v.items():
    print(f"{k} = {x:.2f}")
if "--dry" in sys.argv:
    sys.exit()
f1 = lambda x: f"{x:.1f}"
def edit(rel, pairs):
    p = P / rel; s = p.read_text()
    for a, b in pairs:
        n = len(re.findall(a, s, flags=re.S)); assert n == 1, (rel, a[:50], n)
        s = re.sub(a, lambda m: b, s, flags=re.S)
    p.write_text(s)
title = "A Fast ODE Solver for Atmospheric Chemistry: ROS3 with H211b Control and Retained Steps"
edit("main.tex", [
 (r"\\title\{An Efficient ODE Solver for Atmospheric Chemistry: ROS3 with H211b Control and Retained Steps\}", "\\title{" + title + "}"),
 (r"pdftitle=\{An Efficient ODE Solver for Atmospheric Chemistry: ROS3 with H211b Control and Retained Steps\}", "pdftitle={" + title + "}"),
 (r"\\item At 0\.1\\% pooled and 1\\% per-species RMS relative error, the solver needs one third of a restarting solver's arithmetic on MOZART-TS1\.",
  "\\item At 0.1\\% pooled and 1\\% per-species RMS relative error, the solver needs a third of a restarting solver's arithmetic and time on MOZART-TS1."),
 (r"The step sequences of Figure~\\ref\{fig:trace\} were\nrecorded with the archived code on the archived cells of Diurnal 10~min\.",
  "The step sequences of Figure~\\ref{fig:trace} were\nrecorded with the archived code on the archived cells of Diurnal 10~min,\nand the timing session of Section~\\ref{sec:timing} is recorded in the\nresults directory of the TurboMICM repository \\cite{turbomicm}."),
])
edit("sections/00_abstract.tex", [
 (r"These are reductions in arithmetic\s+work; elapsed-time improvements remain unmeasured\.",
  f"In a compiled library the elapsed time against the restarting loop falls by factors of {f1(v['rC2'])} and {f1(v['rC3'])}, matching the arithmetic."),
 (r"Operator-split atmospheric models integrate stiff chemistry over many\s+short host intervals, and solvers that restart each interval from a small\s+step repeat initialization work that the state does not require\.",
  "Solvers that restart every host interval of an operator-split atmospheric\nmodel from a small step repeat initialization work that the state does not\nrequire."),
])
edit("sections/00_plain_language_summary.tex", [
 (r"The saving in computer\s+time remains to be measured\.", "In a compiled implementation the saving in computer time matches the\nsaving in arithmetic."),
])
edit("sections/06_scans.tex", [
 (r"The scans timed the compiled\s+library; this paper counts arithmetic in the Python implementation on\s+all eleven scenarios, and the elapsed-time comparison of the final\s+solver under this protocol remains to be made\.",
  "The scans timed the compiled library on the six disclosed scenarios;\nthis paper counts arithmetic in the Python implementation on all eleven,\nand Section~\\ref{sec:timing} times the final solver in the same library."),
])
edit("sections/09_discussion.tex", [
 (r"Operation counts,\s+even weighted, do not establish elapsed-time savings; a compiled\s+implementation is needed for that\.",
  "The elapsed-time ratios come from one compiled implementation on one\nlaptop and from box-model scenarios; a coupled model with transport\nremains untested."),
 (r"Execution-time gains, transfer to other mechanisms,\s+and error bounds for the transient after a rate change remain open\.",
  f"In a compiled library the elapsed time falls by factors of {f1(v['rC2'])} and\n{f1(v['rC3'])}. Transfer to other mechanisms, a coupled model, and error\nbounds for the transient after a rate change remain open."),
])
timing = f"""
\\subsection{{Elapsed time in a compiled library}}
\\label{{sec:timing}}

To measure time we ran the retaining solver, the restarting loop, and
MICM's own loop in the compiled library of Section~\\ref{{sec:scans}}
\\cite{{reaper2026}}, whose loops reproduce this paper's counts exactly, on
the 10,000-cell versions of the eleven scenarios at the tolerances of
Table~\\ref{{tab:levels}}, with 32-cell groups and ten threads on a laptop
(Apple M5 Pro), each time the median of three rotating rounds.
Table~\\ref{{tab:timing}} sums the medians over the scenarios. The retaining
solver's time is {f1(v['rC2'])} and {f1(v['rC3'])} times below the restarting loop's at
criteria 2 and 3, against arithmetic ratios of 3.6 and 3.0, and the ratio
follows the attempt ratio scenario by scenario, from {f1(v['minC3'])} on
Diurnal 30~min to {f1(v['maxC3'])} on Benchmark 30~s at criterion 3: with
identical kernels on both sides, time tracks the number of attempts.
Against MICM's own loop the factors are {f1(v['mC2'])} and {f1(v['mC3'])}, of which the
generated kernels account for {f1(v['kC2'])} and {f1(v['kC3'])}, the ratio between MICM's
loop and the same loop in the library.

\\begin{{table}}[t]
\\centering
\\caption{{Elapsed time summed over the eleven scenarios at 10,000 cells and
ten threads (medians of three rounds, in s), at the criterion-2 and
criterion-3 tolerances of Table~\\ref{{tab:levels}}, and the ratio to the
retaining solver.}}
\\label{{tab:timing}}
\\small
\\begin{{tabular}}{{lrrrr}}
\\toprule
& \\multicolumn{{2}}{{c}}{{Criterion 2}} & \\multicolumn{{2}}{{c}}{{Criterion 3}} \\\\
Solver & Time (s) & Ratio & Time (s) & Ratio \\\\
\\midrule
MICM's loop & {v['micmC2']:.1f} & {f1(v['mC2'])} & {v['micmC3']:.1f} & {f1(v['mC3'])} \\\\
Restarting loop, library kernels & {v['loopC2']:.1f} & {f1(v['rC2'])} & {v['loopC3']:.1f} & {f1(v['rC3'])} \\\\
Retaining solver, library kernels & {v['retC2']:.1f} & 1 & {v['retC3']:.1f} & 1 \\\\
\\bottomrule
\\end{{tabular}}
\\end{{table}}
"""
p = P / "sections/08_results.tex"; s = p.read_text()
assert "sec:timing" not in s
s = s.rstrip("\n") + "\n" + timing
p.write_text(s)
p = P / "tools/check_build.py"; s = p.read_text()
s = s.replace('TABLES = ["1", "2", "3"]', 'TABLES = ["1", "2", "3", "4"]').replace("three tables in order", "four tables in order")
p.write_text(s)
for rel in ("README.md",):
    p = P / rel; s = p.read_text()
    s = s.replace("# An Efficient ODE Solver for Atmospheric Chemistry: ROS3 with H211b Control and Retained Steps", "# " + title)
    p.write_text(s)
print("paper sources updated")
