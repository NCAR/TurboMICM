"""Run reproducible C experiments, scoring with the supplied score.py.

Configurations are JSON objects: id, method, options. Every attempted run is
appended immediately, including failures. Subsets preserve the first cells and
their matching reference rows; full experiments use the original input paths.
"""
import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = ROOT.parent / "scenarios"
FIELDS = ["scan", "config", "method", "options", "scenario", "cells", "repeat", "threads",
          "status", "overall_digits", "worst_species", "worst_digits", "gate", "elapsed_ms",
          "internal_steps", "accepted", "rejected", "function_calls", "jacobian_updates",
          "decompositions", "solves", "jacobian_reuses", "factor_reuses",
          "positivity_rejections", "nonfinite_rejections", "clipped_values", "clipped_mass",
          "completed_steps", "score_json", "run_log"]


def subset(name, cells):
    source = SCENARIOS / name
    if not cells:
        return source, sum(1 for _ in (source / "cells.csv").open()) - 1
    dest = ROOT / "build" / "scan-inputs" / f"{name}-{cells}"
    dest.mkdir(parents=True, exist_ok=True)
    for filename in ["cells.csv", "reference_final.csv", "reference_snapshots.csv"]:
        path = source / filename
        if not path.exists():
            continue
        with path.open() as f, (dest / filename).open("w", newline="") as out:
            reader, writer = csv.reader(f), csv.writer(out)
            header = next(reader); writer.writerow(header)
            for i, row in enumerate(reader):
                if (int(row[1]) if header[0] == "step" else i) < cells:
                    writer.writerow(row)
    for filename in ["scenario.txt", "photolysis.csv"]:
        (dest / filename).write_text((source / filename).read_text())
    return dest, cells


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configs", type=Path)
    parser.add_argument("--scan", required=True)
    parser.add_argument("--cells", type=int, default=0)
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--scenarios", nargs="*", default=["tier1-30s", "diurnal-30min", "diurnal-10min",
                        "diurnal-5min", "polluted-10min", "clean-10min"])
    parser.add_argument("--timeout", type=float, default=1800)
    a = parser.parse_args()
    configs = json.loads(a.configs.read_text())
    output = ROOT / "scans" / f"{a.scan}.csv"
    if output.exists():
        raise SystemExit(f"Refusing to overwrite {output}")
    with output.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS); writer.writeheader(); out.flush()
        for config in configs:
            for name in a.scenarios:
                scenario, cells = subset(name, a.cells)
                for repeat in range(a.repeats):
                    work = ROOT / "build" / "scan-runs" / a.scan / config["id"] / name / str(repeat)
                    work.mkdir(parents=True, exist_ok=True)
                    command = [str(ROOT / "bin" / "ts1_scenario"), "--threads", str(a.threads),
                               "--method", config["method"]]
                    for key, value in config.get("options", {}).items():
                        command += ["--option", f"{key}={value}"]
                    command += [str(scenario / "scenario.txt")]
                    print(a.scan, config["id"], name, repeat, flush=True)
                    try:
                        run = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=a.timeout)
                        log = run.stdout + run.stderr
                        status = "ok" if run.returncode == 0 else f"exit-{run.returncode}"
                    except subprocess.TimeoutExpired as e:
                        log = str(e); status = "timeout"
                    (work / "run.txt").write_text(" ".join(command) + "\n" + log)
                    values = dict(re.findall(r"(\w+)=([^\s]+)", log))
                    row = {k: values[k] for k in FIELDS if k in values}
                    row.update(scan=a.scan, config=config["id"], method=config["method"],
                               options=json.dumps(config.get("options", {}), sort_keys=True),
                               scenario=name, cells=cells, repeat=repeat, threads=a.threads, status=status,
                               run_log=str((work / "run.txt").relative_to(ROOT)))
                    if status == "ok":
                        snapshots = "snapshot_every" in (scenario / "scenario.txt").read_text()
                        suffix = "snapshots.csv" if snapshots else "final.csv"
                        score = subprocess.run([sys.executable, str(SCENARIOS / "score.py"),
                                                str(work / ("ts1_" + suffix)),
                                                str(scenario / ("reference_" + suffix)), "--json"],
                                               capture_output=True, text=True, check=False)
                        (work / "score.json").write_text(score.stdout + score.stderr)
                        try:
                            report = json.loads(score.stdout)
                            worst = report["worst_species"][0]
                            row.update(overall_digits=report["overall_digits"], worst_species=worst[0],
                                       worst_digits=worst[1], gate=report["gate"], score_json=score.stdout.strip())
                        except (ValueError, KeyError, IndexError):
                            row["status"] = "score-failed"
                    writer.writerow(row); out.flush()
                    print({k: row.get(k) for k in ["status", "overall_digits", "worst_digits", "elapsed_ms"]}, flush=True)


if __name__ == "__main__":
    main()
