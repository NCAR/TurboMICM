"""Run reproducible C experiments, scoring with the supplied score.py.

Configurations are JSON objects: id, method, options. Every attempted run is
appended immediately, including failures. Subsets preserve the first cells and
their matching reference rows; full experiments use the original input paths.
"""
import argparse
import csv
import json
import hashlib
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = ROOT.parent / "scenarios"
FIELDS = ["scan", "config", "method", "options", "scenario", "cells", "repeat", "threads",
          "utc", "binary_sha256", "compiler", "output_sha256", "score_reused",
          "status", "overall_digits", "worst_species", "worst_digits", "gate", "level", "level_result", "elapsed_ms",
          "internal_steps", "accepted", "rejected", "function_calls", "jacobian_updates",
          "decompositions", "solves", "jacobian_reuses", "factor_reuses",
          "positivity_rejections", "nonfinite_rejections", "clipped_values", "clipped_mass",
          "nonlinear_iterations", "nonlinear_rejections",
          "extrap_order2", "extrap_order3", "extrap_order4",
          "rate_updates", "rate_cache",
          "completed_steps", "score_json", "run_log"]


def subset(name, cells, scenario_root=SCENARIOS):
    source = scenario_root / name
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
    parser.add_argument("--scenario-root", type=Path, default=SCENARIOS)
    parser.add_argument("--discard-repeats", action="store_true", help="Remove state CSVs after scoring repetitions after the first")
    parser.add_argument("--reuse-identical-scores", action="store_true", help="Reuse the scorer result only for bit-identical output and the same reference/level")
    a = parser.parse_args()
    configs = json.loads(a.configs.read_text())
    compiler = (ROOT / "build" / "compiler.conf").read_text().strip().replace("\n", " | ")
    score_cache = {}
    output = ROOT / "scans" / f"{a.scan}.csv"
    if output.exists():
        raise SystemExit(f"Refusing to overwrite {output}")
    with output.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS); writer.writeheader(); out.flush()
        for config in configs:
            for name in a.scenarios:
                scenario, cells = subset(name, a.cells, a.scenario_root.resolve())
                for repeat in range(a.repeats):
                    work = ROOT / "build" / "scan-runs" / a.scan / config["id"] / name / str(repeat)
                    work.mkdir(parents=True, exist_ok=True)
                    binary = ROOT / "bin" / ("micm_scenario" if config.get("backend") == "micm" else "ts1_scenario")
                    binary_hash = hashlib.sha256(binary.read_bytes()).hexdigest()
                    command = [str(binary), "--threads", str(a.threads),
                               "--method", config["method"]]
                    for key, value in config.get("options", {}).items():
                        command += [f"--{key}", str(value)] if config.get("backend") == "micm" else ["--option", f"{key}={value}"]
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
                    actual_compiler = ('clang++ -std=c++20 ' + (ROOT/'build/reference-compiler.conf').read_text().splitlines()[1]
                                       if config.get('backend') == 'micm' else compiler)
                    row.update(utc=datetime.now(timezone.utc).isoformat(), binary_sha256=binary_hash, compiler=actual_compiler)
                    if status == "ok":
                        settings = {}
                        for line in (scenario / "scenario.txt").read_text().splitlines():
                            line = line.split('#',1)[0]
                            if '=' in line:
                                key,value=line.split('=',1); settings[key.strip()]=value.strip()
                        snapshots = int(settings.get("snapshot_every",0)) > 0
                        suffix = "snapshots.csv" if snapshots else "final.csv"
                        state_file=work/("ts1_"+suffix)
                        with state_file.open('rb') as state:
                            output_hash=hashlib.file_digest(state,'sha256').hexdigest()
                        cache_key=(str(scenario),suffix,config.get('level'),output_hash)
                        cached=score_cache.get(cache_key) if a.reuse_identical_scores else None
                        score = subprocess.run([sys.executable, str(SCENARIOS / "score.py"),
                                                str(work / ("ts1_" + suffix)),
                                                str(scenario / ("reference_" + suffix)), "--json"] +
                                               (["--level", str(config["level"])] if "level" in config else []),
                                               capture_output=True, text=True, check=False) if cached is None else None
                        score_text=score.stdout if score is not None else cached
                        row.update(output_sha256=output_hash,score_reused=int(cached is not None))
                        (work / "score.json").write_text(score_text + (score.stderr if score is not None else ''))
                        try:
                            report = json.loads(score_text)
                            worst = report["worst_species"][0]
                            row.update(overall_digits=report["overall_digits"], worst_species=worst[0],
                                       worst_digits=worst[1], gate=report["gate"], score_json=score_text.strip())
                            row.update(level=report.get("level", ""), level_result=report.get("level_result", ""))
                            score_cache[cache_key]=score_text
                        except (ValueError, KeyError, IndexError):
                            row["status"] = "score-failed"
                    writer.writerow(row); out.flush()
                    print({k: row.get(k) for k in ["status", "overall_digits", "worst_digits", "elapsed_ms"]}, flush=True)
                    if a.discard_repeats and repeat:
                        for filename in ["ts1_final.csv", "ts1_snapshots.csv"]:
                            (work/filename).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
