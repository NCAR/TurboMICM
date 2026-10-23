#!/usr/bin/env python3
"""Check generated chemistry schedules against generic evaluation and solves."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MICM_CONFIGS = ROOT.parent / "micm" / "test" / "tutorial" / "configs"
CC = os.environ.get("CC", "clang")
OPTFLAGS = os.environ.get("OPTFLAGS", "-O3 -DNDEBUG").split()
OMPFLAGS = os.environ.get("OMPFLAGS", "").split()
OMPLIBS = os.environ.get("OMPLIBS", "").split()
LIBSRC = ["mechanism.c", "solver.c", "state.c", "constraints.c", "matrix.c",
          "sparse.c", "lu.c", "linear.c", "system.c", "version.c"]


def run(args, **kwargs):
    subprocess.run([str(arg) for arg in args], check=True, **kwargs)


with tempfile.TemporaryDirectory() as temporary:
    temporary = Path(temporary)
    input_csv = ROOT.parent / "reference" / "ts1_cells_100_seed1.csv"
    generic = ROOT / "bin" / "ts1_bench_cpu"
    specialized = ROOT / "bin" / "ts1_bench_specialized"
    for name, executable in (("generic", generic), ("specialized", specialized)):
        work = temporary / name
        work.mkdir()
        run([executable, "--threads", "1", "100", "100", "30", input_csv],
            cwd=work, stdout=subprocess.DEVNULL)
    run([sys.executable, ROOT.parent / "reference" / "compare_final.py",
         temporary / "specialized" / "ts1_final.csv",
         temporary / "generic" / "ts1_final.csv", "--rtol", "1e-10",
         "--atol", "1e-30"])

    harness = ROOT / "tests" / "specialized_compare.c"
    for directory in ("rate_constants_no_user_defined",
                      "rate_constants_user_defined"):
        config_dir = MICM_CONFIGS / directory
        config = config_dir / ("config.yaml" if directory.endswith("no_user_defined")
                               else "config.json")
        normalized = temporary / f"{directory}.json"
        generated = temporary / f"{directory}.c"
        executable = temporary / directory
        run([ROOT / ".conda" / "bin" / "python", ROOT / "tools" / "read_config.py",
             config, normalized])
        run([ROOT / ".conda" / "bin" / "python", ROOT / "tools" / "generate_c.py",
             normalized, generated])
        run([CC, "-std=c11", *OPTFLAGS, *OMPFLAGS, "-DCHEM_SPECIALIZED",
             "-I", ROOT / "include", "-I", ROOT / "src",
             *[ROOT / "src" / source for source in LIBSRC], generated,
             harness, "-lm", *OMPLIBS, "-o", executable])
        run([executable])

print("specialized solver and generated schedules agree with generic evaluation")
