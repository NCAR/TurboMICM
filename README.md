# TurboMICM

Optimized ODE solvers for the Model Independent Chemistry Module.

This repository tracks the MICM-to-C experiment, a controlled comparison of LLM
coding agents. Every agent gets the same prompt: port MICM's CPU solver to
portable C (C11, libm, and OpenMP) and make it fast without depending on one chip,
compiler, or operating system.

## Layout

`main` holds the experiment:

    prompt.v3     the prompt template; prompt and prompt.v2 are earlier versions
    prompts/      make_prompts.py and the 18 prompts it generated, one per run
    reference/    MICM's test-case list, the MICM TS1 reference driver, the TS1
                  cell generator, the accuracy gate, and reference data
                  (see reference/README.md)
    admin/        make_case_list.py, setup_runs.py, and runs.md, the run log

Each run's solution is on its own branch, `runs/<model>/<effort>`, holding the
agent's repository history unchanged, so commit hashes match the run log:

    runs/sol/high    runs/sol/xhigh    runs/sol/max    ...

A repeat instance of the same model and effort is run `<model>-<effort>-<n>`
(`prompts/make_prompts.py --instance 2 sol-high` writes its prompt) and is kept on
branch `runs/<model>/<effort>-<n>`, for example `runs/sol/high-2`.

## Where the runs happen

The runs execute in /Users/fillmore/NSF-NCAR/ODE/benchmarks, which every prompt
names as ROOT: this layout, plus `micm/` (MICM at e15898ce with its test build)
and one git repository per run. A finished run is fetched into its branch here:

    git fetch /Users/fillmore/NSF-NCAR/ODE/benchmarks/sol-high main:runs/sol/high

The scripts find files relative to their own location, so they work in the run
root; regenerate prompts there, not in this checkout.
