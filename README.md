# TurboMICM

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23618864.svg)](https://doi.org/10.5281/zenodo.23618864)

Optimized ODE solvers for the Model Independent Chemistry Module.

This repository tracks the MICM-to-C experiment, a controlled comparison of LLM
coding agents. Every agent gets the same prompt: port MICM's CPU solver to
portable C (C11, libm, and OpenMP) and make it fast without depending on one chip,
compiler, or operating system.

## Citation

Please cite the archived TurboMICM release when using this software or experiment:

Fillmore, D. (2026). *NCAR/TurboMICM: v2026.10.001-beta* [Software]. Zenodo.
[https://doi.org/10.5281/zenodo.23618864](https://doi.org/10.5281/zenodo.23618864).

This DOI identifies release `v2026.10.001-beta`, published on 2026-10-24.
Machine-readable citation metadata is in [CITATION.cff](CITATION.cff).

## REAPER

The framework that the contest feeds into takes mechanism-specific kernels to a
host model by its run-time path, REAPER: Run-time: Emit, Assemble, Plug in,
Execute, Reuse. Emit straight-line C from the mechanism the model loaded;
assemble it with the system compiler; plug the kernels in through function
pointers, after a fingerprint check, with the generic solver as the fallback;
execute; and reuse the shared library, cached under a hash of the mechanism,
the compiler, and the flags. A reaper cuts a field into straight rows; REAPER
cuts a mechanism's reaction graph into straight-line code. The design is on the
wiki: [From contest to framework](https://github.com/NCAR/TurboMICM/wiki/From-Contest-to-Framework).

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
