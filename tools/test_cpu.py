"""Run ported cases and report unported cases as failures, never as passes."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CASE_LIST = ROOT / '../reference/micm_cases.txt'
# Each listed executable must check every assertion in each mapped case.
CASES = {
    **{case: ('build/test_rates', 'rate calculation') for case in [
        'arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithSystem',
        'arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithPrescribedArguments',
        'tunneling_rate_constant:TunnelingRateConstant.CalculateWithMinimalArguments',
        'tunneling_rate_constant:TunnelingRateConstant.CalculateWithAllArguments',
        'user_defined_rate_constant:UserDefinedRateConstant.CalculateWithSystem',
        'user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRate',
        'user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRateAndName',
        'user_defined_rate_constant:UserDefinedRateConstant.CustomScalingFactor',
    ]},
    'backward_euler:BackwardEuler.CanCallSolve': ('build/test_backward_euler', 'chain solve'),
    'backward_euler:BackwardEuler.IsConverged': ('build/test_backward_euler', 'convergence thresholds'),
}


def main():
    results = {}
    failures = 0
    for line in CASE_LIST.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        case, _kind = line.split()
        if case not in CASES:
            status, reason = 'FAIL', 'not yet ported'
        else:
            executable, reason = CASES[case]
            if executable not in results:
                results[executable] = subprocess.run([str(ROOT / executable)],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).returncode
            status = 'PASS' if results[executable] == 0 else 'FAIL'
            if status == 'FAIL': reason = 'ported check failed'
        failures += status == 'FAIL'
        print(case, status, reason)
    return bool(failures)

if __name__ == '__main__':
    sys.exit(main())
