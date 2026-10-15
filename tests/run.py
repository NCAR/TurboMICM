"""Report every required case; missing ports are failures, never implicit passes."""
import argparse
import sys
import traceback
from collections import Counter
from pathlib import Path

from harness import CASES, NotApplicable
import rate_cases
import state_cases
import system_cases
import matrix_cases
import sparse_cases
import grouped_cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--filter',default='')
    parser.add_argument('--traceback',action='store_true')
    args = parser.parse_args()
    totals = Counter()
    for line in (Path(__file__).resolve().parents[1]/'case-list.txt').read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        identifier,kind = line.split()
        if args.filter not in identifier:
            continue
        function = CASES.get(identifier)
        if function is None:
            state,reason = 'FAIL','case port is not implemented'
        else:
            try:
                reason = function(identifier) or 'assertions passed'
                state = 'PASS'
            except NotApplicable as error:
                state,reason = 'N/A',str(error)
            except Exception as error:
                state,reason = 'FAIL',str(error).replace('\n',' ')
                if args.traceback:
                    traceback.print_exc(file=sys.stderr)
        print(f'{identifier} {state} {reason}',flush=True)
        totals[state] += 1
    print(' '.join(f'{key}={value}' for key,value in sorted(totals.items())),file=sys.stderr)
    return int(bool(totals['FAIL']))


if __name__ == '__main__':
    sys.exit(main())
