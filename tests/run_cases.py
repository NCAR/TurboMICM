"""Run all named MICM case ports, preserving case identifiers and failure status."""
import argparse
import fnmatch
import importlib
from pathlib import Path
import sys
import traceback
from harness import CASE_LIST,REGISTRY,NOT_APPLICABLE

for file in sorted(Path(__file__).parent.glob('cases_*.py')):
    importlib.import_module(file.stem)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--filter',default='*')
parser.add_argument('--traceback',action='store_true')
args = parser.parse_args()
counts = {'PASS':0,'FAIL':0,'N/A':0}
for id_ in CASE_LIST:
    if not fnmatch.fnmatchcase(id_,args.filter):
        continue
    if id_ in NOT_APPLICABLE:
        status,reason = 'N/A',NOT_APPLICABLE[id_]
    elif id_ not in REGISTRY:
        status,reason = 'FAIL','case port not yet implemented'
    else:
        try:
            REGISTRY[id_](id_)
            status,reason = 'PASS',''
        except Exception as error:
            status,reason = 'FAIL',' '.join(str(error).split())[:500] or type(error).__name__
            if args.traceback:
                traceback.print_exc(file=sys.stderr)
    counts[status] += 1
    print(f'{id_} {status} {reason}'.rstrip(),flush=True)
print('Cases: '+', '.join(f'{v} {k}' for k,v in counts.items()),file=sys.stderr)
sys.exit(bool(counts['FAIL']))
