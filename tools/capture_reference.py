"""Run supplied test executables only from a private copy of their build tree."""
from pathlib import Path
import argparse
import shutil
import subprocess

p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,default=Path('../micm/build'))
p.add_argument('--copy',type=Path,default=Path('build/reference-tests'))
p.add_argument('--output',type=Path,default=Path('tests/reference'))
args=p.parse_args()
if not args.copy.exists():
    shutil.copytree(args.source,args.copy,symlinks=True)
args.output.mkdir(parents=True,exist_ok=True)
for line in Path('tests/cases.txt').read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    identifier,kind=line.split()
    if kind!='program':
        continue
    executable=(args.copy/('test_'+identifier)).resolve()
    result=subprocess.run([str(executable)],cwd=args.copy,capture_output=True,text=True,check=True)
    (args.output/(identifier+'.txt')).write_text(result.stdout)
    if result.stderr:
        (args.output/(identifier+'.stderr.txt')).write_text(result.stderr)
    print(identifier, len(result.stdout.splitlines()), 'lines')
