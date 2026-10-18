"""Check two tight independent references and install the tighter local files."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['winter-night-48h','summer-polluted-48h','tropical-clean-48h']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    log = ROOT/'scans/stress-reference.csv'
    with log.open() as file:
        reader = csv.DictReader(file); fields = reader.fieldnames; rows = list(reader)
    assert len(rows) == 6 and all(row['status'] == '0' for row in rows)
    manifest = {}
    for name in NAMES:
        paths = {level:ROOT/'runs/scans'/f'stress-reference-{level}'/name/'1' for level in [9,10]}
        comparisons = {}
        for level,other in [(9,10),(10,9)]:
            result = subprocess.run([sys.executable,str(ROOT.parent/'scenarios/score.py'),
                str(paths[level]/'ts1_snapshots.csv'),str(paths[other]/'ts1_snapshots.csv'),
                '--json','--level','8'],capture_output=True,text=True,check=True)
            summary = json.loads(result.stdout)
            comparisons[str(level)] = summary
            row = next(r for r in rows if r['id'] == f'stress-reference-{level}' and r['scenario'] == name)
            row.update(overall_digits=summary['overall_digits'],worst_species=summary['worst_species'][0][0],
                       worst_digits=summary['worst_species'][0][1],gate=summary['gate'],score_json=result.stdout.strip(),
                       accuracy_scope=f'convergence against tolerance 1e-{other}; not an independent truth error')
        directory = ROOT/'runs/stress'/name
        for basename in ['final','snapshots']:
            shutil.copyfile(paths[10]/f'ts1_{basename}.csv',directory/f'reference_{basename}.csv')
        manifest[name] = dict(scenario=(directory/'scenario.txt').read_text(),
            hashes={file:sha(directory/file) for file in ['cells.csv','photolysis.csv','reference_final.csv','reference_snapshots.csv']},
            convergence=comparisons)
        print(name,comparisons['9']['overall_digits'],comparisons['9']['worst_species'][0],flush=True)
    with log.open('w',newline='') as file:
        writer = csv.DictWriter(file,fieldnames=fields,lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    (ROOT/'validation/stress-references.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__ == '__main__':
    main()
