"""Create 10,000-cell versions of the disclosed cases inside this checkout.

Only the supplied 1,000-cell accuracy reference is available here. score.py
checks accuracy on that prefix and the finite/nonnegative gate on ALL output
rows. These runs establish scaling and domain-wide gate behavior, not an
independent accuracy reference for the additional 9,000 cells.
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

root=Path(__file__).resolve().parents[1]
source=root.parent/'scenarios'
destination=root/'build/scaling'
cases=[('tier1-30s',2026,30,101,0,1),('diurnal-30min',101,1800,48,12,1),
       ('diurnal-10min',101,600,144,36,1),('diurnal-5min',101,300,288,72,1),
       ('polluted-10min',102,600,144,36,10),('clean-10min',103,600,144,36,.1)]
manifest=[]
for name,seed,dt,steps,snapshot,nox in cases:
    directory=destination/name
    command=[sys.executable,str(source/'make_scenario.py'),str(directory),'--cells','10000',
             '--seed',str(seed),'--dt',str(dt),'--steps',str(steps),
             '--constant' if not snapshot else '--diurnal','--nox-scale',str(nox)]
    if snapshot:command+=['--snapshot-every',str(snapshot)]
    subprocess.run(command,check=True)
    # Verify the generated prefix instead of merely assuming the recipe matches.
    with (directory/'cells.csv').open() as actual,(source/name/'cells.csv').open() as expected:
        for line in expected:
            assert next(actual)==line,(name,'generated cells differ from supplied prefix')
    for filename in ['reference_final.csv','reference_snapshots.csv']:
        original=source/name/filename
        if original.exists():
            link=directory/filename
            if not link.exists():link.symlink_to(original)
    manifest.append(dict(scenario=name,cells=10000,accuracy_reference_cells=1000,seed=seed,
                         dt=dt,steps=steps,nox=nox,snapshot_every=snapshot,
                         cells_sha256=hashlib.sha256((directory/'cells.csv').read_bytes()).hexdigest()))
(root/'results/scaling_inputs.json').write_text(json.dumps(manifest,indent=2)+'\n')
