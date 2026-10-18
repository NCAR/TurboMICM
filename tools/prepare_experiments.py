"""Generate larger and independent stress inputs, inside this repository only."""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT.parent/'scenarios'/'make_scenario.py'


def create(directory,name,cells,seed,dt,steps,extra):
    subprocess.run([sys.executable,str(GENERATOR),str(directory/name),
                    '--cells',str(cells),'--seed',str(seed),'--dt',str(dt),
                    '--steps',str(steps),*extra],check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--large',action='store_true')
    parser.add_argument('--stress',action='store_true')
    args = parser.parse_args()
    if args.large:
        directory = ROOT/'runs'/'large'
        for name,seed,dt,steps,scale in [
            ('tier1-30s',2026,30,101,1),('diurnal-30min',101,1800,48,1),
            ('diurnal-10min',101,600,144,1),('diurnal-5min',101,300,288,1),
            ('polluted-10min',102,600,144,10),('clean-10min',103,600,144,.1)]:
            extra = ['--constant'] if name == 'tier1-30s' else [
                '--diurnal','--nox-scale',str(scale),'--snapshot-every',str(21600//dt)]
            create(directory,name,10000,seed,dt,steps,extra)
            # Verify the disclosed prefix before using it for accuracy scoring.
            generated = (directory/name/'cells.csv').read_text().splitlines()
            supplied = (ROOT.parent/'scenarios'/name/'cells.csv').read_text().splitlines()
            assert generated[:len(supplied)] == supplied,name
    if args.stress:
        directory = ROOT/'runs'/'stress'
        settings = [
            ('winter-night-48h',807,1800,96,60,10,21,1),
            ('summer-polluted-48h',808,600,288,65,172,12,20),
            ('tropical-clean-48h',809,300,576,-10,355,4,.03)]
        for name,seed,dt,steps,lat,doy,start,nox in settings:
            create(directory,name,96,seed,dt,steps,[
                '--diurnal','--lat',str(lat),'--doy',str(doy),'--start',str(start),
                '--nox-scale',str(nox),'--snapshot-every',str(21600//dt)])


if __name__ == '__main__':
    main()
