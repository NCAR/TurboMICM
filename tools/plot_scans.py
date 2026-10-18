#!/usr/bin/env python3
"""Generate static scientific work–precision plots from recorded C trials."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.environ['MPLCONFIGDIR']=str(ROOT/'.cache/matplotlib')
import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAMES=['tier1-30s','diurnal-30min','diurnal-10min','diurnal-5min','polluted-10min','clean-10min']

def main():
    rows=[]
    for path in sorted((ROOT/'scans').glob('*.csv')):
        if path.stem in ('smoke','baseline') or path.stem.startswith('diagnostic'): continue
        for row in csv.DictReader(path.open()):
            if row['returncode']=='0' and row.get('overall_digits') and row['scenario'] in NAMES:
                row['config']=json.loads(row['configuration']);rows.append(row)
    fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
    colors=dict(zip(['ros2','ros3','ros4','rodas3','rodas4','extrap2','extrap3','micm3'],plt.get_cmap('tab10').colors))
    for ax,name in zip(axes.flat,NAMES):
        for method,color in colors.items():
            series=[r for r in rows if r['scenario']==name and r['config']['method']==method and r['suite'] in ('family','loose','baseline')]
            series.sort(key=lambda r:float(r['overall_digits']))
            if series: ax.plot([float(r['overall_digits']) for r in series],[float(r['elapsed_ms']) for r in series],'.-',color=color,label=method)
        points=[r for r in rows if r['scenario']==name and r['suite']=='operating-points']
        if points: ax.scatter([float(r['overall_digits']) for r in points],[float(r['elapsed_ms']) for r in points],marker='*',s=90,c='black',label='declared points',zorder=5)
        ax.set_title(name);ax.set_yscale('log');ax.grid(alpha=.25)
        ax.set_xlabel('Overall significant digits');ax.set_ylabel('C elapsed time (ms, log scale)')
        ax.axvline(2,color='gray',ls=':',lw=1);ax.axvline(3,color='gray',ls=':',lw=1)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=9,frameon=False)
    out=ROOT/'plots';out.mkdir(exist_ok=True)
    fig.savefig(out/'work_precision.png',dpi=160)
    fig.savefig(out/'work_precision.pdf')
    print(out/'work_precision.png')
if __name__=='__main__': main()
