#!/usr/bin/env python3
"""Export scientific work–precision figures from C run records (never Python time)."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.environ['MPLCONFIGDIR']=str(ROOT/'.cache/matplotlib')
import csv
import json
from collections import defaultdict
from statistics import median
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scan import NAMES

METHODS=['ros2','ros3','ros4','rodas3','rodas4','extrap2','extrap3','rodas31','micm3','MICM-three-stage']
LABELS=['ROS2','ROS3','ROS4','RODAS3','RODAS4','Euler extrap. 2','Euler extrap. 3','RODAS3.1','C MICM3','C++ MICM3']
SUITES={'family','family-corrected','loose','mid','rodas31','baseline','tolerance-isolation','micm-baseline'}

def main():
    rows=[]
    for path in sorted((ROOT/'scans').glob('*.csv')):
        if path.stem in ('smoke','baseline') or path.stem.startswith('diagnostic'): continue
        for row in csv.DictReader(path.open()):
            if row.get('message','').startswith('Provisional:') or row.get('compiler')=='gcc-15': continue
            if row.get('returncode')!='0' or not row.get('overall_digits') or row.get('scenario') not in NAMES: continue
            row['config']=json.loads(row['configuration']);rows.append(row)
    colors=dict(zip(METHODS,plt.get_cmap('tab10').colors))
    out=ROOT/'plots';out.mkdir(exist_ok=True)
    for metric,suffix,axis,levels in [('overall_digits','','Overall significant digits',[2,3]),
                                     ('worst_digits','_species','Worst-species significant digits',[1,2])]:
        fig,axes=plt.subplots(2,3,figsize=(15,8.7),constrained_layout=True)
        for ax,name in zip(axes.flat,NAMES):
            for method,label in zip(METHODS,LABELS):
                groups=defaultdict(list)
                for r in rows:
                    if r['scenario']==name and r['config']['method']==method and r['suite'] in SUITES:
                        if float(r[metric])<0: continue
                        groups[r['configuration']].append(r)
                series=[]
                for records in groups.values():
                    series.append((float(records[0][metric]),median(float(r['elapsed_ms']) for r in records)))
                series.sort()
                if series:
                    ax.plot([x for x,y in series],[y for x,y in series],'.--' if method=='MICM-three-stage' else '.-',
                            color=colors[method],label=label,lw=1.4,ms=5)
            points=defaultdict(list)
            for r in rows:
                if r['scenario']==name and r['suite']=='operating-points' and r['compiler']=='clang':
                    points[r['point']].append(r)
            if points:
                pairs=[(float(rs[0][metric]),median(float(r['elapsed_ms']) for r in rs)) for rs in points.values()]
                ax.scatter([x for x,y in pairs],[y for x,y in pairs],marker='*',s=90,c='black',
                           label='declared points (median 5)',zorder=5)
            ax.set_title(name);ax.set_yscale('log');ax.grid(alpha=.25)
            ax.set_xlabel(axis);ax.set_ylabel('C/C++ elapsed time (ms, log scale)')
            for level in levels: ax.axvline(level,color='gray',ls=':',lw=1)
        handles,labels=axes[0,0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='outside lower center',ncol=6,frameon=False,fontsize=9)
        fig.savefig(out/f'work_precision{suffix}.png',dpi=160)
        fig.savefig(out/f'work_precision{suffix}.pdf')
        plt.close(fig)
    print('Overall and worst-species work–precision PNG/PDF figures generated')
if __name__=='__main__': main()
