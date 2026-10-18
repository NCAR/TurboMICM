"""Render work-precision figures from recorded native runs (no solver timing)."""
import csv
import json
import os
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
os.environ['MPLCONFIGDIR'] = str(ROOT/'build'/'matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCENARIOS = ['tier1-30s','diurnal-30min','diurnal-10min','diurnal-5min','polluted-10min','clean-10min']
COLORS = dict(ros2='#8b6d2c',ros3='#7d53a3',ros4='#a35772',rodas3='#186cab',
              rodas4='#228163',rodas31='#d45d28',sdirk2='#7c7c7c')


def rows(path):
    with path.open() as file:
        return list(csv.DictReader(file))


def accuracy(row):
    return min(float(row['overall_digits']),float(row['worst_digits'])+1)


def save(fig,path):
    vector = path.with_suffix('.svg')
    fig.savefig(vector,metadata={'Date':None})
    vector.write_text('\n'.join(line.rstrip() for line in vector.read_text().splitlines())+'\n')
    fig.savefig(path.with_suffix('.png'),dpi=160)


def main():
    destination = ROOT/'figures'; destination.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.titlesize':12,'axes.spines.top':False,
                         'axes.spines.right':False,'svg.fonttype':'none',
                         'svg.hashsalt':'chemistry-work-precision'})
    data = []
    for path in sorted((ROOT/'scans').glob('[0-9][0-9]-*.csv')):
        data += [r for r in rows(path) if r.get('cells','1000') == '1000' and not
                 (path.name == '04-workprecision.csv' and r['id'] == 'rodas3-r0.0002-a18') and not
                 (path.name == '01-families.csv' and r['method'] == 'sdirk2')]
    fig,axes = plt.subplots(2,3,figsize=(13,7.2),layout='constrained')
    for ax,name in zip(axes.flat,SCENARIOS):
        for method,color in COLORS.items():
            selected = [r for r in data if r['scenario'] == name and r['method'] == method and
                        r['gate'] == 'PASS' and r['overall_digits']]
            if selected:
                ax.scatter([accuracy(r) for r in selected],[float(r['elapsed_ms']) for r in selected],
                           color=color,s=20,alpha=.65,label=method,edgecolors='none')
        for level in [2,3]:
            ax.axvline(level,color='#888888',ls='--',lw=.8,zorder=0)
        ax.set_title(name); ax.set_yscale('log'); ax.grid(axis='y',alpha=.2)
        ax.set_xlim(left=1)
    axes[0,0].legend(fontsize=8,ncol=2,frameon=False)
    fig.supxlabel('Ranking accuracy: min(overall digits, worst species digits + 1)')
    fig.supylabel('C elapsed time for 1,000 cells (ms, logarithmic scale)')
    fig.suptitle('Controlled numerical-method scans — individual exploratory measurements',fontsize=15)
    save(fig,destination/'work-precision-scans')
    plt.close(fig)

    config_path = ROOT/'scans'/'operating-points.json'
    large_path = ROOT/'scans'/'final-10000.csv'
    if not config_path.exists() or not large_path.exists():
        return
    points = json.loads(config_path.read_text())
    final = rows(large_path)
    baseline = rows(ROOT/'scans'/'final-baseline-10000.csv')
    fig,axes = plt.subplots(2,3,figsize=(13,7.2),layout='constrained')
    palette = ['#186cab','#d45d28','#228163','#7d53a3','#8b6d2c']
    for ax,scenario in zip(axes.flat,SCENARIOS):
        for point,color in zip(points,palette):
            selected = [r for r in final if r['id'] == point['id'] and r['scenario'] == scenario]
            assert len(selected) == 5 and all(r['gate'] == 'PASS' for r in selected),(scenario,point['id'])
            times = [float(r['elapsed_ms']) for r in selected]
            median = statistics.median(times)
            x = min(accuracy(r) for r in selected)
            ax.errorbar(x,median,yerr=[[median-min(times)],[max(times)-median]],
                        fmt='o',color=color,capsize=3,markersize=6,label=point['id'])
        reference = [r for r in baseline if r['scenario'] == scenario]
        assert len(reference) == 5,scenario
        failures = sum(r['status'] != '0' or r['gate'] != 'PASS' for r in reference)
        if failures:
            ax.text(.04,.96,f'MICM default: failed {failures}/5',transform=ax.transAxes,
                    va='top',fontsize=8,color='#555555',
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':2})
        else:
            times = [float(r['elapsed_ms']) for r in reference]
            median = statistics.median(times)
            ax.errorbar(min(accuracy(r) for r in reference),median,
                        yerr=[[median-min(times)],[max(times)-median]],
                        fmt='x',color='#333333',capsize=3,markersize=7,label='MICM default')
        for level in [2,3]:
            ax.axvline(level,color='#888888',ls='--',lw=.8,zorder=0)
        ax.set_title(scenario); ax.set_yscale('log'); ax.grid(axis='y',alpha=.2)
    axes[0,0].legend(fontsize=8,frameon=False)
    fig.supxlabel('Ranking accuracy on disclosed 1,000-cell prefix (digits)')
    fig.supylabel('Native elapsed time for 10,000 cells (ms, median and full range of 5)')
    fig.suptitle('Declared operating points — sequential native timing runs',fontsize=15)
    save(fig,destination/'work-precision-operating-points')
    plt.close(fig)


if __name__ == '__main__':
    main()
