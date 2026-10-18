"""Publication-style work–precision figures from measured native runs."""
import csv
from collections import defaultdict
import os
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
os.environ['MPLCONFIGDIR'] = str(ROOT / 'build' / 'matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAMES = ['tier1-30s', 'diurnal-30min', 'diurnal-10min', 'diurnal-5min', 'polluted-10min', 'clean-10min']
POINTS = ['fast', 'balanced', 'accurate', 'recommended', 'precise']


def read(name):
    with (ROOT / 'scans' / name).open() as f:
        return list(csv.DictReader(f))


def main():
    rows = read('final-clang-1000.csv')
    reference = read('micm-timings.csv')
    dest = ROOT / 'plots'; dest.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none', 'svg.hashsalt': 'work-precision', 'savefig.facecolor': 'white'})
    for metric, label, filename in [('overall_digits', 'Overall significant digits', 'work_precision'),
                                    ('worst_digits', 'Worst-species significant digits', 'worst_species_precision')]:
        fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout='constrained')
        for name, ax in zip(NAMES, axes.flat):
            group = defaultdict(list)
            for row in rows:
                if row['scenario'] == name and row['status'] == 'OK':
                    group[row['configuration']].append(row)
            x, y, lower, upper = [], [], [], []
            for point in POINTS:
                rs = group[point]; times = [float(r['elapsed_ms']) for r in rs]
                x.append(float(rs[0][metric])); y.append(median(times))
                lower.append(y[-1]-min(times)); upper.append(max(times)-y[-1])
            ax.errorbar(x, y, yerr=[lower, upper], color='#126682', marker='o', lw=1.7,
                        markersize=4, capsize=2, label='C RODAS3: five operating points')
            index = POINTS.index('recommended')
            ax.plot(x[index], y[index], 'o', color='#d57a1f', markersize=7, label='Recommended')
            for config, color, marker, text in [('MICM3-default', '#b84540', 'x', 'MICM3 defaults'),
                                               ('MICM3-controlled', '#333333', '*', 'MICM3 with accuracy control')]:
                rs = [r for r in reference if r['scenario'] == name and r['configuration'] == config and r['status'] == 'OK']
                if rs and float(rs[0][metric]) >= (1 if metric == 'overall_digits' else 0):
                    ax.plot(float(rs[0][metric]), median(float(r['elapsed_ms']) for r in rs),
                            marker, color=color, markersize=8, label=text)
                elif config == 'MICM3-default':
                    ax.text(.03, .95, 'Default MICM3: fails or inaccurate', transform=ax.transAxes,
                            fontsize=8, color=color, va='top')
            for level in ([2, 3] if metric == 'overall_digits' else [1, 2]):
                ax.axvline(level, color='#b7b7b7', lw=.8, ls='--', zorder=0)
            ax.set_title(name, loc='left', fontweight='bold', fontsize=11)
            ax.set_yscale('log'); ax.grid(axis='y', alpha=.15)
            ax.set_xlabel(label); ax.set_ylabel('Native elapsed time (ms)')
        handles, labels = axes[0, 0].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        fig.legend(by_label.values(), by_label.keys(), loc='outside lower center', ncol=2,
                   frameon=False, fontsize=10)
        fig.suptitle('TS1: 1,000 cells, 10 workers, clang -O3 -DNDEBUG\n'
                     'Median of five native runs; bars span observed minimum–maximum', fontsize=13)
        fig.savefig(dest / f'{filename}.png', dpi=180)
        svg = dest / f'{filename}.svg'
        fig.savefig(svg, metadata={'Date': None})
        svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
        plt.close(fig)
    print('wrote plots/work_precision and plots/worst_species_precision (.png, .svg)')


if __name__ == '__main__':
    main()
