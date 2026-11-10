#!/usr/bin/env python3
"""Render standalone performance and scaling plots from compare_solvers results."""
import argparse
import csv
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'nbody-matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results',type=Path)
    args=parser.parse_args()
    rows=list(csv.DictReader((args.results/'comparison.csv').open()))
    status=json.loads((args.results/'status.json').read_text())
    integrator='leapfrog' if 'leapfrog' in status.get('integrator','').lower() else 'Euler'
    scaling=list(csv.DictReader((args.results/'scaling.csv').open()))
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    keys=sorted(set((r['solver'],int(r['threads'])) for r in rows))
    for solver,threads in keys:
        group=sorted((r for r in rows if r['solver']==solver and int(r['threads'])==threads),key=lambda r:int(r['bodies']))
        label=solver+(f' ({threads} threads)' if solver=='parallel-barnes-hut' else '')
        ax.errorbar([int(r['bodies']) for r in group],[float(r['mean_ms_per_step']) for r in group],
                    yerr=[float(r['stddev_ms_per_step']) for r in group],marker='o',capsize=3,label=label)
    ax.set(xscale='log',yscale='log',xlabel='Number of bodies',ylabel=f'Mean runtime (ms / {integrator} step)',title='Four solver performance; error bars = sample standard deviation')
    ax.grid(True,alpha=.2); ax.legend(fontsize=8)
    for ext in ('png','pdf'): fig.savefig(args.results/f'performance.{ext}',dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,5),layout='constrained')
    for count in sorted(set(int(r['bodies']) for r in scaling)):
        group=sorted((r for r in scaling if int(r['bodies'])==count),key=lambda r:int(r['threads']))
        ax.plot([int(r['threads']) for r in group],[float(r['speedup_vs_1thread']) for r in group],marker='o',label=f'{count:,} bodies')
    threads=sorted(set(int(r['threads']) for r in scaling))
    ax.plot(threads,threads,'k--',alpha=.4,label='Ideal linear scaling')
    ax.set(xlabel='OpenMP threads',ylabel='Speedup relative to one thread',title='Parallel Barnes-Hut strong scaling',xticks=threads)
    ax.grid(True,alpha=.2); ax.legend()
    for ext in ('png','pdf'): fig.savefig(args.results/f'scaling.{ext}',dpi=180)
    plt.close(fig)
    print(f'Plots saved to {args.results}')


if __name__=='__main__': main()
