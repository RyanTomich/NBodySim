#!/usr/bin/env python3
"""Plot runtime, force accuracy, and their tradeoff for a theta sweep."""
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
    sizes=sorted(set(int(r['bodies']) for r in rows))
    fig,axes=plt.subplots(len(sizes),3,figsize=(15,4.2*len(sizes)),squeeze=False,layout='constrained')
    keys=sorted(set((r['solver'],int(r['threads'])) for r in rows if r['theta']!=''))
    for index,count in enumerate(sizes):
        selected=[r for r in rows if int(r['bodies'])==count]
        direct=next(float(r['mean_ms_per_step']) for r in selected if r['solver']=='n2')
        ax_time,ax_error,ax_tradeoff=axes[index]
        for solver,threads in keys:
            group=sorted((r for r in selected if r['solver']==solver and int(r['threads'])==threads),key=lambda r:float(r['theta']))
            if not group: continue
            theta=[float(r['theta']) for r in group]
            times=[float(r['mean_ms_per_step']) for r in group]
            error=[float(r['force_relative_l2']) for r in group]
            label=('Serial BH' if solver=='serial-barnes-hut' else f'Parallel BH ({threads} threads)')
            ax_time.errorbar(theta,times,yerr=[float(r['stddev_ms_per_step']) for r in group],marker='o',capsize=2,label=label)
            if threads==1: ax_error.plot(theta,error,marker='o',label=label)
            ax_tradeoff.plot(error,times,marker='o',label=label)
            if threads==1:
                for x,y,t in zip(error,times,theta): ax_tradeoff.annotate(f'{t:g}',(x,y),xytext=(3,3),textcoords='offset points',fontsize=7)
        ax_time.axhline(direct,color='gray',linestyle='--',label='Direct n2')
        ax_time.set(xlabel='Theta',ylabel=f'Mean ms / {integrator} step',yscale='log',title=f'{count:,} bodies: runtime (± sample SD)')
        ax_error.set(xlabel='Theta',ylabel='Relative force L2 error',yscale='log',title=f'{count:,} bodies: force accuracy')
        ax_tradeoff.set(xlabel='Relative force L2 error',ylabel=f'Mean ms / {integrator} step',xscale='log',yscale='log',title=f'{count:,} bodies: speed vs error (labels = theta)')
        for ax in axes[index]: ax.grid(True,alpha=.2); ax.legend(fontsize=7)
    for extension in ('png','pdf'): fig.savefig(args.results/f'theta_sweep.{extension}',dpi=180)
    plt.close(fig)
    # Best measured settings on the sampled grid, without interpolation.
    budgets=[1e-4,1e-3,1e-2]
    report=['# Best sampled settings at shared force-error budgets','',
            'Each row is the fastest measured setting satisfying the budget; no interpolation.',
            'These are grid choices, not an optimized or universal accuracy guarantee.','',
            '| Bodies | Force relative L2 budget | Solver | Threads | Selected theta | ms/step | Actual force error |',
            '|---|---|---|---|---|---|---|']
    for count in sizes:
        for budget in budgets:
            for solver,threads in keys:
                candidates=[r for r in rows if int(r['bodies'])==count and r['solver']==solver and int(r['threads'])==threads and float(r['force_relative_l2'])<=budget]
                if not candidates: continue
                r=min(candidates,key=lambda r:float(r['mean_ms_per_step']))
                report.append(f'| {count} | {budget:g} | {solver} | {threads} | {r["theta"]} | {float(r["mean_ms_per_step"]):.4f} | {float(r["force_relative_l2"]):.3g} |')
    (args.results/'accuracy_budgets.md').write_text('\n'.join(report)+'\n')
    print(f'Plots and shared-error-budget comparison saved to {args.results}')


if __name__=='__main__': main()
