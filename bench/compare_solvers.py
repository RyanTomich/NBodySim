#!/usr/bin/env python3
"""Benchmark all four registered solvers and parallel Barnes-Hut strong scaling.
Uses the production build/bench timer and shared forward Euler integrator.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import sys

from run import ROOT, read_state, trajectory_errors, timing_summary, write_csv


def execute(command, folder, env=None, stem='run'):
    result = subprocess.run([str(x) for x in command], cwd=ROOT, env=env,
                            capture_output=True, text=True)
    (folder/(stem+'.stdout.log')).write_text(result.stdout)
    (folder/(stem+'.stderr.log')).write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}); see {folder}/{stem}.stderr.log')
    return result.stdout, result.stderr


def metadata(stderr, prefix):
    for line in stderr.splitlines():
        if line.startswith(prefix+': '):
            return json.loads(line[len(prefix)+2:])
    return {}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bodies',type=int,nargs='+',default=[1000,4000,16000])
    parser.add_argument('--threads',type=int,nargs='+',default=[1,2,4,8])
    parser.add_argument('--samples',type=int,default=30)
    parser.add_argument('--warmups',type=int,default=3)
    parser.add_argument('--steps',type=int,default=20)
    parser.add_argument('--dt',type=float,default=0.0001)
    parser.add_argument('--theta',type=float,default=0.6)
    parser.add_argument('--order',type=int,default=4)
    parser.add_argument('--seed',type=int,default=42)
    args=parser.parse_args()
    if (min(args.bodies)<2 or min(args.threads)<1 or args.samples<2 or args.warmups<0
            or args.steps<2 or not math.isfinite(args.dt) or args.dt<=0
            or not math.isfinite(args.theta) or args.theta<0 or args.order not in range(1,7)):
        parser.error('invalid workload/settings; require at least two steps')
    args.threads=sorted(set([1]+args.threads))
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out=ROOT/'results'/f'four-solvers-{stamp}'
    out.mkdir(parents=True)
    status={'status':'running','settings':vars(args),'integrator':'production forward Euler',
            'scope':'simulate_with_solver only; includes integration and full solver calls',
            'cpu_count':os.cpu_count(),'platform':os.uname().sysname,'machine':os.uname().machine,
            'initialization':'generator uniform positions, equal total mass 1, v=0; direct n2 initial acceleration',
            'reference':'n2 with identical initial state, dt and steps',
            'thread_environment':{'OMP_DYNAMIC':'FALSE','OMP_PROC_BIND':'close','OMP_PLACES':'cores'}}
    status_path=out/'status.json'
    def save_status(): status_path.write_text(json.dumps(status,indent=2,allow_nan=False)+'\n')
    save_status()
    print(f'Results: {out}',flush=True)
    rows=[]
    try:
        execute(['make','-B','all'],out,stem='build')
        shutil.copy2(ROOT/'build/bench',out/'bench')
        config,_=execute(['make','print-config'],out,stem='build-config')
        status['build_config']=config
        for count in args.bodies:
            case=out/str(count); case.mkdir()
            # Generator output stays in ignored bench/data; retain the actual input in results.
            raw=ROOT/'bench/data'/f'instance-suite-n{count}-seed{args.seed}-{stamp}.txt'
            execute([sys.executable,'generate_instance.py','--n',count,'--L',1,
                     '--w-min',1/count,'--w-max',1/count,'--v-max',0,'--a-max',0,
                     '--seed',args.seed,'--output',raw],case,stem='generate')
            initial_text,_=execute([ROOT/'build/app','--solver','n2',raw,args.dt,1],case,stem='initialize')
            initial=case/'input.txt'; initial.write_text(initial_text)
            raw_state=read_state(raw); initial_state=read_state(initial)
            assert [r[:7] for r in raw_state]==[r[:7] for r in initial_state]
            configurations=[('serial-barnes-hut',1,'theta',args.theta),('fmm',1,'order',args.order)]
            configurations += [('parallel-barnes-hut',t,'theta',args.theta) for t in args.threads]
            random.Random(args.seed+count).shuffle(configurations)
            configurations.insert(0,('n2',1,None,None))
            reference=None; force_reference=None; parallel_reference=None
            for solver,threads,key,value in configurations:
                folder=case/f'{solver}-{threads}threads'; folder.mkdir()
                env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',
                         OMP_PROC_BIND='close',OMP_PLACES='cores')
                option=[] if key is None else [key,value]
                samples_text,stderr=execute([out/'bench',initial,args.dt,args.steps,args.samples,
                     args.warmups,'timing',folder/'final.txt',solver,*option],folder,env)
                (folder/'samples.csv').write_text(samples_text)
                samples=[float(r['seconds']) for r in csv.DictReader(samples_text.splitlines())]
                if len(samples)!=args.samples or any(not math.isfinite(x) or x<=0 for x in samples):
                    raise ValueError('Invalid timing samples')
                final=read_state(folder/'final.txt')
                if solver=='n2': reference=final
                if solver=='parallel-barnes-hut':
                    if parallel_reference is None: parallel_reference=final
                    elif final!=parallel_reference: raise ValueError('Thread count changed final state')
                errors=trajectory_errors(final,reference)
                write_csv(folder/'trajectory-error.csv',list(errors[0]),errors)
                # v=0 means first Euler drift leaves positions unchanged: compare forces at shared positions.
                force_text,_=execute([ROOT/'build/app','--solver',solver,
                    *([] if key is None else ['--'+key,value]),initial,args.dt,1],folder,env,stem='force')
                force_file=folder/'force-state.txt'; force_file.write_text(force_text)
                force_state=read_state(force_file)
                assert [r[1:4] for r in force_state]==[r[1:4] for r in initial_state]
                if solver=='n2': force_reference=force_state
                force_error=trajectory_errors(force_state,force_reference)[2]
                summary=timing_summary(samples,args.steps)
                mean=summary['mean_seconds_per_step']*1000
                row=dict(bodies=count,solver=solver,threads=threads,steps=args.steps,dt=args.dt,
                         theta=args.theta if key=='theta' else '',order=args.order if key=='order' else '',
                         samples=args.samples,mean_ms_per_step=mean,
                         stddev_ms_per_step=summary['stddev_seconds_per_step']*1000,
                         force_relative_l2=force_error['relative_l2_error'],force_rms=force_error['rms_error'],
                         position_rms=errors[0]['rms_error'],velocity_rms=errors[1]['rms_error'])
                totals=metadata(stderr,'solver_totals')
                if solver=='parallel-barnes-hut' and totals.get('actual_threads')!=threads:
                    raise ValueError('OpenMP did not use the requested thread count')
                rows.append(row)
                details=dict(status='success',solver=solver,threads=threads,
                             solver_config=metadata(stderr,'solver_config'),
                             solver_totals=metadata(stderr,'solver_totals'),**summary)
                (folder/'status.json').write_text(json.dumps(details,indent=2)+'\n')
                write_csv(out/'comparison.csv',list(rows[0]),rows)
                print(f'{count:6d} {solver:21s} threads={threads:2d}: {mean:9.4f} ms/step; '
                      f'force rel.L2={row["force_relative_l2"]:.3g}; position RMS={row["position_rms"]:.3g}',flush=True)
        scaling=[]
        for count in args.bodies:
            group=[r for r in rows if r['bodies']==count and r['solver']=='parallel-barnes-hut']
            baseline=next(r['mean_ms_per_step'] for r in group if r['threads']==1)
            for r in sorted(group,key=lambda r:r['threads']):
                speedup=baseline/r['mean_ms_per_step']
                scaling.append(dict(bodies=count,threads=r['threads'],mean_ms_per_step=r['mean_ms_per_step'],
                                    speedup_vs_1thread=speedup,parallel_efficiency=speedup/r['threads']))
        write_csv(out/'scaling.csv',list(scaling[0]),scaling)
        report=['# Four-solver comparison','',
                f'Theta {args.theta} for both Barnes-Hut solvers; FMM order {args.order}.',
                f'{args.samples} samples, {args.warmups} warmups, {args.steps} Euler steps, dt={args.dt}.',
                'Serial solvers use one thread. Only parallel Barnes-Hut uses multiple threads.',
                'Timing includes integration and all tree/force work; excludes setup/reset/I/O.',
                'Force error is measured at shared initial positions; trajectory error against direct n2.',
                '', '| Bodies | Solver | Threads | ms/step | Force relative L2 | Position RMS |',
                '|---|---|---|---|---|---|']
        for r in sorted(rows,key=lambda r:(r['bodies'],r['solver'],r['threads'])):
            report.append(f'| {r["bodies"]} | {r["solver"]} | {r["threads"]} | '
                          f'{r["mean_ms_per_step"]:.4f} | {r["force_relative_l2"]:.3g} | {r["position_rms"]:.3g} |')
        report+=['','## Strong scaling','','| Bodies | Threads | Speedup vs 1 thread | Efficiency |',
                 '|---|---|---|---|']
        for r in scaling:
            report.append(f'| {r["bodies"]} | {r["threads"]} | {r["speedup_vs_1thread"]:.2f} | {r["parallel_efficiency"]:.1%} |')
        report+=['','These results describe this uniform distribution and machine, not a universal ranking.',
                 'Equal theta does not imply equal accuracy. Scaling may saturate or regress at higher thread counts.',
                 'Workspace is prepared once and reused across warmups and samples, as in the existing benchmark driver.']
        (out/'summary.md').write_text('\n'.join(report)+'\n')
        status['status']='success'
    except Exception as exc:
        status.update(status='failed',error=str(exc)); raise
    finally:
        save_status()


if __name__=='__main__':
    main()
