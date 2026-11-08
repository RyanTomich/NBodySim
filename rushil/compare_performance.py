#!/usr/bin/env python3
"""Compare src and rushil Barnes-Hut with identical saved inputs and leapfrog.
Run from any directory. Results and binaries are retained in results/.
"""
import argparse
import csv
import datetime
import json
import math
import os
from pathlib import Path
import random
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def run(command, output, env=None):
    result = subprocess.run([str(x) for x in command], cwd=ROOT, env=env,
                            capture_output=True, text=True)
    output.write_text(result.stdout)
    output.with_suffix(output.suffix + '.stderr').write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}): {command}; see {output}')
    return result.stdout


def read_state(path):
    lines = path.read_text().splitlines()
    data = [[float(x) for x in row.split()] for row in lines[1:]]
    if len(data) != int(lines[0]) or any(len(r) != 10 or not all(map(math.isfinite, r)) for r in data):
        raise ValueError(f'Invalid final state: {path}')
    return data


def errors(actual, reference):
    def rms(indices):
        return math.sqrt(sum((a[k]-b[k])**2 for a,b in zip(actual, reference)
                             for k in indices)/len(actual))
    mass = sum(r[0] for r in actual)
    def cm(state, indices):
        return [sum(r[0]*r[k] for r in state)/mass for k in indices]
    return dict(position_rms=rms(range(1,4)), velocity_rms=rms(range(4,7)),
                cm_position_error=math.dist(cm(actual,range(1,4)),cm(reference,range(1,4))),
                cm_velocity_error=math.dist(cm(actual,range(4,7)),cm(reference,range(4,7))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bodies', nargs='+', type=int, default=[1000,4000])
    parser.add_argument('--steps', type=int, default=50)
    parser.add_argument('--dt', type=float, default=0.0001)
    parser.add_argument('--theta', type=float, default=0.6)
    parser.add_argument('--samples', type=int, default=5)
    parser.add_argument('--warmups', type=int, default=1)
    parser.add_argument('--threads', nargs='+', type=int, default=[1,4])
    parser.add_argument('--base-cases', nargs='+', type=int, default=[2000],
                        help='summarize task cutoffs; 0 tasks every internal node; 2147483647 is serial')
    args = parser.parse_args()
    if (min(args.base_cases)<0 or max(args.base_cases)>2147483647 or min(args.bodies)<2 or min(args.threads)<1 or args.steps<1 or args.samples<2
            or args.warmups<0 or not math.isfinite(args.dt) or args.dt<=0
            or not math.isfinite(args.theta) or args.theta<=0):
        parser.error('Invalid count, timestep, theta, samples, or warmups')
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out=ROOT/'results'/f'barnes-hut-src-vs-rushil-{stamp}'
    out.mkdir(parents=True)
    status=dict(status='running', settings=vars(args), integrator='matched leapfrog',
                scope='integration + initial acceleration + per-step tree maintenance/forces; no setup/reset/I/O',
                compiler_flags='-O3 -std=c17 -fno-fast-math; no LTO; OpenMP enabled for both drivers',
                platform=os.uname().sysname, machine=os.uname().machine)
    status_path=out/'status.json'
    status_path.write_text(json.dumps(status,indent=2))
    print(f'Results: {out}',flush=True)
    try:
        compiler=os.environ.get('CC','clang')
        omp=(['-Xclang','-fopenmp','-isystem','/opt/homebrew/opt/libomp/include']
             if os.uname().sysname=='Darwin' else ['-fopenmp'])
        libs=(['-L/opt/homebrew/opt/libomp/lib','-lomp'] if os.uname().sysname=='Darwin' else ['-fopenmp'])
        flags=[compiler,'-O3','-std=c17','-Wall','-Wextra','-Wpedantic','-fno-fast-math',*omp]
        driver=ROOT/'rushil/bench/compare_driver.c'
        run([*flags,'-Isrc',driver,*sorted(p for p in (ROOT/'src').glob('*.c') if p.name!='main.c'),
             *libs,'-lm','-o',out/'src-driver'],out/'build-src.log')
        run([*flags,'-DNAIVE_BACKEND',driver,ROOT/'rushil/n2.c',ROOT/'rushil/state.c',
             ROOT/'rushil/physics.c',*libs,'-lm','-o',out/'naive-driver'],out/'build-naive.log')
        for cutoff in args.base_cases:
            run([*flags,'-DRUSHIL_BACKEND',f'-DSUMMARIZE_BASE_CASE={cutoff}',driver,
                 ROOT/'rushil/state.c',ROOT/'rushil/physics.c',*libs,'-lm',
                 '-o',out/f'rushil-driver-{cutoff}'],out/f'build-rushil-{cutoff}.log')
        summaries=[]
        for n in args.bodies:
            case=out/str(n); case.mkdir()
            rng=random.Random(42)
            input_path=case/'input.txt'
            with input_path.open('w') as f:
                f.write(f'{n}\n')
                for _ in range(n):
                    row=[1/n,*[rng.uniform(-1,1) for _ in range(3)],0,0,0,0,0,0]
                    f.write(' '.join(format(x,'.17g') for x in row)+'\n')
            configurations=[('n2',1,'naive-driver','n2',''),('src',1,'src-driver','barnes-hut','')]
            variants=[('rushil',t,f'rushil-driver-{cutoff}','barnes-hut',cutoff)
                      for t in args.threads for cutoff in args.base_cases]
            # Avoid always measuring thresholds in ascending order.
            random.Random(42+n).shuffle(variants)
            configurations += variants
            cutoff_reference=None
            reference=None
            baseline=None
            n2_baseline=None
            for name,threads,binary,backend,cutoff in configurations:
                label=f'{name}-{threads}threads'+(f'-base{cutoff}' if name=='rushil' else '')
                env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE')
                final=case/f'{label}-final.txt'
                stdout=run([out/binary,input_path,args.dt,args.steps,args.samples,args.warmups,
                            args.theta,final,backend],case/f'{label}-samples.csv',env)
                samples=[float(row['seconds']) for row in csv.DictReader(stdout.splitlines())]
                if len(samples)!=args.samples or any(not math.isfinite(s) or s<=0 for s in samples):
                    raise ValueError('Invalid timing samples')
                state=read_state(final)
                if reference is None: reference=state
                if name=='rushil':
                    if not cutoff_reference:
                        cutoff_reference=state
                    elif state != cutoff_reference:
                        raise ValueError('Task cutoffs/thread counts changed final state')
                mean=statistics.mean(samples)/args.steps*1000
                if name=='src': baseline=mean
                if name=='n2': n2_baseline=mean
                record=dict(bodies=n,implementation=name,threads=threads,summarize_base_case=cutoff,steps=args.steps,
                            dt=args.dt,theta=args.theta,samples=args.samples,mean_ms_per_step=mean,
                            stddev_ms_per_step=statistics.stdev(samples)/args.steps*1000,
                            src_speedup=(baseline/mean if baseline is not None else ''),
                            n2_speedup=n2_baseline/mean,
                            **errors(state,reference))
                summaries.append(record)
                with (out/'comparison.csv').open('w',newline='') as f:
                    writer=csv.DictWriter(f,fieldnames=list(summaries[0]))
                    writer.writeheader(); writer.writerows(summaries)
                print(f'{n:5d} {label:18s} {mean:9.4f} ms/step; position RMS vs n2={record["position_rms"]:.3g}',flush=True)
        with (out/'comparison.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(summaries[0])); writer.writeheader(); writer.writerows(summaries)
        status['status']='success'
    except Exception as exc:
        status.update(status='failed',error=str(exc))
        raise
    finally:
        status_path.write_text(json.dumps(status,indent=2))


if __name__=='__main__':
    main()
