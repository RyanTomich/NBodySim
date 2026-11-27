#!/bin/bash
#SBATCH --job-name=nbody-sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --threads-per-core=1

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load gcc/12.2.0
export LC_ALL=C OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
unset GOMP_CPU_AFFINITY
hostname
lscpu
lscpu -e=CPU,CORE,SOCKET,NODE
exec srun --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    --distribution=block:block --cpu-bind=verbose,threads --mem-bind=local --cpu-freq=Performance \
    python3 -u bench/engaging.py worker "$1" "$SLURM_ARRAY_TASK_ID"
