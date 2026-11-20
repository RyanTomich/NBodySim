#!/bin/bash
#SBATCH --job-name=nbody-sweep
#SBATCH --partition=mit_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
echo "After completion: sacct -j $SLURM_JOB_ID --units=G -o JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS"
module load gcc/12.2.0
export CC=gcc OPENMP=1 LTO=1
export OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
hostname
git rev-parse HEAD
gcc --version
python3 --version
lscpu
lscpu -e=CPU,CORE,SOCKET,NODE

exec srun --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    --distribution=block:block --hint=nomultithread --cpu-bind=verbose --cpu-freq=Performance \
    python3 -u bench/compare_solvers.py --engaging --check-placement "$@"
