#!/bin/bash
#SBATCH --job-name=nbody-prepare
#SBATCH --partition=mit_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load gcc/12.2.0
export CC=gcc OPENMP=1 LTO=1
git rev-parse HEAD
gcc --version
python3 --version
exec python3 -u bench/engaging.py submit "$@"
