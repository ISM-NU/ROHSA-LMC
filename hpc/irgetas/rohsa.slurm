#!/bin/bash
#SBATCH --job-name=rohsa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=validation_output/rohsa_irgetas_logs/%x-%j.out
#SBATCH --error=validation_output/rohsa_irgetas_logs/%x-%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

: "${VGT_REPO:?submit_rohsa.sh must set VGT_REPO}"
: "${ROHSA_CONFIG:?submit_rohsa.sh must set ROHSA_CONFIG}"

cd "$VGT_REPO"
mkdir -p validation_output/rohsa_irgetas_logs
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$VGT_REPO/.julia-irgetas}"

echo "job_id=$SLURM_JOB_ID"
echo "host=$(hostname)"
echo "partition=$SLURM_JOB_PARTITION"
echo "config=$ROHSA_CONFIG"
echo "julia_depot=$JULIA_DEPOT_PATH"
date --iso-8601=seconds

# Build the bundled upstream source on Irgetas, but only when the binary is
# missing or older than a source file. Concurrent jobs share this source tree,
# and rebuilding replaces the binary while another job may be starting it.
# After copying a binary from another machine, delete it so it is rebuilt here.
module purge
module load GCC
rohsa_src=rohsa/upstream/ROHSA/src
flock "$rohsa_src/.build.lock" bash -c '
    if [ ! -x "$1/ROHSA" ] || [ -n "$(find "$1" -maxdepth 1 \( -name "*.f" -o -name "*.f90" -o -name Makefile \) -newer "$1/ROHSA")" ]; then
        make -C "$1" clean all
    fi' _ "$rohsa_src"

command -v julia >/dev/null || {
    echo "Julia is not on PATH. Follow hpc/irgetas/README.md before submitting." >&2
    exit 2
}

# ROHSA's cost/gradient loop is OpenMP-parallel. Match threads to the Slurm
# allocation; an unset OMP_NUM_THREADS would start one thread per node core.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
# Headroom for automatic arrays on the stack under -fopenmp. Not the cause of
# the earlier segfault at "Start updating last level", which was the 32-bit
# L-BFGS-B workspace overflow (fixed, docs/upstream-prs/rohsa/01-*.md).
ulimit -s unlimited
export OMP_STACKSIZE=256M
echo "omp_threads=$OMP_NUM_THREADS"

srun --ntasks=1 julia --project=. rohsa/cli/preprocess.jl --config "$ROHSA_CONFIG"

date --iso-8601=seconds
