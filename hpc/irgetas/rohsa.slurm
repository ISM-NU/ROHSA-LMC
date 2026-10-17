#!/bin/bash
#SBATCH --job-name=rohsa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --output=validation_output/rohsa_irgetas_logs/%x-%j.out
#SBATCH --error=validation_output/rohsa_irgetas_logs/%x-%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

: "${VGT_REPO:?submit.sh must set VGT_REPO}"
: "${ROHSA_CONFIG:?submit.sh must set ROHSA_CONFIG}"

cd "$VGT_REPO"
mkdir -p validation_output/rohsa_irgetas_logs
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$VGT_REPO/.julia-irgetas}"

echo "job_id=$SLURM_JOB_ID"
echo "host=$(hostname)"
echo "partition=$SLURM_JOB_PARTITION"
echo "config=$ROHSA_CONFIG"
echo "julia_depot=$JULIA_DEPOT_PATH"
date --iso-8601=seconds

# The checked-in executable may have been built on another computer. Rebuild
# from the unchanged bundled upstream source on Irgetas for a compatible binary.
module purge
module load GCC
make -C rohsa/upstream/ROHSA/src clean all

command -v julia >/dev/null || {
    echo "Julia is not on PATH. Follow rohsa/hpc/irgetas/README.md before submitting." >&2
    exit 2
}

srun --ntasks=1 julia --project=. rohsa/cli/preprocess.jl --config "$ROHSA_CONFIG"

date --iso-8601=seconds
