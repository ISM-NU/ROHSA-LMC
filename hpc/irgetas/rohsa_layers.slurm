#!/bin/bash
#SBATCH --job-name=rohsa_layers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=0-06:00:00
#SBATCH --output=validation_output/vgt_irgetas_logs/%x-%j.out
#SBATCH --error=validation_output/vgt_irgetas_logs/%x-%j.err

# Build ROHSA velocity-layer cubes and their VGT configs (rohsa/cli/layer_cubes.jl).
#   sbatch --export=ALL,ROHSA_REPO=$PWD,REGION_DIR=<rohsa region dir>,N_GAUSS=<n>,LAYER_DIR=<out> hpc/irgetas/rohsa_layers.slurm

set -euo pipefail
: "${ROHSA_REPO:?}" "${REGION_DIR:?}" "${N_GAUSS:?}" "${LAYER_DIR:?}"
cd "$ROHSA_REPO"
mkdir -p validation_output/vgt_irgetas_logs
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$ROHSA_REPO/.julia-irgetas}"
threads="${SLURM_CPUS_PER_TASK:-1}"
echo "job_id=$SLURM_JOB_ID region=$REGION_DIR n_gauss=$N_GAUSS out=$LAYER_DIR threads=$threads"
srun --ntasks=1 julia --project=. --threads="$threads" rohsa/cli/layer_cubes.jl "$REGION_DIR" "$N_GAUSS" "$LAYER_DIR" ${LAYER_ARGS:-}
