#!/bin/bash
#SBATCH --job-name=derotate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=0-03:00:00
#SBATCH --output=validation_output/vgt_irgetas_logs/%x-%j.out
#SBATCH --error=validation_output/vgt_irgetas_logs/%x-%j.err

# Shift cubes by the local disk velocity and write a VGT config per cube
# (rohsa/cli/derotate_cube.jl).
#   sbatch --export=ALL,ROHSA_REPO=$PWD,REFERENCE=<reference_velocity_channels.fits>,OUT_DIR=<out>,CONFIG_DIR=<dir with vgt_<name>.toml>,CUBES="name=cube.fits ..." hpc/irgetas/derotate.slurm

set -euo pipefail
: "${ROHSA_REPO:?}" "${REFERENCE:?}" "${OUT_DIR:?}" "${CONFIG_DIR:?}" "${CUBES:?}"
cd "$ROHSA_REPO"
mkdir -p validation_output/vgt_irgetas_logs
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$ROHSA_REPO/.julia-irgetas}"
threads="${SLURM_CPUS_PER_TASK:-1}"
echo "job_id=$SLURM_JOB_ID reference=$REFERENCE out=$OUT_DIR threads=$threads"
# shellcheck disable=SC2086
srun --ntasks=1 julia --project=. --threads="$threads" rohsa/cli/derotate_cube.jl "$REFERENCE" "$OUT_DIR" "$CONFIG_DIR" $CUBES
