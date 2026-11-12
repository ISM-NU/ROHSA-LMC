#!/bin/bash
# Submit one ROHSA run config to Irgetas, sized from its crop and N.
#
#   bash hpc/irgetas/submit_rohsa.sh configs/rohsa/<name>.toml
#
# Override the sizing with CPUS=, MEM= or TIME= in the environment. The sizes
# follow the measured table in hpc/irgetas/README.md.
set -euo pipefail

config=${1:-}
if [ -z "$config" ] || [ ! -f "$config" ]; then
    echo "Usage: $0 configs/rohsa/<name>.toml" >&2
    exit 2
fi

toml_int() {
    # First integer value of key $1 in the config ("n_gauss = [8]" -> 8).
    sed -n "s/^$1 *= *\[\{0,1\}\([0-9][0-9]*\).*/\1/p" "$config" | head -n 1
}
crop=$(toml_int crop_size)
n_gauss=$(toml_int n_gauss)
crop=${crop:-0}
n_gauss=${n_gauss:-6}

partition="ZEN4"
if [ "$crop" -gt 0 ] && [ "$crop" -le 512 ]; then
    cpus=16; memory="12G"; time_limit="03:00:00"
elif [ "$crop" -gt 0 ] && [ "$crop" -le 1024 ]; then
    cpus=16; memory="40G"; time_limit="12:00:00"
elif [ "$crop" -gt 0 ] && [ "$crop" -le 2048 ]; then
    cpus=32
    if [ "$n_gauss" -le 6 ]; then memory="100G"; time_limit="1-00:00:00"
    else memory="140G"; time_limit="1-12:00:00"; fi
else
    # Whole cube (crop_size = 0), e.g. the 3624^2 foreground cube.
    cpus=32; memory="180G"; time_limit="3-00:00:00"
fi
cpus=${CPUS:-$cpus}
memory=${MEM:-$memory}
time_limit=${TIME:-$time_limit}

# Paths inside the config (input, output_dir) resolve from the directory you
# submit from, so an experiment repo can keep its configs and data and call
# this script through its deps/ROHSA-LMC submodule.
repo=$(cd "$(dirname "$0")/../.." && pwd)
run_dir=$PWD
config=$(realpath "$config")
mkdir -p validation_output/rohsa_irgetas_logs

julia --project="$repo" "$repo/hpc/irgetas/rohsa_preflight.jl" "$config"

name=$(basename "$config" .toml)
echo "submitting $name: crop=$crop N=$n_gauss cpus=$cpus mem=$memory time=$time_limit"
sbatch \
    --job-name="$name" \
    --partition="$partition" \
    --cpus-per-task="$cpus" \
    --mem="$memory" \
    --time="$time_limit" \
    --export=ALL,ROHSA_REPO="$repo",ROHSA_CONFIG="$config",RUN_DIR="$run_dir" \
    "$repo/hpc/irgetas/rohsa.slurm"
