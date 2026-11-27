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
crop=${crop:-0}

# The preflight resolves the candidates (an [auto_config] n_gauss axis can sweep
# several N) and reports the largest, which sets the memory request.
repo=$(cd "$(dirname "$0")/../.." && pwd)
mkdir -p validation_output/rohsa_irgetas_logs
preflight_log=$(mktemp)
trap 'rm -f "$preflight_log"' EXIT
julia --project="$repo" "$repo/hpc/irgetas/rohsa_preflight.jl" "$(realpath "$config")" | tee "$preflight_log"
n_gauss=$(sed -n 's/^ROHSA largest candidate N: \([0-9][0-9]*\)$/\1/p' "$preflight_log")
[ -n "$n_gauss" ] || { echo "preflight did not report the largest N" >&2; exit 1; }

partition="ZEN4"
if [ "$crop" -gt 0 ] && [ "$crop" -le 512 ]; then
    cpus=16; memory="12G"; time_limit="03:00:00"
elif [ "$crop" -gt 0 ] && [ "$crop" -le 1024 ]; then
    cpus=16; memory="40G"; time_limit="12:00:00"
elif [ "$crop" -gt 0 ] && [ "$crop" -le 2048 ]; then
    cpus=32
    # Memory grows with N, so N 7-8 must fit in what the measured N = 16 job used (89 GB, f34c7ec). N 7-8 itself was never measured.
    if [ "$n_gauss" -le 6 ]; then memory="100G"; time_limit="1-00:00:00"
    elif [ "$n_gauss" -le 8 ]; then memory="120G"; time_limit="1-12:00:00"
    else memory="120G"; time_limit="08:00:00"; fi  # N = 16 measured: 89 GB, under 4 h
else
    # Whole cube (crop_size = 0), e.g. the 3624^2 foreground cube.
    cpus=32; memory="180G"; time_limit="3-00:00:00"
fi
cpus=${CPUS:-$cpus}
memory=${MEM:-$memory}
time_limit=${TIME:-$time_limit}

# Paths inside the config (input, output_dir) resolve from the directory you
# submit from, so an experiment repo can keep its configs and data and call
# this script from its sibling checkout (../ROHSA-LMC).
run_dir=$PWD
config=$(realpath "$config")

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
