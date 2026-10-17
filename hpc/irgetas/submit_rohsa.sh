#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 pilot|full"
}

mode=${1:-}
case "$mode" in
    pilot)
        config="rohsa/hpc/irgetas/config.pilot.toml"
        partition="ZEN4"
        memory="64G"
        time_limit="1-00:00:00"
        ;;
    full)
        config="rohsa/config.toml"
        partition="ZEN4"
        memory="180G"
        time_limit="7-00:00:00"
        ;;
    *)
        usage >&2
        exit 2
        ;;
esac

repo=$(cd "$(dirname "$0")/../../.." && pwd)
cd "$repo"
mkdir -p validation_output/rohsa_irgetas_logs

julia --project=. rohsa/hpc/irgetas/preflight.jl "$config"

sbatch \
    --partition="$partition" \
    --mem="$memory" \
    --time="$time_limit" \
    --export=ALL,VGT_REPO="$repo",ROHSA_CONFIG="$config" \
    rohsa/hpc/irgetas/run.slurm
