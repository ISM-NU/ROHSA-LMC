#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 pilot|pilot1024|crop2048|crop2048_noise|full"
}

mode=${1:-}
case "$mode" in
    pilot)
        config="configs/rohsa/lmc512_pilot.toml"
        partition="ZEN4"
        memory="64G"
        time_limit="1-00:00:00"
        ;;
    pilot1024)
        config="configs/rohsa/lmc1024_pilot.toml"
        partition="ZEN4"
        memory="64G"
        time_limit="1-00:00:00"
        ;;
    crop2048)
        config="configs/rohsa/lmc2048_n6_defaultnoise.toml"
        partition="ZEN4"
        memory="100G"
        time_limit="7-00:00:00"
        ;;
    crop2048_noise)
        config="configs/rohsa/lmc2048_n6_oldnoiseestimate.toml"
        partition="ZEN4"
        memory="100G"
        time_limit="7-00:00:00"
        ;;
    full)
        config="configs/rohsa/full_cube.toml"
        partition="ZEN4"
        memory="180G"
        time_limit="7-00:00:00"
        ;;
    *)
        usage >&2
        exit 2
        ;;
esac

repo=$(cd "$(dirname "$0")/../.." && pwd)
cd "$repo"
mkdir -p validation_output/rohsa_irgetas_logs

julia --project=. hpc/irgetas/rohsa_preflight.jl "$config"

sbatch \
    --partition="$partition" \
    --mem="$memory" \
    --time="$time_limit" \
    --export=ALL,VGT_REPO="$repo",ROHSA_CONFIG="$config" \
    hpc/irgetas/rohsa.slurm
