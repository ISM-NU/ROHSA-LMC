#!/bin/bash
# run_one.sh <name> <binary> <input> <rms> <n_gauss> <threads>
set -e
cd "$(dirname "$0")"
name=$1; bin=$2; input=$3; rms=$4; n=$5; thr=$6
mkdir -p runs
sed -e "s#@INPUT@#$PWD/$input#" -e "s#@OUT@#$PWD/runs/$name.out#" -e "s#@RMS@#$PWD/$rms#" \
    -e "s#@N@#$n#" -e "s#@VARSIG@#0.0d0#" params_template.txt > runs/$name.params
# real inputs live outside this directory
sed -i "s#$PWD//#/#g" runs/$name.params
OMP_NUM_THREADS=$thr /usr/bin/time -f "$name wall %e s cpu %P maxrss %M KB" $bin runs/$name.params > runs/$name.log 2> runs/$name.time
tail -n 1 runs/$name.time
