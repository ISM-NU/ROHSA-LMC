#!/bin/bash -l
#SBATCH --partition=ZEN4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --time=16:00:00
#SBATCH --output=timing/%x-%j.out
# One timed ROHSA run: sbatch --job-name=<name> --cpus-per-task=<threads> timing_job.sh <name> <binary> <input> <rms> <n_gauss> <var_sig>
# Each run gets its own allocation with exactly <threads> cores, so runs do not share cores.
set -e
module load GCC
name=$1; bin=$2; input=$3; rms=$4; n=$5; varsig=$6
thr=${SLURM_CPUS_PER_TASK:-1}
mkdir -p timing
sed -e "s#@INPUT@#$input#" -e "s#@OUT@#$PWD/timing/$name.out#" -e "s#@RMS@#$rms#" \
    -e "s#@N@#$n#" -e "s#@VARSIG@#$varsig#" params_template.txt > timing/$name.params
echo "name=$name bin=$bin threads=$thr node=$(hostname) cpu=$(lscpu | sed -n 's/^Model name: *//p')"
t0=$(date +%s.%N)
OMP_NUM_THREADS=$thr OMP_PROC_BIND=close OMP_PLACES=cores $bin timing/$name.params > timing/$name.log 2>&1
t1=$(date +%s.%N)
# compute nodes have no /usr/bin/time; peak memory is in sacct (MaxRSS)
echo "$name threads $thr wall $(awk -v a=$t0 -v b=$t1 "BEGIN{printf \"%.1f\", b-a}") s" | tee timing/$name.time
