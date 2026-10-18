#!/bin/bash -l
#SBATCH --job-name=t_serial_paired
#SBATCH --partition=ZEN4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=timing/%x-%j.out
# Serial speed of master, #2, #3 and all four combined, measured fairly: the four
# builds run at the same time on four cores of one node (same load, same memory
# bandwidth), and each round rotates which core each build gets. Three rounds.
# Usage: sbatch serial_paired_job.sh <input.dat> <rms.dat> <n_gauss> <var_sig>
set -e
module load GCC
input=$1; rms=$2; n=$3; varsig=$4
builds=(master pr2 pr3 all)
cores=($(taskset -cp $$ | sed 's/.*: //' | tr ',' '\n' | while IFS=- read a b; do seq $a ${b:-$a}; done | head -4))
echo "node=$(hostname) cores=${cores[*]} cpu=$(lscpu | sed -n 's/^Model name: *//p')"
mkdir -p timing/paired
for round in 1 2 3; do
  pids=()
  for i in 0 1 2 3; do
    b=${builds[$i]}; core=${cores[$(( (i + round) % 4 ))]}; name=paired/${b}_r$round
    sed -e "s#@INPUT@#$input#" -e "s#@OUT@#$PWD/timing/$name.out#" -e "s#@RMS@#$rms#" \
        -e "s#@N@#$n#" -e "s#@VARSIG@#$varsig#" params_template.txt > timing/$name.params
    ( t0=$(date +%s.%N)
      OMP_NUM_THREADS=1 taskset -c $core src/$b/ROHSA_serial timing/$name.params > timing/$name.log 2>&1
      t1=$(date +%s.%N)
      echo "$b round $round core $core wall $(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.1f", b-a}') s" | tee timing/$name.time ) &
    pids+=($!)
  done
  wait "${pids[@]}"
done
md5sum timing/paired/*.out
