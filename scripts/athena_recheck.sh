#!/bin/bash
#SBATCH --job-name=findrna-recheck
#SBATCH --partition=hmq
#SBATCH --exclude=athena-h03
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
reference="${ATHENA_REFERENCE_ROOT:-$HOME/findrna-20261115}"
export PATH="$HOME/julia-1.12.7/bin:$reference/tools/bin:$HOME/miniforge3/envs/barrnap/bin:$PATH"
export BARRNAP="$reference/tools/barrnap/bin/barrnap"
export NHMMER="$reference/tools/bin/nhmmer"
export JULIA_NUM_THREADS=4 BENCH_THREADS=4 BENCH_REPEATS=3
julia --startup-file=no --project=. scripts/develop.jl
julia --startup-file=no --project=. scripts/test.jl
export FINDRNA_RESULTS="$PWD/benchmark/results/synthetic"
julia --startup-file=no --project=. benchmark/compare.jl
failed=0
for genome in "$reference"/genomes/*.fasta; do
    export FINDRNA_RESULTS="$PWD/benchmark/results/$(basename "$genome" .fasta)"
    julia --startup-file=no --project=. benchmark/compare.jl "$genome" || failed=1
done
exit "$failed"
