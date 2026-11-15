#!/bin/bash
#SBATCH --job-name=findrna
#SBATCH --partition=hmq
#SBATCH --exclude=athena-h03
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
# Submit from the prepared repository after loading Julia, nhmmer 3.4 and
# barrnap 0.9 (or set BARRNAP_MODERN=true and BARRNAP to modern barrnap).
cd "${SLURM_SUBMIT_DIR:?}"
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export BENCH_THREADS="${SLURM_CPUS_PER_TASK}"
export FINDRNA_RESULTS="benchmark/results/athena-${SLURM_JOB_ID}"
julia --startup-file=no --project=. benchmark/compare.jl "$@"
