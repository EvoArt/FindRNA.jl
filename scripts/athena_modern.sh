#!/bin/bash
#SBATCH --job-name=findrna-modern
#SBATCH --partition=hmq
#SBATCH --exclude=athena-h03
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
if [ ! -x tools/modern/bin/barrnap ]; then
    "$HOME/miniforge3/bin/conda" create -y --override-channels -c conda-forge -c bioconda \
        -p "$PWD/tools/modern" barrnap=1.10.6 hmmer=3.4
fi
export PATH="$HOME/julia-1.12.7/bin:$PWD/tools/modern/bin:$PATH"
export BARRNAP="$PWD/tools/modern/bin/barrnap"
reference="${ATHENA_REFERENCE_ROOT:-$HOME/findrna-20261115}"
export NHMMER="$reference/tools/bin/nhmmer"
export BARRNAP_MODERN=true
export JULIA_NUM_THREADS=4 BENCH_THREADS=4 BENCH_REPEATS=3
export FINDRNA_RESULTS="$PWD/benchmark/results/modern-synthetic"
julia --startup-file=no --project=. benchmark/compare.jl
failed=0
for genome in "$reference"/genomes/*.fasta; do
    [ -f "$genome" ] || continue
    export FINDRNA_RESULTS="$PWD/benchmark/results/modern-$(basename "$genome" .fasta)"
    julia --startup-file=no --project=. benchmark/compare.jl "$genome" || failed=1
done
exit "$failed"
