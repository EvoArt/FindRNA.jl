#!/bin/bash
#SBATCH --job-name=findrna-validate
#SBATCH --partition=hmq
#SBATCH --exclude=athena-h03
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export PATH="$HOME/julia-1.12.7/bin:$PWD/tools/bin:$HOME/miniforge3/envs/barrnap/bin:$PATH"
export JULIA_NUM_THREADS=4
export BENCH_THREADS=4
export BENCH_REPEATS=3
mkdir -p tools build benchmark/results
if [ ! -x tools/bin/nhmmer ]; then
    curl -fL --retry 3 https://eddylab.org/software/hmmer/hmmer-3.4.tar.gz -o hmmer.tar.gz
    echo 'ca70d94fd0cf271bd7063423aabb116d42de533117343a9b27a65c17ff06fbf3  hmmer.tar.gz' | sha256sum -c -
    tar -xzf hmmer.tar.gz -C build
    (cd build/hmmer-3.4 && ./configure --prefix="$SLURM_SUBMIT_DIR/tools" && make -j4 && make install)
fi
if [ ! -d tools/barrnap ]; then
    git clone --depth 1 --branch 0.9 https://github.com/tseemann/barrnap.git tools/barrnap
fi
export BARRNAP="$PWD/tools/barrnap/bin/barrnap"
export NHMMER="$PWD/tools/bin/nhmmer"
julia --startup-file=no --project=. scripts/develop.jl
julia --startup-file=no --project=. scripts/test.jl
export FINDRNA_RESULTS="$PWD/benchmark/results/synthetic"
julia --startup-file=no --project=. benchmark/compare.jl
# Preserve each failure and finish the panel so discrepancies are inspectable.
# Download public accessions directly; never upload local research genomes.
mkdir -p genomes
for accession in CP148753.1 NZ_JBBHKQ010000001.1 CP148752.1 JBBHKO010000001.1 NZ_JBBHKP010000001.1; do
    curl -fL --retry 3 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$accession&rettype=fasta&retmode=text" -o "genomes/$accession.fasta"
    head -c 1 "genomes/$accession.fasta" | grep -q '>'
    sleep 1
done
failed=0
for genome in genomes/*.fasta; do
    export FINDRNA_RESULTS="$PWD/benchmark/results/$(basename "$genome" .fasta)"
    julia --startup-file=no --project=. benchmark/compare.jl "$genome" || failed=1
done
echo "REAL_GENOME_PARITY_FAILURE=$failed"
exit "$failed"
