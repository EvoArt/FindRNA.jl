#!/bin/bash
# Download the accessions in benchmark/genomes.txt from NCBI into $1 (default genomes/).
set -euo pipefail
out="${1:-genomes}"
mkdir -p "$out"
grep -v '^#' "$(dirname "$0")/genomes.txt" | while read -r acc; do
    [ -n "$acc" ] || continue
    [ -s "$out/$acc.fasta" ] && continue
    curl -fsSL --retry 5 --retry-delay 5 \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text" \
        -o "$out/$acc.fasta"
    head -c 1 "$out/$acc.fasta" | grep -q '>'
    sleep 1
done
