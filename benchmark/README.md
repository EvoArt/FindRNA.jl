# Reproducible comparisons

Run with HMMER **3.4** and **barrnap 0.9** on PATH:

```sh
julia --project=. scripts/develop.jl
JULIA_NUM_THREADS=1 BENCH_THREADS=1 julia --project=. benchmark/compare.jl
```

The default deterministic fixture contains model-derived positive sequences on
both strands, two copies in one contig, a repetitive negative and all-N input.
Raw FASTA, GFF, nhmmer table, tool versions, SHA-256 hashes and repeated timings
are saved under `benchmark/results/`. Supply a real genome as the first argument.
Strict parity checks on real genomes may fail; retain the evidence and investigate
rather than widening tolerances to force success. False positives, missing hits
and strand mismatches fail, as do endpoints differing by more than 30 bases.

Use `NHMMER`, `BARRNAP`, `BENCH_THREADS`, `BENCH_REPEATS`, and `FINDRNA_RESULTS`
to select executables and settings. Use a fresh output directory for each run.
Julia's first-call time includes compilation; warm Julia and subprocess timings
are reported separately. Subprocess timings include startup and output writing.
This is an end-to-end benchmark, not a claim about isolated DP kernel speed.

For **modern Barrnap 1.10.6** (the current Bioconda recipe/tag), use a separate installation and set
`BARRNAP_MODERN=true`. Its covariance models and different pipeline mean raw
calls and timings are reported separately; HMMER parity remains checked.
Modern and 0.9 runs must use different result directories.

## Athena

`athena.sbatch` uses the `hmq` Slurm partition and excludes `athena-h03`, following
the existing TikzMakie Athena configuration. Adjust for the current cluster.
Stage this repository, resolve dependencies on the login node, load the reference
executables, and submit `sbatch benchmark/athena.sbatch [genome.fna]`.
The script runs on a compute node and records the job ID in the output directory.
The first Athena run is staged in `~/findrna-20261115`; results are retained per input.

HMMER 3.4: https://hmmer.org/archive.html

Barrnap: https://github.com/tseemann/barrnap/releases (latest labelled release 1.10.5;
the newer 1.10.6 tag is packaged by Bioconda).
