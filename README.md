# FindRNA.jl

Find bacterial ribosomal RNA in FASTA files from a native Julia session on
Windows, macOS or Linux. Intended for small interactive analyses in VS Code
as well as scripts, with the display and progress conventions of DenoiseDNA
and AssignTaxonomy.

**Development version; not yet registered.** The search engine is the separate
Julia package [ProfileHMMs](packages/ProfileHMMs), with its own project, tests
and license. This repository contains both packages; neither needs a HMMER
executable at runtime.

## Try it locally

From this checkout, run `julia --project=. scripts/develop.jl` once. Then:

```julia
using FindRNA
hits = findrna("assembly.fna")
hits16s = findrna("assembly.fna"; models=["16S_rRNA"])
sequences = extract_sequences(hits16s)  # Vector{LongDNA{4}}, oriented to the RNA
write_fasta("16S.fasta", hits16s)
write_gff("rna.gff", hits)
FindRNA.setprogress!(false)             # or progress=false on a single call
```

Results support iteration, indexing and the Tables.jl interface:
`DataFrame(hits)` works when DataFrames is loaded. Pass the extracted sequences
to AssignTaxonomy's `assign_taxonomy(sequences, reference)` workflow; FindRNA
does not force a dependency on the classifier or DenoiseDNA. It operates on
DNA FASTA input, including sequences exported by other tools.

The bundled models cover bacterial 16S, 23S and 5S. Supply `profiles="custom.hmm"`
for another HMMER3 nucleotide profile database. Protein HMMs are rejected.
The name leaves room for other RNA profiles; this is not a tRNA predictor or
a replacement for modern Barrnap's complete annotation pipeline.

## Search behaviour and limits

Coordinates are 1-based, inclusive, on the original FASTA sequence. Both strands
are searched. Repeated non-overlapping copies are recovered by searching the
remaining flanks. Long sequences use overlapping windows; each window owns its
results, allowing parallel scans without indexing by a thread-pool ID.
Default window sizes scale with model length; custom windows must satisfy
`window > overlap >= 0`. Choose an overlap long enough to contain the expected
RNA. Overlapping calls from adjacent windows retain the best-scoring alignment.

The engine computes local single-hit Viterbi alignments and uses an ungapped
prefilter. It does **not** implement HMMER's Forward/domain pipeline, composition
bias corrections or Infernal covariance models. Reported E-values are approximate
Viterbi tail probabilities scaled by `db_size` (default 1), not calibrated
whole-database nhmmer E-values. Ambiguous bases have neutral emission odds.
Set `msv_threshold=Inf` to bypass the heuristic prefilter when investigating
sensitivity. Do not assume differences from HMMER are merely floating-point
rounding or implementation-language effects.

## Tests and comparisons

`julia --project=. --threads=2,1 scripts/test.jl` tests both packages. CI runs
Julia 1.10 and current stable on Windows, Linux and macOS. A separate Linux job
compares against HMMER 3.4 and barrnap 0.9 using identical nucleotide models.
It checks model/sequence identity, strand, one-to-one hit counts and endpoints
within 30 bases. This is an explicit regression tolerance for the deterministic
fixture, not a claim that all biological sequences agree within 30 bases.

See [benchmark instructions](benchmark/README.md). Reference comparison failures
are errors, never optional skips. Synthetic fixtures test predictable behaviour;
real-genome and modern Barrnap comparisons are needed before a release claims
biological parity. Timing reports distinguish first Julia call (including JIT),
warmed Julia calls, and reference subprocess wall time.

## Licensing

Julia code: MIT. Rfam model data: CC0, with source and checksum in
[data/README.md](data/README.md). Research data, downloaded tools, old experiments
and local manifests are excluded from version control.
