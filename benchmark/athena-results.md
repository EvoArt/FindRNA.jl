# Athena results (job 4052391, 2026-11-15)

Skylake-AVX512 node, 4 threads; Julia 1.12.7; HMMER 3.4; barrnap 0.9 and 1.10.6.
All tools searched with E-value 1e-5; nhmmer and barrnap 0.9 used the same profiles
as FindRNA. Times are medians of 3 runs (wall clock, whole process for nhmmer
and barrnap). Counts are reference calls.

Rows marked `modern` compare against barrnap 1.10.6 (Infernal covariance
models); its 23S calls are absent for three genomes and its 5S boundaries
differ by up to 12 bp, so these rows are reported, not asserted.

| run | FindRNA first / warm (s) | nhmmer (s) | barrnap (s) | vs nhmmer: matched / endpoints >10 bp / missing / extra | vs barrnap |
|---|---|---|---|---|---|
| 0.9-CP148752.1 | 2.22 / 1.96 | 4.49 | 4.16 | 25 / 0 / 0 / 0 (≤2 bp) | 25 / 0 / 0 / 0 (≤2 bp) |
| 0.9-CP148753.1 | 1.81 / 1.64 | 3.65 | 3.44 | 13 / 0 / 0 / 0 (≤3 bp) | 13 / 0 / 0 / 0 (≤3 bp) |
| 0.9-JBBHKO010000001.1 | 1.64 / 1.36 | 2.85 | 2.85 | 9 / 0 / 0 / 0 (≤3 bp) | 9 / 0 / 0 / 0 (≤3 bp) |
| 0.9-NZ_JBBHKP010000001.1 | 1.33 / 1.0 | 2.12 | 2.11 | 10 / 0 / 0 / 0 (≤6 bp) | 10 / 0 / 0 / 0 (≤5 bp) |
| 0.9-NZ_JBBHKQ010000001.1 | 1.3 / 1.08 | 2.56 | 2.29 | 12 / 0 / 0 / 0 (≤4 bp) | 12 / 0 / 0 / 0 (≤4 bp) |
| 0.9-synthetic | 0.37 / 0.23 | 0.72 | 0.89 | 9 / 0 / 0 / 0 (≤3 bp) | 9 / 0 / 0 / 0 (≤3 bp) |
| modern-CP148752.1 | 1.94 / 1.86 | 4.53 | 16.25 | 25 / 0 / 0 / 0 (≤2 bp) | 25 / 0 / 0 / 0 (≤4 bp) |
| modern-CP148753.1 | 1.73 / 1.6 | 3.6 | 10.57 | 13 / 0 / 0 / 0 (≤3 bp) | 9 / 0 / 0 / 4 (≤4 bp) |
| modern-JBBHKO010000001.1 | 1.54 / 1.35 | 2.84 | 10.44 | 9 / 0 / 0 / 0 (≤3 bp) | 6 / 0 / 0 / 3 (≤7 bp) |
| modern-NZ_JBBHKP010000001.1 | 1.14 / 1.1 | 2.14 | 9.08 | 10 / 0 / 0 / 0 (≤6 bp) | 6 / 4 / 0 / 0 (≤12 bp) |
| modern-NZ_JBBHKQ010000001.1 | 1.21 / 1.05 | 2.56 | 10.43 | 12 / 0 / 0 / 0 (≤4 bp) | 8 / 0 / 0 / 4 (≤4 bp) |
| modern-synthetic | 0.36 / 0.23 | 0.72 | 27.13 | 9 / 0 / 0 / 0 (≤3 bp) | 9 / 0 / 0 / 0 (≤4 bp) |
