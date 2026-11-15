# Registration

ProfileHMMs (github.com/EvoArt/ProfileHMMs.jl) must be in General before
FindRNA, which depends on it.

1. Register ProfileHMMs: comment `@JuliaRegistrator register` on the tested
   commit of ProfileHMMs.jl. New packages wait three days in General.
2. After ProfileHMMs is merged, check that this repository's CI passes against
   the registered version (`scripts/develop.jl` is then unnecessary; delete it
   and run `Pkg.test()` directly), then comment `@JuliaRegistrator register` on
   that FindRNA commit.
3. AutoMerge will flag the name FindRNA as close to the registered FoldRNA
   (Damerau–Levenshtein distance 2). Reply on the General PR explaining the
   difference (rRNA gene finding vs RNA folding) and ask for a manual merge.
4. TagBot tags each release once its registry PR is merged.
