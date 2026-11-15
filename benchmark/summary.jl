# Markdown summary of the result directories written by compare.jl:
#   julia benchmark/summary.jl benchmark/results/<run>
using TOML

med(xs) = (s = sort(Float64.(xs)); n = length(s); isodd(n) ? s[cld(n, 2)] : (s[n÷2] + s[n÷2+1]) / 2)

function tally(path)
    isfile(path) || return "–"
    st = [split(l, '\t')[end] for l in readlines(path)[2:end]]
    e = [parse(Int, split(l, '\t')[8]) for l in readlines(path)[2:end] if split(l, '\t')[end] == "matched"]
    "$(count(==("matched"), st)) / $(count(==("missing"), st)) / $(count(==("extra"), st))" *
        (isempty(e) ? "" : " (≤$(maximum(e)) bp)")
end

root = isempty(ARGS) ? joinpath(@__DIR__, "results") : ARGS[1]
println("| run | FindRNA first / warm (s) | nhmmer (s) | barrnap (s) | vs nhmmer: matched / missing / extra | vs barrnap |")
println("|---|---|---|---|---|---|")
for d in sort(readdir(root; join=true))
    f = joinpath(d, "timings.toml"); isfile(f) || continue
    t = TOML.parsefile(f)
    r(x) = string(round(x; digits=2))
    println("| ", basename(d), " | ", r(t["findrna_first_call_seconds"]), " / ", r(med(t["findrna_warm_seconds"])),
            " | ", r(med(t["nhmmer_process_seconds"])), " | ", r(med(t["barrnap_process_seconds"])),
            " | ", tally(joinpath(d, "nhmmer-differences.tsv")), " | ", tally(joinpath(d, "barrnap-differences.tsv")), " |")
end
