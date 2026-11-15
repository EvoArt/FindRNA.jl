using FindRNA
include(joinpath(@__DIR__, "..", "benchmark", "compare.jl"))
mktempdir() do dir
    path = joinpath(dir,"fixture.fna")
    make_fixture(path)
    r = findrna(path;progress=false)
    @assert length(r) == 9
    write_gff(joinpath(dir,"out.gff"),r)
    @assert length(parse_gff(joinpath(dir,"out.gff"))) == 9
    @testset "Comparison round trip" begin
        compare_hits(r, parse_gff(joinpath(dir,"out.gff")))
    end
end
