using Test, FindRNA, ProfileHMMs, BioSequences, Tables
include("fixtures.jl")

@testset "FindRNA" begin
    ps = parse_hmm_file(default_profiles())
    @test getfield.(ps,:name) == ["16S_rRNA","23S_rRNA","5S_rRNA"]
    @test getfield.(ps,:length) == [1533,2906,119]
    mktempdir() do dir
        path = joinpath(dir,"fixture.fna")
        _, seqs = make_fixture(path)
        r = findrna(path; progress=false)
        @test r isa RNAResult
        @test length(r) == 8
        @test all(h -> h.target_name ∉ ("negative","ambiguous"),r)
        for p in ps, (suffix,strand) in (("forward",'+'),("reverse",'-'))
            hits = filter(h -> h.target_name == p.name * "_" * suffix && h.query_name == p.name,r)
            @test length(hits) == 1
            if length(hits) == 1
                @test hits[1].strand == strand
                @test abs(hits[1].target_start - 81) <= 5
                @test abs(hits[1].target_end - (80+p.length)) <= 5
            end
        end
        @test length(filter(h -> h.target_name == "two_16S",r)) == 2
        @test Tables.istable(typeof(r))
        @test length(Tables.columntable(r).target_name) == length(r)
        extracted = extract_sequences(r)
        @test length(extracted) == length(r)
        @test all(length(s) == h.target_end-h.target_start+1 for (s,h) in zip(extracted,r))
        @test !occursin('\e',sprint(show,MIME"text/plain"(),r))
        write_gff(joinpath(dir,"hits.gff"),r)
        @test startswith(read(joinpath(dir,"hits.gff"),String),"##gff-version 3\n")
        write_fasta(joinpath(dir,"hits.fa"),r)
        @test count(startswith(">"),readlines(joinpath(dir,"hits.fa"))) == length(r)
        @test_throws ArgumentError findrna(path; models=["no_such_model"],progress=false)
        @test_throws ArgumentError findrna(path; models="16S_rRNA",progress=false)
        empty = findrna(path; models=String[],progress=false)
        @test isempty(empty)
        @test isempty(extract_sequences(empty))
    end
end
