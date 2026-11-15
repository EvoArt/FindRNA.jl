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
        @test length(r) == 9
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

@testset "Short real RNA inside a long window" begin
    # NCBI CP148752.1:116377-116486, a 5S locus missed without null-length scoring.
    lines = readlines(joinpath(@__DIR__,"fixtures","pseudomonas-5s.fna"))
    s = LongDNA{4}(join(filter(x -> !startswith(x,">"),lines)))
    @test length(s) == 110
    h = build_plan7(parse_hmm_file(default_profiles())[3]; target_length=2000)
    padded = LongDNA{4}("N"^2400 * String(s) * "N"^2490)
    hits = search_sequence(h,padded)
    @test length(hits) == 1
    if length(hits) == 1
        @test hits[1].strand == '+'
        @test (hits[1].target_start,hits[1].target_end) == (2401,2510)
    end
end
