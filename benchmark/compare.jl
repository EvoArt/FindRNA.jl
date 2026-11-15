# Required references are never silently skipped. Output files retain raw evidence.
using FindRNA, ProfileHMMs, Test, Statistics, SHA, TOML
include(joinpath(@__DIR__, "..", "test", "fixtures.jl"))

function parse_nhmmer(path)
    rows = NamedTuple[]
    for line in eachline(path)
        (startswith(line,"#") || isempty(strip(line))) && continue
        p = split(line)
        a,b = parse.(Int,p[7:8])
        push!(rows,(target=String(p[1]), model=String(p[3]), strand=only(p[12]),
                    start=min(a,b), stop=max(a,b), bits=parse(Float64,p[14])))
    end
    rows
end

function parse_gff(path)
    rows = NamedTuple[]
    for line in eachline(path)
        (startswith(line,"#") || isempty(strip(line))) && continue
        p = split(line,'\t')
        p[3] == "rRNA" || continue
        m = match(r"(?:^|;)Name=([^;]+)",p[9])
        m === nothing && error("GFF record lacks model Name: $line")
        push!(rows,(target=String(p[1]),model=String(m[1]),strand=only(p[7]),
                    start=parse(Int,p[4]),stop=parse(Int,p[5]),bits=NaN))
    end
    rows
end

function compare_hits(got, ref; coordinate_tolerance=30)
    # One-to-one matching: an extra call fails just as a missing call does.
    unmatched = collect(eachindex(got))
    deltas = Int[]
    for r in ref
        candidates = filter(unmatched) do i
            h = got[i]
            h.target_name == r.target && h.query_name == r.model && h.strand == r.strand
        end
        @test !isempty(candidates)
        isempty(candidates) && continue
        errors = [max(abs(got[i].target_start-r.start),abs(got[i].target_end-r.stop)) for i in candidates]
        k = argmin(errors); i = candidates[k]
        @test errors[k] <= coordinate_tolerance
        push!(deltas,errors[k]); filter!(!=(i),unmatched)
    end
    @test isempty(unmatched)
    deltas
end

function benchmark_main()
    out = abspath(get(ENV,"FINDRNA_RESULTS",joinpath(@__DIR__,"results")))
    mkpath(out)
    fasta = isempty(ARGS) ? joinpath(out,"fixture.fna") : abspath(ARGS[1])
    isempty(ARGS) && make_fixture(fasta)
    hmmer = get(ENV,"NHMMER","nhmmer")
    barrnap = get(ENV,"BARRNAP","barrnap")
    threads = parse(Int,get(ENV,"BENCH_THREADS","1"))
    repeats = parse(Int,get(ENV,"BENCH_REPEATS","3"))
    repeats > 0 || error("BENCH_REPEATS must be positive")
    hv = read(`$hmmer -h`,String)
    occursin("HMMER 3.4",hv) || error("HMMER 3.4 required; got: $hv")
    bv = read(pipeline(`$barrnap --version`,stderr=stdout),String)
    write(joinpath(out,"versions.txt"),"Julia $VERSION\n$hv\n$bv")
    tbl = joinpath(out,"nhmmer.tbl")
    gff = joinpath(out,"barrnap.gff")
    hmm = default_profiles()
    # Equal input and profiles for HMMER/0.9. Modern Barrnap uses its own CMs.
    nhcmd = `$hmmer --cpu $threads --noali -E 1e-5 --tblout $tbl $hmm $fasta`
    bacmd = `$barrnap --threads $threads --kingdom bac --evalue 1e-5 $fasta`
    function external(cmd, output, err)
        run(pipeline(cmd; stdout=output, stderr=err))
    end
    opts = (progress=false,threaded=threads>1)
    cold = @elapsed hits = findrna(fasta;opts...)
    jt = [@elapsed(findrna(fasta;opts...)) for _ in 1:repeats]
    nt = [@elapsed(external(nhcmd,joinpath(out,"nhmmer.txt"),joinpath(out,"nhmmer.stderr"))) for _ in 1:repeats]
    bt = [@elapsed(external(bacmd,gff,joinpath(out,"barrnap.stderr"))) for _ in 1:repeats]
    write_gff(joinpath(out,"findrna.gff"),hits)
    refs = parse_nhmmer(tbl); br = parse_gff(gff)
    report = Dict("julia"=>string(VERSION),"threads"=>threads,"repeats"=>repeats,
        "input_sha256"=>bytes2hex(sha256(read(fasta))),"model_sha256"=>bytes2hex(sha256(read(hmm))),
        "findrna_first_call_seconds"=>cold,"findrna_warm_seconds"=>jt,
        "nhmmer_process_seconds"=>nt,"barrnap_process_seconds"=>bt,
        "findrna_hits"=>length(hits),"nhmmer_hits"=>length(refs),"barrnap_hits"=>length(br))
    open(joinpath(out,"timings.toml"),"w") do io; TOML.print(io,report); end
    println("Median seconds: FindRNA warm=",median(jt),"; nhmmer=",median(nt),"; barrnap=",median(bt))
    if get(ENV,"BARRNAP_MODERN","false") == "true"
        # Different models/algorithms: report all raw calls; no fabricated parity.
        @testset "HMMER reference" begin
            @test !isempty(refs)
            compare_hits(hits,refs)
        end
        println("Modern Barrnap is reported separately; CM/HMM parity is not asserted.")
    else
        occursin("0.9",bv) || error("Expected barrnap 0.9; set BARRNAP_MODERN=true for CM comparison")
        @testset "Reference hit identity, strand and endpoints" begin
            @test !isempty(refs)
            @test !isempty(br)
            compare_hits(hits,refs)
            compare_hits(hits,br)
        end
    end
end
abspath(PROGRAM_FILE) == (@__FILE__) && benchmark_main()
