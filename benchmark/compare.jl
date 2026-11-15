# Required references are never silently skipped. Output files retain raw evidence.
using FindRNA, ProfileHMMs, Test, SHA, TOML
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

function write_differences(path, got, refs)
    remaining = collect(eachindex(got))
    open(path,"w") do io
        println(io,"target\tmodel\tstrand\treference_start\treference_end\tjulia_start\tjulia_end\tendpoint_error\tjulia_bits\treference_bits\tstatus")
        for r in refs
            candidates = filter(i -> (got[i].target_name,got[i].query_name,got[i].strand) ==
                                     (r.target,r.model,r.strand),remaining)
            if isempty(candidates)
                println(io,join((r.target,r.model,r.strand,r.start,r.stop,"","","","",r.bits,"missing"),'\t'))
                continue
            end
            errs = [max(abs(got[i].target_start-r.start),abs(got[i].target_end-r.stop)) for i in candidates]
            k = argmin(errs); i = candidates[k]; h = got[i]
            filter!(!=(i),remaining)
            println(io,join((r.target,r.model,r.strand,r.start,r.stop,h.target_start,h.target_end,
                            errs[k],h.bit_score,r.bits,errs[k]<=30 ? "matched" : "coordinates"),'\t'))
        end
        for i in remaining
            h = got[i]
            println(io,join((h.target_name,h.query_name,h.strand,"","",h.target_start,h.target_end,
                            "",h.bit_score,"","extra"),'\t'))
        end
    end
end

function median_time(xs)
    s = sort(xs); n = length(s)
    isodd(n) ? s[cld(n,2)] : (s[n÷2]+s[n÷2+1])/2
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
    version_path = joinpath(out,"barrnap-version.txt")
    open(version_path,"w") do io
        run(pipeline(`$barrnap --version`; stdout=io,stderr=io))
    end
    bv = read(version_path,String)
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
    write_differences(joinpath(out,"nhmmer-differences.tsv"),hits,refs)
    write_differences(joinpath(out,"barrnap-differences.tsv"),hits,br)
    report = Dict("julia"=>string(VERSION),"threads"=>threads,"repeats"=>repeats,
        "cpu"=>Sys.CPU_NAME,"machine"=>Sys.MACHINE,
        "engine_source_sha256"=>bytes2hex(sha256(read(joinpath(dirname(pathof(ProfileHMMs)),"Algorithms.jl")))),
        "input_sha256"=>bytes2hex(sha256(read(fasta))),"model_sha256"=>bytes2hex(sha256(read(hmm))),
        "findrna_first_call_seconds"=>cold,"findrna_warm_seconds"=>jt,
        "nhmmer_process_seconds"=>nt,"barrnap_process_seconds"=>bt,
        "findrna_hits"=>length(hits),"nhmmer_hits"=>length(refs),"barrnap_hits"=>length(br))
    open(joinpath(out,"timings.toml"),"w") do io; TOML.print(io,report); end
    println("Median seconds: FindRNA warm=",median_time(jt),"; nhmmer=",median_time(nt),"; barrnap=",median_time(bt))
    if get(ENV,"BARRNAP_MODERN","false") == "true"
        # Different models/algorithms: report all raw calls; no fabricated parity.
        @testset "HMMER reference" begin
            @test !isempty(refs)
            compare_hits(hits,refs)
        end
        if isempty(ARGS)
            @testset "Modern Barrnap synthetic fixture" begin
                @test !isempty(br)
                compare_hits(hits,br)
            end
        else
            println("Modern Barrnap real-genome calls are reported separately; CM/HMM parity is not asserted.")
        end
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
