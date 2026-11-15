module FindRNA

using BioSequences, FASTX, Tables
using ProfileHMMs
import ProfileHMMs: _title, _kv, _table, _more, _nshow, _g

export findrna, RNAResult, extract_sequences, write_fasta, write_gff, default_profiles
export setprogress!

"RNA hits with vector indexing and Tables.jl row access. Coordinates are 1-based inclusive."
struct RNAResult <: AbstractVector{SearchHit}
    hits::Vector{SearchHit}
    source::String
end
Base.size(r::RNAResult) = size(r.hits)
Base.getindex(r::RNAResult, i::Int) = r.hits[i]
Base.IndexStyle(::Type{RNAResult}) = IndexLinear()
Tables.istable(::Type{RNAResult}) = true
Tables.rowaccess(::Type{RNAResult}) = true
Tables.rows(r::RNAResult) = r.hits
Tables.schema(::RNAResult) = Tables.Schema(fieldnames(SearchHit), fieldtypes(SearchHit))

"Path to the bundled bacterial 16S, 23S and 5S rRNA profiles (Rfam via barrnap 0.9)."
default_profiles() = joinpath(@__DIR__, "..", "data", "bac.hmm")

"""
    findrna(fasta; profiles=default_profiles(), models=nothing, progress, kwargs...)

Find nucleotide profile matches on both strands. `models=["16S_rRNA"]` selects
16S only. Returns a vector-like, Tables.jl-compatible `RNAResult`.
The current engine uses Viterbi scoring, not HMMER's Forward/domain pipeline;
its approximate E-values must not be interpreted as calibrated nhmmer E-values.
See `ProfileHMMs.nhmmscan` for search options.
"""
function findrna(fasta::AbstractString; profiles=default_profiles(), kwargs...)
    RNAResult(nhmmscan(profiles, fasta; kwargs...), abspath(fasta))
end

"Extract strand-oriented LongDNA sequences, in result order, for AssignTaxonomy."
function extract_sequences(r::RNAResult; fasta=r.source)
    wanted = Set(h.target_name for h in r)
    seqs = Dict{String,LongDNA{4}}()
    open(FASTA.Reader, fasta) do reader
        for record in reader
            id = String(identifier(record))
            id in wanted || continue
            haskey(seqs, id) && throw(ArgumentError("duplicate FASTA identifier: $id"))
            seqs[id] = LongDNA{4}(sequence(record))
        end
    end
    map(r.hits) do h
        haskey(seqs, h.target_name) || throw(ArgumentError("missing sequence: $(h.target_name)"))
        seq = seqs[h.target_name][h.target_start:h.target_end]
        h.strand == '-' ? reverse_complement(seq) : seq
    end
end

function write_fasta(io::IO, r::RNAResult; kwargs...)
    for (i, (h, seq)) in enumerate(zip(r, extract_sequences(r; kwargs...)))
        println(io, ">RNA_", i, " ", h.target_name, ":", h.target_start, "-", h.target_end,
                "(", h.strand, ") ", h.query_name)
        println(io, seq)
    end
    nothing
end

# Escape reserved characters byte-wise, including non-ASCII identifiers.
function _gff_escape(s)
    join((b in UInt8.(vcat(collect('a':'z'), collect('A':'Z'), collect('0':'9'),
                          ['.', ':', '^', '*', '$', '@', '!', '+', '_', '?', '-', '|'])) ?
          string(Char(b)) : "%" * uppercase(string(b; base=16, pad=2))) for b in codeunits(s))
end

"Write GFF3; the score column contains the approximate Viterbi E-value."
function write_gff(io::IO, r::RNAResult)
    println(io, "##gff-version 3")
    for (i, h) in enumerate(r)
        feature = endswith(h.query_name, "rRNA") ? "rRNA" : "RNA"
        println(io, join((_gff_escape(h.target_name), "FindRNA", feature,
                h.target_start, h.target_end, h.evalue, h.strand, '.',
                "ID=RNA_$i;Name=$(_gff_escape(h.query_name))"), '\t'))
    end
    nothing
end

for f in (:write_fasta, :write_gff)
    @eval function $f(path::AbstractString, r::RNAResult; kwargs...)
        open(path, "w") do io
            $f(io, r; kwargs...)
        end
        path
    end
end

Base.show(io::IO, r::RNAResult) = print(io, "RNAResult(", length(r), " hits)")
function Base.show(io::IO, ::MIME"text/plain", r::RNAResult)
    _title(io, "RNAResult", "$(length(r)) RNA hits")
    _kv(io, "source", basename(r.source))
    isempty(r) && return
    n = _nshow(io, length(r), 8)
    rows = [[first(h.target_name, 24), h.query_name, "$(h.target_start):$(h.target_end)",
             string(h.strand), _g(h.bit_score)] for h in r.hits[1:n]]
    _table(io, ["sequence", "model", "coordinates", "strand", "bits"], rows;
           align=[:l, :l, :r, :r, :r])
    _more(io, length(r), n, "hits")
end
end
