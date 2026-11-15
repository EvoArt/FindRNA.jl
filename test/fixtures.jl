# Fully deterministic synthetic fixtures derived from the bundled profiles.
# Both reference tools and FindRNA receive exactly these same sequences.
using BioSequences, ProfileHMMs
function make_fixture(path)
    ps = parse_hmm_file(FindRNA.default_profiles())
    seqs = [LongDNA{4}(join("ACGT"[argmin(p.match_emissions[:,k])] for k in 1:p.length)) for p in ps]
    open(path,"w") do io
        for (p,s) in zip(ps,seqs)
            println(io,">",p.name,"_forward\n", "N"^80,s,"N"^80)
            println(io,">",p.name,"_reverse\n", "N"^80,reverse_complement(s),"N"^80)
        end
        println(io,">two_16S\n", "N"^80,seqs[1],"N"^150,seqs[1],"N"^80)
        println(io,">negative\n", "ACGT"^500)
        println(io,">ambiguous\n", "N"^2000)
    end
    ps,seqs
end
