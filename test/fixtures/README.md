# Real-sequence regression fixture

`pseudomonas-5s.fna` is NCBI accession **CP148752.1**, coordinates **116377–116486**
(1-based inclusive, plus strand). Retrieved from NCBI EFetch on 2026-11-15:

https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP148752.1&seq_start=116377&seq_stop=116486&rettype=fasta&retmode=text

The 110-base sequence exercises a real 5S locus rejected by the ungapped filter
when the geometric length term was missing from the null score. The test embeds
it in a 5,000-base window with neutral flanks.
