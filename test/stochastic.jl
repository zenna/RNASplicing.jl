using RNASplicing
Lib = RNASplicing.Library
using FASTX
using BioSequences

# Generate 1000-length RNA with 4% chance of N, 24% for A, C, G, or U
sp = SamplerWeighted(rna"ACGU", fill(0.25, 3))
seq = randseq(RNAAlphabet{4}(), sp, 50)
splice(seq, Lib.majorsos, Lib.coremachinery)