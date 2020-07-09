module Library

using BioSequences
using ..RNASplicing: OpaqueTransElement, Interaction, Spliceosome, only
# The core splicing machinery can be expressed in this framework

# # Core Machinery

"5’ splice site - 9-nucleotide motif that generally conforms to the sequence"
donorsite = biore"(C|A)AG|GU(A|G)AGU"rna
# donorsite(seq) = occursin(biore"(C|A)AG|GU(A|G)AGU"rna, seq) ? 1.0 : 0.0
u1 = OpaqueTransElement(:u1)
u2 = OpaqueTransElement(:u2)
u2af = OpaqueTransElement(:u2af)
e1 = Interaction(donorsite, only(u1), 0, 10.0) 
coremachinery = [e1]

acceptsite = biore"A(C|U)AG"rna

# # SREs

# `seq` have a sequence of at least three Gs?
grun = biore"GGG+"rna

# Major splicesosome
majorsos = Spliceosome(Dict(u1 => 1.0))
end