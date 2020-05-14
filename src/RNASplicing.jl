module RNASplicing

using BioSequences

# The biology 
# A gene is a DNA
# DNA is am molecule composed of two helixes coiled around one another.
# Each DNA stand is a polynucleotide -- composed of nucleiotides.
# Each nucleotide is is composod of one of four nitrogen-containing nucleuobases (cytosine [C], guanine [G], adenine [A] or thymine [T]), a sugar called deoxyribose, and a phosphate group
# The bases between stands are bound, (A on 1 stand implies that T is on the other in the corresponding position. C with G pair in the same way.
# Stands of DNA have a notion of directionality
# The 5′-end (pronounced "five prime end") designates the end of the DNA or RNA strand that has the fifth carbon
# The 3′-end (three prime end) of a strand is so named due to it terminating at the hydroxyl group of the third carbon in the sugar-ring, and is known as the tail end.  

# RNA stands are created from DNA stands in transcription
# RNA is single straded
# Splicing is the removal of introns

# Overall model
# Input is a sequence of RNA
# RNA interacts with spliceosome -- a complex of small nuclear ribonucleo proteins (snRNPs).
# Each snRNP is represented as a Motif
# A Motif 
# Splicing is regulated mainly by RNA-binding proteins (RBPs)
# RBPs often act according to positional principles defined by an RNA splicing map to enhance or repress exon inclusion  

# sNRPS can
# Bind to a splice site
# "Help" another protein bind, e.g U2AF, when bound, helps (allows?) BBP to bind to branch site
  # Q does it promote particular branch site?
  # Q Constraints on branch site? A base? 
# Replace another protein.  U2 Replaced BBP.  Then whats the point of BBP?
# Change some state of the molecule (Adenine bulges out when bound to by U2)
# "Recruit" other snRNPS
# U6 pushed out U1 from branch site
# When U6 tries to interact with U2 (U4) which are in the middle, is release
# interaction of u6 and u2 can cause rna to split, and loop to form
"""
A spliceosome is composed of five small nuclear RNAs (snRNA) and a range (around 150) of associated protein factors.
snRNA is usually about 100-300 nucleotides in length
RNAs combine with the protein factors to make-protein complexes called snRNPs (small nuclear ribonucleo proteins, pronounced "snurps").
snRNA gives specificity to individual introns by "recognizing" the sequences of critical splicing signals at the 5' and 3' ends and branch site of introns
"""
struct Splicesosome
end

"""
A sequence is a discrete sequence or combination of base juxtapositions
found in naturally occurring RNAs in unexpectedly high abundance
"""
struct Motif
  pattern::Pattern
end

"Splice sequence `seq`"
function splice(seq)
  motifloc = matchmotifs(seq)
  updateprobs!(motifloc)
  boundaries(motifloc)
end

"""
"""
function matchmotifs(seq, motifs::AbstractVector{<:Motif})
  motifloc = Dict{Int, Motif} # Mapping from site to SRE bound to protein bound to that site
  n = length(seq)
  for i = 1:n
    for m in motifs
      # Q: What if multiple matches?
      if match(seq, i, m.pattern)
        motifloc[i] = (id = m.id, pos = i, prob = prob(seq, i, m.pattern))
      end
    end
  end
  motifloc
end

"Probability that seq[i:end] matches pattern"
function prob(seq, i, pattern)
end

function initconverged(motif)
  1
end

function converged(motifloc, state; n = 10)
  state > n, state + 1
end


function updateprobs!(motifloc)
  # The presence of one splicesosome may inhibit or surpress others
  state = initconverged(motifloc)
  while !converged(motifloc, state)
    for i in motifloc
      for j in motifloc 
        pi = motifloc[i]
        pj = motifloc[j]
        motifloc[j].prob *= effect(i, pi.id, pi.prob, j, pj.id, pj.prob)
      end
    end
  end
  motifloc
end

"""
Multiplicative effect that motif `id1` has on the probability that of motif `id2`
given that motif1 is in position `p1` and motif `p2` is in position `p2`.
"""
function effect(p1, id1, prob1, p2, id2, prob2, T)
  t = T[id1, id2]
  f(t, abs(p2 - p1), prob1)
end

"Influnce "
function f(t, δ)

end

"Compute the boundaries."
function boundaries(motifloc; thresh = 0.7)
  startx = []
  endx = []
  for i in motifloc
    pi = motifloc[i]
    if (pi.id == U1 && pi.prob > thresh)
      push!(startx, pi)
    end
    if (pi.id == U2 && pi.prob > thresh)
      push!(endx, pj)
    end
  end
  startx, endx
end

end # module
