# A stochastic model of splice simulation
import BioSequences: biores, rna, BioSequence, RNAAlphabet

# Slicing is a process by which subsequences of preRNA (called introns, by definition) are removed split into matureRNA 

# The splicing process is stochastic.  Given a preRNA sequence, there is a
# distribution over corresponding mRNA sequences.
# Sources of stochasticity:
#   Structure: RNA modifications /bulged nucleotides make pairing of rNA and snRNPs variable
#   Mutation: Single nucleotide variation (Q: At what point do mutations happen?)
#   Composition: Variability in RNA-Binding Protein concentrations of cell

# Trans-elements -- small nucleus RNA, proteins and their combination small-
# nuclear ribonuclear proteins (snRNP) perform splicing by binding to the RNA
# in a sytematic fashion.

abstract type TransElement end

"Symbolic TransElement -- represents only its name"
struct OpaqueTransElement <: TransElement
  name::Symbol
end

# A spliceosome is a complex of transelements.

"Collection of transelements at different concentrations"
struct Spliceosome{T <: TransElement, R <: Real}
  elems::Dict{T, R}   # Maps from transelement to its concentration in spliceosome
end

"concentration of transelement `te` in spliceosome `so`"
concentration(so, te) = so.elems[te]

# The primary transcript is a sequence of RNA to be spliced by the spliceosome
const RNASequence = BioSequence{RNAAlphabet{4}}

# Trans-elements can bind to sequences in primary transcript

"Primary transcript with tran-elements bound to it"
struct BoundSeq{T<:TransElement}
  seq::RNASequence                        # Primary transcriipt
  binds::Vector{Tuple{Int, T}} # position -> transelement bound at that position
end

# Each (potential) bound sequence has an associated energy ℓ
# We assume that ℓ has a particular structure:
# A pattern (or motif) is a property of a RNA sequence
# We assume that the extent to which a pattern has occured is graded
# Formally a pattern is a function from a subsequence to the unit interval, where
# p(seq) = 1 indicates a perfect match and p(seq) = 0 indicates no match at all

# An example of a hard pattern is a G-run, a sequence of Gs

"Does `seq` have a sequence of at least three Gs?"
grun(seq) = occursin(biore"GGG+"rna, seq) ? 1.0 : 0.0

# Splicing literature makes a distinction between the core splicing machinary
# and splicing regulary elements.
# The core splicing machinary is carried out by the major splicesosome, which
# consists of exactly 5 rnSNAs and around 150 proteins, which combine to make
# 5 rnSNPs: U1, U2, U4, U5, U6.
# In the core splicing mechanism is an ordered, stepwise assembly of discrete snRNP particles
# on the pre-mRNA primary transcript.  First U1 binds to the 5' splice site, U2 binds
# to a branch site, and so on until a split occurs.

# Splicing regulartory elements (SREs) are patterns on RNA
# In particular, they inhibit or promote splicing machinery

# We will model both of these within the same framework
# We model SQQs as tuples of the form $(p, e)$, where $p$ is a pattern (described above)
# and $e: \to mathbb{R}$ is an effect -- a function that maps
# The existance of a pattern can increase or decrease the binding energy

struct SRE{P, E}
  pattern::P
  effect::E
end

# The core splicing machinary can be expressed in this framework

"5’ splice site - 9-nucleotide motif that generally conforms to the sequence"
donorsite(seq) = occursin(biores"(C/A)AG/GU(A/G)AGU", seq) ? 1.0 : 0.0
const e1 = SQQ(donor, Dict(TransElement(:u1) => 0.1)) 
const coremachinary = [e1]

"The energy of a configuration `bseq` assuming splicing regulary elements `sres`"
function ℓ(bseq::BoundSeq, sqqs)
  while !converged() # FIXME
    for sqq in sqqs
      # Now what? Bind the sqq 
    end
  end
end

# The splicosome can be simulated
"Simulate the spliceosome"
function simspliceosome(rng, bseq, so)
  while !converged()
    bseq = stableconfig(rng, ℓ, bseq, so)
    bseq = maybesplice(bseq)
  end
  bseq
  #FIXME: Might want introns too
end

"Simulate core splicing machinary"
function maybesplice(bseq)
  # Check for presence of u1 and u2
  bseq
end

"Stochastic optimization to find stable (energetically minimal wrt ℓ) configuration of binding"
function stableconfig(rng, ℓ, bseq, so)
end

# Model
# A splicing round occurs when U1 and U2AF proteins reach a global minimum