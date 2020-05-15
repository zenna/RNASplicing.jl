# # A Stochastic Splicing Simulator

# ## Overview
# Slicing is a process by which subsequences of preRNA (called introns, by definition) are removed.
# The remaining subsequences (called exons) are rejoined together.
# The resulting RNA sequence with introns removed is called mature RNA (mRNA)

# The splicing process is stochastic.  Given a preRNA sequence, there is a
# distribution over corresponding mRNA sequences.
# Sources of stochasticity:
# - Structure: RNA modifications /bulged nucleotides make pairing of rNA and snRNPs variable
# - Mutation: Single nucleotide variation (Q: At what point do mutations happen?)
# - Composition: Variability in RNA-Binding Protein concentrations of cell

# ## Model

# Trans-elements are small nucleus RNA, proteins and their combination called small-
# nuclear ribonuclear proteins (snRNP).  Trans elements perform splicing by binding to the RNA
# in a sytematic fashion.

using BioSequences

"A Trans-Element / Factor: Proten or snRNP that can bind to primanry RNA transcript"
abstract type TransElement end

# There are many properties of a trans element that may affect its role in splicing.
# The simplest representation associated a transelement only with a name:

"Symbolic TransElement -- represents only its name"
struct OpaqueTransElement <: TransElement
  name::Symbol
end

# A __spliceosome__ is a complex of transelements.
# We model a spliceosome as a partial function from trans-elements to real values representing their
# relative concentrations:

"Collection of transelements at different concentrations"
struct Spliceosome{T <: TransElement, R <: Real}
  elems::Dict{T, R}   # Maps from transelement to its concentration in spliceosome
end

"concentration of transelement `te` in spliceosome `so`"
concentration(so, te) = so.elems[te]

# The __primary transcript__ is a sequence of RNA to be spliced by the spliceosome
const RNASequence = BioSequence{RNAAlphabet{4}}

# Trans-elements can bind to sequences in primary transcript.
# The following data structure represents a transcript bound with any number of trans-elements:

"Primary transcript with tran-elements bound to it"
struct BoundSeq{T <: TransElement}
  seq::RNASequence                        # Primary transcriipt
  binds::Vector{Tuple{Int, T}} # position -> transelement bound at that position
end

# Each bound sequence `bseq` has an associated energy `ℓ(bseq)`
# We assume that ℓ has a particular structure: proteins bound to particular sites on
# the primary transcript modulate the energy function locally

# A pattern (or motif) is a property of an RNA sequence.
# We assume that the extent to which a pattern has occured is graded.
# Formally a pattern is a function from a subsequence to the unit interval, where
# p(seq) = 1 indicates a perfect match and p(seq) = 0 indicates no match at all.

# An example of a hard pattern is a G-run, a sequence of Gs:

"Does `seq` have a sequence of at least three Gs?"
grun(seq) = occursin(biore"GGG+"rna, seq) ? 1.0 : 0.0

# Splicing literature makes a distinction between the core splicing machinary
# and __splicing regulary elements__.
# The core splicing machinary is carried out by the __major splicesosome__, which
# consists of exactly 5 rnSNAs and around 150 proteins, which combine to make
# 5 rnSNPs: U1, U2, U4, U5, U6.
# The core splicing mechanism is an ordered, stepwise assembly of discrete snRNP particles
# on the pre-mRNA primary transcript.  First U1 binds to the 5' splice site, U2 binds
# to a branch site, and so on until a split occurs.

# Splicing regulartory elements (SREs) are patterns on RNA that inhibit or promote splicing machinery.
# They are believed to be the main mechanism for alternative splcing.

# We will model both the core machinery and SREs within the same framework.
# We model SQQs as tuples of the form $(p, e)$, where $p$ is a pattern (described above)
# and $e: \to \mathbb{R}$ is an effect -- a function that maps
# The existance of a pattern can increase or decrease the binding energy

struct SQQ{P, E}
  pattern::P
  effect::E
end

# The core splicing machinary can be expressed in this framework

"5’ splice site - 9-nucleotide motif that generally conforms to the sequence"
donorsite(seq) = occursin(biore"(C|A)AG|GU(A|G)AGU"rna, seq) ? 1.0 : 0.0
const e1 = SQQ(donorsite, Dict(OpaqueTransElement(:u1) => 0.1)) 
const coremachinary = [e1]

"The energy of a configuration `bseq` assuming splicing regulary elements `sres`"
function ℓ(bseq::BoundSeq, sqqs)
  while !converged() # FIXME
    for sqq in sqqs
      sqq = 3
    end
  end
end

# The splicosome can be simulated
"Simulate the splicing"
function simsplice(rng, bseq, so)
  while !converged()
    bseq = stableconfig(rng, ℓ, bseq, so)
    bseq = maybesplice(bseq)
  end
  bseq
  #FIXME: Might want introns too
end

"Simulate core splicing machinary"
function maybesplice(bseq)
  bseq   # Check for presence of u1 and u2
end

"Stochastic optimization to find stable (energetically minimal wrt ℓ) configuration of binding"
function stableconfig(rng, ℓ, bseq, so)
end

# Model
# A splicing round occurs when U1 and U2AF proteins reach a global minimum

# # Measurement (Sequencing)

# Unfortunately, we do not have exact data, and hence we must model the stochasticity in the sequencing process.
# The following model of sequencing assumes that from a given `rna`, `n` (possibly differing, due to alternative splicing) mRNAs
# are produced, and the sequencing measures a single sunsequence (of random length) from each

"Sequence rng process"
function sequence(rng, seq, n)
  sequences = RNASequence[]

  # Distribution over start point
  splice_start(rng) = uniform(rng, 1:length(seq))

  # Distribution over length
  splce_length(rng) = uniform(rng, 1:100)

  for i = 1:n
    mRNA = simsplice(rng, seq, n)
    lb = splice_start(rng)
    ub = lb + splice_length(rng)
    push!(sequences, mRNA[lb:ub]) # FIXME: Might go out of bounds
  end
  sequences
end