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

using BioSequences, Random

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

"Is transelement `te` in the splicesosome?"
insos(te, sos::Spliceosome) = te in keys(sos.elems)
# Spliceosome(elems::Dict{T, R}) where {T, R} = Spliceosome{T, R}(elems)

"concentration of transelement `te` in spliceosome `so`"
concentration(sos, te) = sos.elems[te]

"Number of trans-elements in splicesosome"
ntranselems(sos::Spliceosome) = length(sos.elems)

# The __primary transcript__ is a sequence of RNA to be spliced by the spliceosome
const RNASequence = BioSequence{RNAAlphabet{4}}

# Trans-elements can bind to sequences in primary transcript.
# The following data structure represents a transcript bound with any number of trans-elements:

# Position -- Natural number 1:length(seq)
const Pos = Int

"Primary transcript with tran-elements bound to it"
struct BoundSeq{T <: TransElement}
  seq::RNASequence             # Primary transcriipt
  binds::Vector{Tuple{Pos, T}} # position -> transelement bound at that position
end

BoundSeq(seq::RNASequence, ::Type{T} = TransElement) where T =
  BoundSeq(seq, Tuple{Pos, T}[])

# Each bound sequence `bseq` has an associated energy `ℓ(bseq)`
# We assume that ℓ has a particular structure: proteins bound to particular sites on
# the primary transcript modulate the energy function locally

# ### Interactions 
# Splicing literature makes a distinction between the core splicing machinery
# and __splicing regulary elements__.
# The core splicing machinery is carried out by the __major splicesosome__, which
# consists of exactly 5 rnSNAs and around 150 proteins, which combine to make
# 5 rnSNPs: U1, U2, U4, U5, U6.
# The core splicing mechanism is an ordered, stepwise assembly of discrete snRNP particles
# on the pre-mRNA primary transcript.  First U1 binds to the 5' splice site, U2 binds
# to a branch site, and so on until a split occurs.

# Splicing regulartory elements (SREs) are patterns on RNA that inhibit or promote splicing machinery.
# They are believed to be the main mechanism for alternative splcing.

# We will model both the core machinery and SREs within the same framework.
# We model Interactions as tuples of the form $(p, e)$, where $p$ is a __pattern__ (described above)
# and $e is a set of __effects__:

struct Interaction{P, E}
  pattern::P
  effects::E
end

# #### Patterns

# A pattern (or motif) is a property of an RNA sequence.
# Formally a pattern is a set of strings, i.e., a language.
# For instance, a g-run is the set of strings $\{G, GG, GGG, GGGG, GGGGG, \dots\}$.
# A pattern $p$ __matches__ a string $s$ if any element of $p$ is a subsequence of $s$.

# We assume that the extent to which a pattern has occured is graded.
# A __pattern weight function__ maps patterns to a non-negative Real number, denoting the degree to which it has matched.
# p(seq) = 1 indicates a perfect match and p(seq) = 0 indicates no match at all.

# #### EFfects

# An effect $E \subseteq TransElements \times Position \times \mathbb{R}$ is a ternary relation:
# a set of tuples of the form $(t, l, f)$, where $t$ is a transelement, $l$ is an integer valued location, and $f$ is a factor.
# Informally,. the semantics of an interaction $(p, \{(t_1, l_1, f_1), (t_2, l_2, f_2), \dots\}$.
# are that if a pattern $p$ is found, then the probability that transelement $t_i$ will bind to position $p_i$ is increated (or decreased) multiplically by factor $f_i$

"Returns a protein extractor that returns the singleton `te` iff its in the splicesosome"
only(te::T) where T = sos -> insos(te, sos) ? [te] : T[]

"Initial probabilities for (TransElement, Position)"
initT(sos::Spliceosome{T, R}, seq) where {T, R} = Dict{Tuple{T, Pos}, R}() 

"Absolute position where pattern matched"
matchpos(m::BioSequences.RE.RegexMatch) = m.captured[1]

"Simulate one binding step"
function step(rng, bseq, sos, intrs)
  T = initT(sos, bseq.seq)

  for intr in intrs
    for x in eachmatch(intr.pattern, bseq.seq)
      for (p, relpos, mul) in intr.effects
        abspos = matchpos(x) + relpos
        for protein in p(sos)
          if (protein, abspos) ∉ keys(T)
            T[protein, abspos] = 1.0
          end
          T[protein, abspos] *= mul
        end
      end
    end
  end
  bseq = bind(rng, T, bseq)
  maybesplice(bseq)
end

"Simulate a binding process"
bind(rng, T, bseq) = ..

# The splicosome can be simulated
"Simulate the splicing"
function splice(rng, bseq::BoundSeq, sos, intrs; n = 1000)
  for i = 1:n
    bseq = step(rng, bseq, sos, intrs)
    bseq = maybesplice(bseq)
  end
  bseq
end

splice(bseq, sos, intrs) = splice(Random.GLOBAL_RNG, bseq, sos, intrs)
splice(rng, rnaseq::RNASequence, sos, intrs) = splice(rng, BoundSeq(rnaseq), sos, intrs)

"Simulate core splicing machinery -- checks for presence of u1 and u2"
maybesplice(bseq) = bseq


# Model
# A splicing round occurs when U1 and U2AF proteins reach a global minimum

# # Measurement (Sequencing)

# Unfortunately, we do not have exact data, and hence we must model the stochasticity in the sequencing process.
# The following model of sequencing assumes that from a given `rna`, `n` (possibly differing, due to alternative splicing) mRNAs
# are produced, and the sequencing measures a single sunsequence (of random length) from each

"Sequence RNA process"
function sequence(rng, seq, n)
  sequences = RNASequence[]

  splice_start(rng) = uniform(rng, 1:length(seq))

  splce_length(rng) = uniform(rng, 1:100)

  for i = 1:n
    mRNA = simsplice(rng, seq, n)
    lb = splice_start(rng)
    ub = lb + splice_length(rng)
    push!(sequences, mRNA[lb:ub]) # FIXME: Might go out of bounds
  end
  sequences
end