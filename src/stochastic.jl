# # A Stochastic Splicing Simulator

# ## Overview
# Slicing is a process by which subsequences of preRNA (called introns, by definition) are removed.
# The remaining subsequences (called exons) are rejoined together.
# The resulting RNA sequence with introns removed is called mature RNA (mRNA).

# The splicing process is often stochastic.  Given a preRNA sequence, there is a
# distribution over corresponding mRNA sequences.
# As is almost always the case, the stochasticity results from incomplete information.
# Sources of uncertainty:
# - Structure: RNA modifications / bulged nucleotides make pairing of rNA and snRNPs variable.
# - Mutation: Single nucleotide variation (Q: At what point do mutations happen?)
# - Spliceosome composition: Variability in RNA-Binding Protein concentrations of cell

# ## Model

# At a high level, this is an abstract biophysical model that simulates splicing using a discrete-time simulation.

# __Trans-elements__ are small nucleus RNA, proteins and their combination called small- nuclear
# ribonuclear proteins (snRNP).  Trans elements perform splicing by binding to the RNA
# in a sytematic fashion.

using BioSequences, Random

"A Trans-Element / Factor: Proten or snRNP that can bind to primanry RNA transcript"
abstract type TransElement end

# There are many properties of a trans-element that may affect its role in splicing.
# The simplest representation associates a transelement only with a name:

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
  # FIXME: Ensure concentratios are normalized
end

"Is transelement `te` in the splicesosome `sos`?"
insos(te, sos::Spliceosome) = te in keys(sos.elems)
# Spliceosome(elems::Dict{T, R}) where {T, R} = Spliceosome{T, R}(elems)

"concentration of transelement `te` in spliceosome `so`"
concentration(sos, te) = sos.elems[te]

"Number of trans-elements in splicesosome"
ntranselems(sos::Spliceosome) = length(sos.elems)

# The __primary transcript__ is a sequence of RNA to be spliced by the spliceosome
const RNASequence = BioSequence{RNAAlphabet{4}}

# Trans-elements can bind to sequences in primary transcript.
# A __bound sequence__ or transcript is a primary transcript with any number of trans-elements bound to it.

"Position -- Natural number 1:length(seq)"
const Pos = Int

"Primary transcript with tran-elements bound to it"
struct BoundSeq{T <: TransElement}
  seq::RNASequence             # Primary transcriipt
  binds::Vector{Tuple{Pos, T}} # position -> transelement bound at that position
end

BoundSeq(seq::RNASequence, ::Type{T} = TransElement) where T =
  BoundSeq(seq, Tuple{Pos, T}[])

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
# We model Interactions as tuples of the form $(p, e)$, where $p$ is a __pattern__ (described below)
# and $e$ is a set of __effects__:

struct Interaction{P, TE, F}
  pattern::P
  te::TE
  loc::Pos
  factor::F
end

# #### Patterns

# A pattern (or motif) is a property of an RNA sequence.
# Formally a pattern is a set of strings, i.e., a language.
# For instance, a g-run is the set of strings $\{G, GG, GGG, GGGG, GGGGG, \dots\}$.
# A pattern $p$ __matches__ a string $s$ if any element of $p$ is a subsequence of $s$.

# We assume that the extent to which a pattern has occured is graded.
# A __pattern weight function__ maps patterns to a non-negative real number, denoting the degree to which it has matched.
# p(seq) = 1 indicates a perfect match and p(seq) = 0 indicates no match at all.

# #### EFfects

# An effect $E \subseteq \textrm{TransElements} \times \textrm{Location} \times \mathbb{R}$ is a ternary relation:
# a set of tuples of the form $(t, l, f)$, where $t$ is a transelement, $l$ is an integer valued location, and $f$ is a factor.
# Informally, the semantics of an interaction $(p, \{(t_1, l_1, f_1), (t_2, l_2, f_2), \dots\})$.
# are that if a pattern $p$ is found, then the probability that transelement $t_i$ will bind to location $l_i$ is increated (or decreased) multiplicatively by factor $f_i$.

"Returns a protein extractor that returns the singleton `te` iff its in the splicesosome"
only(te::T) where T = sos -> insos(te, sos) ? [te] : T[]

"Initial probabilities for (TransElement, Position)"
initT(sos::Spliceosome{T, R}, seq) where {T, R} = Dict{Tuple{T, Pos}, R}() 

"Absolute position where pattern matched"
matchpos(m::BioSequences.RE.RegexMatch) = m.captured[1]

"Simulate one binding step"
function step(rng, bseq, sos, intrs)
  T = initT(sos, bseq.seq)

  # For each interaction
  for intr in intrs
    # and for each match of the interaction's pattern
    for x in eachmatch(intr.pattern, bseq.seq)
      # For each  
      abspos = matchpos(x) + intr.loc
      for protein in intr.te(sos)
        if (protein, abspos) ∉ keys(T)
          T[protein, abspos] = 1.0
        end
        T[protein, abspos] *= intr.factor
      end
    end
  end
  bseq = bind(rng, T, bseq, sos)
  maybesplice(bseq)
end

"Simulate a binding process"
function bind(rng, T, bseq, sos)
  @show T
  # Select transelement with probability proportional to concentration
  # bseq

  # How 
end

"Simulate core splicing machinery -- checks for presence of u1 and u2"
maybesplice(bseq) = bseq

# The splicosome is then simulated by simply repeating this process until some criterion of convergence is me
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