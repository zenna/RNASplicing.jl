Files:
 - GCF_000146045.2_R64_genes.fsa
      - fasta file with the sequences of all the annotated genes
      - structure: two lines alternating as follows
        >gene_name GeneID:geneID Chromosome:chromosomeID Start:gene_start Stop:gene_stop Strand:strand
        sequence of gene
      - meanings: 
             - geneID: gene ID number
             - chromosomeID: chromosome ID
             - gene_start: the first base of the gene relative to the chromosome sequence
             - gene_stop: the last base of the gene relative to the chromosome sequence
             - strand: whether this sequence came from the + or - strand. sequences that came from the - have been corrected to match the + strand
 - GCF_000146045.2_R64_genes_with_introns.fsa
      - the same as above, but only has intron containing genes
 - GCF_000146045.2_R64_introns.txt
      - list of all annotated introns
      - structure: four tab-separated columns, as follows: geneID \t gene_name \t intron_start \t intron_end
      - meanings: 
            - geneID: gene ID number
            - gene_name: the systematic name of the gene
            - intron_start: first base of intron in gene. 1-based coordinates for sequences in GCF_000146045.2_R64_genes
            - intron_stop: last base of intron in gene. 1-based coordinates for sequences in GCF_000146045.2_R64_genes
 - GCF_000146045.2_R64_chromosomes.gff
      - conversion between IDs for chromosomes and their numerals 
      - structure: two tab-separated columns, as follows: chromosomeID \t chromosome=chr_numeral
      - meanings:
            - chromosomeID: chromosome ID
            - chr_numeral: the roman numeral for that yeast chromosome
 - splice_sites.txt
      - contains the 5'SS and 3'SS ends for each annotated intron
      - structure: gene_name \t first 6 bases of intron \t last 80 bases of intron or, if shorter, whole intronic sequence (padded with spaces)
      - meanings:
             - gene_name: the systematic name of the gene
      - motifs matching the canonical branch point sequence are capitalized

Example: YBL018C

Fasta entry:
>YBL018C GeneID:852263 Chromosome:NC_001134.8 Start:185998 Stop:186474 Strand:-
ATGGGGAAAAAGACTTTTAGAGAATGGCAATATTTCAAGTTATCAATGTATGTATATTTTTGACTTTTTGAGTCTCAACTACCGAAGAGAAATAAACTACTAACGTACTTTAATATTTATAGTACTTCATTCGATCAAGATGTGGACGATGCACATGCTATTGATCAAATGACATGGAGGCAATGGTTAAATAACGCATTGAAACGTAGTTATGGTATTTTTGGTGAAGGTGTAGAATATTCATTCCTGCATGTTGACGATAAGCTGGCTTACATCAGAGTAAATCATGCGGATAAAGATACATTTTCTTCGTCCATCAGTACATACATATCTACTGATGAACTCGTCGGTTCACCATTAACAGTGTCAATTCTACAAGAATCTTCCAGTTTGAGACTTCTGGAGGTTACTGACGATGACCGCCTATGGCTGAAAAAAGTAATGGAAGAAGAAGAACAAGACTGTAAATGTATATAG

Introns entry:
852263  YBL018C 48      122

Chromosomes entry: 
NC_001134.8	chromosome=II

Splice sites entry:
YBL018C gtatgt  tcaatgtatgtatatttttgactttttgagtctcaactaccgaagagaaataaacTACTAACgtactttaatatttatag


So, this is a gene on chromosome II with one intron between bases 48 and 122 that contains a canonical branch point sequence motif. The intron is represented below with lower case letters within the full sequence 

ATGGGGAAAAAGACTTTTAGAGAATGGCAATATTTCAAGTTATCAATgtatgtatatttttgactttttgagtctcaactaccgaagagaaataaactactaacgtactttaatatttatagTACTTCATTCGATCAAGATGTGGACGATGCACATGCTATTGATCAAATGACATGGAGGCAATGGTTAAATAACGCATTGAAACGTAGTTATGGTATTTTTGGTGAAGGTGTAGAATATTCATTCCTGCATGTTGACGATAAGCTGGCTTACATCAGAGTAAATCATGCGGATAAAGATACATTTTCTTCGTCCATCAGTACATACATATCTACTGATGAACTCGTCGGTTCACCATTAACAGTGTCAATTCTACAAGAATCTTCCAGTTTGAGACTTCTGGAGGTTACTGACGATGACCGCCTATGGCTGAAAAAAGTAATGGAAGAAGAAGAACAAGACTGTAAATGTATATAG


