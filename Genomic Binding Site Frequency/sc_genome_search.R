
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")


# Load genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3
genome
chr1 <- genome$chrI
chr2 <- genome$chrII
chr3 <- genome$chrIII
chr4 <- genome$chrIV
chr5 <- genome$chrV
chr6 <- genome$chrVI
chr7 <- genome$chrVII
chr8 <- genome$chrVIII
chr9 <- genome$chrIX
chr10 <- genome$chrX
chr11 <- genome$chrXI
chr12 <- genome$chrXII
chr13 <- genome$chrXIII
chr14 <- genome$chrXIV
chr15 <- genome$chrXV
chr16 <- genome$chrXVI

# Load ZF library binding sequences
seqmatrix <- c('GAAGATGGT','GACGACGGC',
                'TTAGAAGTG','GAAGACGCT',
                'GAGGACGTG','GACGCTGCT',
                'GAGTGAGGA','TGGGTGGCA',
                'TGGGGTGCC','GCCGAAGAT',
                'GATGTAGCC','TTTGTTGGC',
                'TTATGGGAG','GAAGTGGTC',
                'GGGGACGTC','GTGTAGGGG',
                'GCAGGAGGT','GTAGATGGA',
                'GGAGGGGCT','GATGAAGCT',
                'GCGTGGGCG')

# Sequences correspond to:
# 13-6
# 14-3
# 21-16
# 36-4
# 37-12
# 42-10
# 43-8
# 54-8
# 55-1
# 62-1
# 92-1
# 93-10
# 97-4
# 128-2
# 129-3
# 150-4
# 151-1
# 158-2
# 172-5
# 173-3
# Zif268

# Max mismatch (Change to 0 for zero mismatch search)
mm = 1

# Function that counts instances of binding sequences
sc_genome_search <- function(seqmatrix){
  for (seq in seqmatrix){
    
    # Forward pattern counts
    seq <- DNAString(seq)
    count=countPattern(seq,chr1,max.mismatch=mm)+countPattern(seq,chr2,max.mismatch=mm)+
      countPattern(seq,chr3,max.mismatch=mm)+countPattern(seq,chr4,max.mismatch=mm)+
      countPattern(seq,chr5,max.mismatch=mm)+countPattern(seq,chr6,max.mismatch=mm)+
      countPattern(seq,chr7,max.mismatch=mm)+countPattern(seq,chr8,max.mismatch=mm)+
      countPattern(seq,chr9,max.mismatch=mm)+countPattern(seq,chr10,max.mismatch=mm)+
      countPattern(seq,chr11,max.mismatch=mm)+countPattern(seq,chr12,max.mismatch=mm)+
      countPattern(seq,chr13,max.mismatch=mm)+countPattern(seq,chr14,max.mismatch=mm)+
      countPattern(seq,chr15,max.mismatch=mm)+countPattern(seq,chr16,max.mismatch=mm)
    
    # Reverse complement counts
    seq <- reverseComplement(seq)
    count_rc=countPattern(seq,chr1,max.mismatch=mm)+countPattern(seq,chr2,max.mismatch=mm)+
      countPattern(seq,chr3,max.mismatch=mm)+countPattern(seq,chr4,max.mismatch=mm)+
      countPattern(seq,chr5,max.mismatch=mm)+countPattern(seq,chr6,max.mismatch=mm)+
      countPattern(seq,chr7,max.mismatch=mm)+countPattern(seq,chr8,max.mismatch=mm)+
      countPattern(seq,chr9,max.mismatch=mm)+countPattern(seq,chr10,max.mismatch=mm)+
      countPattern(seq,chr11,max.mismatch=mm)+countPattern(seq,chr12,max.mismatch=mm)+
      countPattern(seq,chr13,max.mismatch=mm)+countPattern(seq,chr14,max.mismatch=mm)+
      countPattern(seq,chr15,max.mismatch=mm)+countPattern(seq,chr16,max.mismatch=mm)
    
    print(count+count_rc)
  }
}

# Run function
sc_genome_search(seqmatrix)
