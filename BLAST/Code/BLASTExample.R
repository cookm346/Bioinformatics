#set working directory
    setwd("D:/BLASTp")
    
    #optional:
        #read data
        d <- read.csv("ProteinsPeptides.csv", header = TRUE, stringsAsFactors = FALSE)
  
        #extract summary of unique peptides frequencies for each protein
        s <- peptides(d)

        #convert .csv protein and peptide file to .fasta file
        fasta(s)
    
#run BLAST+ search
    result <- blast(1)            #top hit alignment for each peptide
    result <- blast(10)           #top ten hit alignments
    result <- blast(e = .001)     #top ten hit alignments and e < .001
    result <- blast(1, e = .001)  #top one hit alignment and e < .001
    