#GENERATE FAKE DATA WITH 

    data <- data.frame(Protein = rep(paste("protein", letters[1:3]), c(2,3,5)), Peptide = replicate(10, paste(sample(toupper(letters[1:2]), 3, replace = TRUE), collapse = "")))
    data     #view data
    
#IF YOU ALREADY HAVE A FASTA FILE, YOU CAN SKIP THESE STEPS
    
    #GET SUMMARY OF EACH PROTEINS UNIQUE PEPTIDES WITH PEPTIDE FREQUENCY
        peptides(data)
        peps <- peptides(data) #save as peps
    
    #CONVERT DATA TO FASTA FILE
        getwd()                                        #this is where the .fasta file will be saved
        setwd("C:/Users/HIV-Proteomics/Desktop")       #you can change this location
        
        fasta(peps)                                    #save fasta file
        #check location of getwd() and there will be a sequences.fasta file
    
    
        
#RUN BLAST SEARCH
    result <- blast(1)            #top hit alignment for each peptide
    result <- blast(10)           #top ten hit alignments
    result <- blast(e = .001)     #top ten hit alignments and e < .001
    result <- blast(1, e = .001)  #top one hit alignment and e < .001