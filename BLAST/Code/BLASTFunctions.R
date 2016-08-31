#get unique peptides with counts for each protein
    peptides <- function(data){
        #add a column of 1s to be used as tally for each peptide
        d <- data.frame(data, Frequency = rep(1, nrow(data)))
        
        #compute peptide frequency for each protein by summing 1s for 
        # unique peptides for each protein
        s <- aggregate(data = d, Frequency ~ Protein + Peptide, sum)
        
        #return the results ordered by protein, then Peptide
        return(s[order(s$Protein, s$Peptide), ])
    }
    
    
    
#convert csv file to .fasta file
    fasta <- function(data){
        
        #clean protein names to remove commas and spaces
        # commas and spaces can mess up the BLAST output
        for(i in 1:nrow(data)){
            data_row <- gsub(pattern = ",", replacement = "/", x = data[i, 1])
            data[i, 1] <- gsub(pattern = " ", replacement = "-", x = data_row) 
        }
        
        #for each row, paste: > sign, protein, new line, peptide
        f <- paste(">", data[,1], "\n", data[,2], sep = "")
        
        #write all sequences to .fasta file
        write(f, "sequences.fasta")
    }
    

    
#run blast search
    blast <- function(nHits = 10, e = 10){
      
        #start time
        start <- Sys.time()
        
        #get .fasta files from working directory
        print("Reading .fasta files...")
        infiles <- list.files(pattern = "*.fasta")
        
        #remove ".fasta" from file names to be used in output names
        print("Extracting file names...")
        infilenames <- unlist(sapply(infiles, strsplit, split = "*.fasta"))
        
        #write BLAST+ command for each in file
        #this is what would normally be typed in the command line
        print("Writing BLASTp commands...")
        
        cmd <- paste('blastp -query ', infiles, ' -db db/nr.00 -out ', infilenames, '-Results.txt -num_alignments ', nHits, ' -evalue ', e,
                     ' -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sacc stitle"', sep = "")
        
        #length of cmd
        lc <- length(cmd)
        
        print("Conducting BLASTp search...")
        
        #one loop iteration for each command
        for(i in 1:lc){
            shell(cmd[i], intern = TRUE)
            
            #print loop percentage counter
            print(paste("Completed .fasta file", i, "of", lc, sep = " "))
            
        }
        
        #end time
        end <- Sys.time()
        
        #print final message with time to complete search
        dt <- round(as.numeric(difftime(end, start, units = "mins")), 2)
        print(paste("BLASTp search complete (",  dt, " minutes).", sep = ""))
        
        return(results())
    }
    
   
     
    cleanResultFiles <- function(){
        
        #read in all results files
        outfiles <- list.files(pattern = "*-Results.txt")
        
        #file number sequence
        f_seq <- 1:length(outfiles)
        
        #if any results files have size of 0 (meaning they are empty)
        #store size 0 file locations
        if(any(file.info(outfiles)$size == 0)){
            #location of empty files
            z <- which(file.info(outfiles)$size == 0)
            
            #remove empty file locations from sequence of files
            f_seq <- f_seq[!(f_seq %in% z)]
        }
        
        for(i in f_seq){
            
            f <- read.table(outfiles[i], header = FALSE, sep = ",", col.names = seq_len(30), fill = TRUE)
            
            #convert NAs to blanks
            f[is.na(f)] <- ""
            
            #paste columns 14+
            f[, 14] <- do.call(paste, f[14:ncol(f)])
            
            #remove unused columns
            f <- f[ , -(15:ncol(f))]
            
            #overwrite file with cleaned up file
            write.table(f, file = outfiles[i], row.names = FALSE, col.names = FALSE, sep = ",")
            
        }
        
    }
    
    mergeResultFiles <- function(){
        
        #read in all results files
        outfiles <- list.files(pattern = "*-Results.txt")
        
        #file number sequence
        f_seq <- 1:length(outfiles)
        
        #if any results files have size of 0 (meaning they are empty)
        #store size 0 file locations
        if(any(file.info(outfiles)$size == 0)){
            #location of empty files
            z <- which(file.info(outfiles)$size == 0)
            
            #remove empty file locations from sequence of files
            f_seq <- f_seq[!(f_seq %in% z)]
        }
        
        for(i in f_seq){
            
            f <- read.table(outfiles[i], header = FALSE, sep = ",", fill = TRUE, stringsAsFactors = FALSE)
            
            if(i == 1){
                all_files <- rbind("", f)
            } else {
                all_files <- rbind(all_files, f)
            }
        }
        
        #remove the first blank row
        all_files <- all_files[-1, ]
        
        colnames(all_files) <- c("Query sequence ID", "Subject sequence ID", "Percentage of identical matches", "Alignment length", 
                                 "Number of mismatches", "Number of gap openings", "Start of alignment in query",
                                 "End of alignment in query", "Start of alignment in subject", "End of alignment in subject",
                                 "e-value", "Bit score", "Subject accession", "Subject title")
        
        write.table(all_files, file = "AllBLASTpResults.csv", row.names = FALSE, col.names = TRUE, sep = ",")
    }
    
    

    
    results <- function(){
        cleanResultFiles()
        mergeResultFiles()
        r <- read.table("AllBLASTpResults.csv", header = TRUE, sep = ",", fill = TRUE, stringsAsFactors = FALSE)
        return(r)
    }
    