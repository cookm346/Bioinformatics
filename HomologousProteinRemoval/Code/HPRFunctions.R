extractNames <- function(data, procol, pattern = "OS="){
    
    #vector to store row locations where no genus and species name are available
    mis <- 0
    
    gen_spe <- matrix("", nrow(data), 2)
    
    for(i in 1:nrow(data)){
        
        #current protein name
        x <- as.character(data[i, procol])
        
        #the protein name split into words
        ss <- strsplit(x, split = " ")
        ss <- unlist(ss)
        
        #the word lcoation for the start of "OS="
        l <- grep(pattern, ss)
        
        if(length(l) == 0){
            
            mis <- c(mis, i)
            
            next
        } else {
            
            #genus    
            genus <- ss[l]
            genus <- substring(genus, 4)
            
            #species
            species <- ss[l + 1]
            
            #add genus and species names to list of names
            gen_spe[i, 1] <- genus
            gen_spe[i, 2] <- species
            
        }
    }
    
    #remove the first element (which is 0) from the rows with no genus or species info
    mis <- mis[-1]
    
    # gen_spe <- paste(gen_spe[,1], gen_spe[,2])
    
    r <- cbind(data, gen_spe[,1])
    r <- cbind(r, gen_spe[,2])
    colnames(r)[(ncol(r) - 1):ncol(r)] <- c("Genus", "Species")
    return(r)
}



HPR <- function(data, procol, grocol, pattern = "OS="){
    
    #run the extract names function
    data <- extractNames(data, procol, pattern)
    
    #order data by protein
    data <- data[order(data[ , grocol]),] 
    
    #extract the grouping column from the data
    protein <- data[ , grocol]
    
    #extract the genus and species column from the data
    genus <- data[ ,(ncol(data) - 1)]
    species <- data[ ,ncol(data)]
    
    
    #results matrices
    result_genus <- matrix("", nrow(data))
    result_species <- matrix("", nrow(data))
    protein_count <- matrix(1, nrow(data))
    
    #set i to 1 for while loop
    i <- 1
    
    while(i <= length(protein)){
        
        #if there is only one grouping protein 
        # add one to the while loop counter
        if (length(which(protein == protein[i])) == 1){
            i <- i + 1
            
        #if there are multiple grouping proteins that are the same
        } else {
            # find the min and max row locations
            min <- min(which(protein == protein[i]))
            max <- max(which(protein == protein[i]))
            
            #if the genus or species are not uniue
            # add an asterisk to the results column
            if(length(unique(genus[min:max])) != 1){
                result_genus[min:max] <- "*"
            }
            if(length(unique(species[min:max])) != 1){
                result_species[min:max] <- "*"
            }
            #then, add to while loop counter
            i <- i + length(max + 1)
        }
    }
    
    #generate unique counts for proteins
    
    #unique protein counter
    count <- 1
    
    #loop through the list of proteins groups and 
    # check if the current protein group is identical to the previous protein group
    for(i in 2:length(protein)){
        #if the current and previous protein groups are identical
        # add the current count to the current protein
        if(protein[i] == protein[i - 1]){
            protein_count[i] <- count
        } else {
            #update the count by one
            # and add that new count
            count <- count + 1
            protein_count[i] <- count
        }
    }
    
    #save results matrix and save to csv file
    results <- data.frame(data, protein_count, result_genus, result_species)
    write.csv(results, "output.csv")
    return(results)
}
