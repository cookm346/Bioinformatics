binByGenus <- function(data,  pattern = "OS=", rm.zero.vars = TRUE){

    #generate empty matrix to store genus
    genus <- matrix("", nrow(data))
    
    
    for(i in 1:nrow(data)){
        #entire protein name
        x <- as.character(rownames(data)[i])
        
        #split protein name by words
        ss <- strsplit(x, split = " ")
        ss <- unlist(ss)
        
        #word location of "OS=" identifier
        l <- grep(pattern, ss)
        
        #if the "OS=" identifier isnt in the string,
        # start next iteration of the loop
        if(length(l) == 0){
            next
        }
        
        #remove the "OS=" and
        #update genus matrix with current genus and characters 4+ (i.e., no "OS=")
        genus[i] <- substring(ss[l], (nchar(pattern) + 1))
        
    }
    
    data <- cbind(genus, data)
    
    #change NA to 0
    data[is.na(data)] <- 0 
    
    #sum each variable by genus
    data <- aggregate(data = data, . ~ genus, sum)
    
    #use genus names for column names
    rownames(data) <- data[,1]
    
    #remove genus names variable
    data <- data[, -1]
    
    #remove subjects with only 0s for counts
    if(rm.zero.vars == TRUE){
        x <- as.numeric(apply(data, 2, sum))
        x <- which(x == 0)
        
        #make sure there are actually variables to remove
        #if x was 0, the entire dataset would be removed
        if(length(l) != 0){
            data <- data[ , -(x)]
        }
    }
    
    write.csv(data, "binned-by-genus.csv")
    return(data)
}
