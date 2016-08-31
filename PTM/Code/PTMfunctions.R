#################################################
#                   FUNCTIONS                   #
#################################################

#################################################
# EXTRACT NUMBER OF PROTEINS
#################################################

  nProteins <- function(data){
      return(length(unique(data$protein)))
  }
  
  
#################################################
# GENERATE SUMMARY OF DATA BY PROTEIN
#################################################

  getProteinSummary <- function(data){
      
      #explanation of aggregate function arguments:
      # aggreate(variable to summarize, by factor, summary statistic)
      # this function is used once for each summary statistic
      #the output of this function is a data frame with the names of the
      # proteins in the first column, and mass difference summary statistic
      # in the second
      # after the first summary statistic (count/length), the [2] index is
      # used to only add the mass difference summary statistic to the data frame
      x <- aggregate(data$massdiff, list(data$protein), length)
      x <- cbind(x, aggregate(data$massdiff, list(data$protein), mean)[2])
      x <- cbind(x, aggregate(data$massdiff, list(data$protein), median)[2])
      x <- cbind(x, aggregate(data$massdiff, list(data$protein), min)[2])
      x <- cbind(x, aggregate(data$massdiff, list(data$protein), max)[2])
      x <- cbind(x, aggregate(data$massdiff, list(data$protein), sd)[2])
      colnames(x) <- c("Protein", "Count", "Mean", "Median", "Min", "Max", "SD")
      return(x)
  }

  
#################################################
# SUBSET DATA BY PROTEIN COUNT
#################################################

  reduceByCounts <- function(data, thres = 1){
      #get a table of number of occurances by protein
      counts <- table(data$protein)
      #return only those proteins that are greater than or equal to
      # the threshold
      return(data[data$protein %in% names(counts[counts >= thres]), ])
  }
 
   
#################################################
# SUBSET DATA BY MASS DIFFERENCE TOP BINS
#################################################

  reduceByBins <- function(data, topbins){
    #return data, but only those proteins that appear in the topbins
    return(data[data$`Bin Center` %in% topbins$`Bin Center`, ])

  }


#################################################
# SUBSET DATA BY MASS DIFFERENCE RANGE
#################################################

  reduceByRange <- function(data, min, max){
    #return data, but only those that are greater than or equal to
    # the minmum specified and less than or equal to the maximum specified
    return(data[which(data$massdiff >= min & data$massdiff <= max), ])
  }
  
  
#################################################
# SUBSET DATA BY PROTEIN
#################################################
  
  reduceByProtein <- function(data, protein){
    #return data, but only those protein names that are specified
    return(data[data$protein %in% protein , ])
  }

  
#################################################
# DETERMINE BINS FOR EACH MASS DIFFERENCE
#################################################

  getBins <- function(data, size = 0.1){
    
    #generate a sequence from the lowest mass difference (rounded down) 
    # to the highest mass difference (rounded up)
    # by the size of the bin multiplied by 2 (for both sides of the bin)
    b <- seq(from = floor(min(data$massdiff)), to = ceiling(max(data$massdiff)), by = (size * 2))
    
    #if 0 is contained in the sequcen of bins generated
    # generate a new sequence with more extreme lower and upper bounds, 
    # and shift that new sequence to contain 0
    #this is done to ensure 0 is contained in the sequence of bins because 0 mass difference
    # will be the most common bin
    if(0 %in% b){
      b <- seq(from = (floor(min(data$massdiff)) - 1), to = (ceiling(max(data$massdiff)) + 1), by = (size * 2))
      b <- b - size
    }
    
    #run .bincode function to get the bins of all mass differences
    bin_results <- .bincode(data$massdiff, b, FALSE, TRUE)
    #return only the sequence of bins that occur in the results
    bins <- b[bin_results]
    #bind the data with the reults of the corresponding mass difference bin
    d <- cbind(data, (bins + size))
    #add some column names
    colnames(d) <- c(colnames(data), paste("Bin Center", "(+/- ", size , ")", sep = ""))
    return(d)
  }


#################################################
# DERIVE TOP MASS DIFFERENCE BINS
#################################################

  topBins <- function(data, size = 0.1, n = 10){
      
      #generate a sequence from the lowest mass difference (rounded down) 
      # to the highest mass difference (rounded up)
      # by the size of the bin multiplied by 2 (for both sides of the bin)
      b <- seq(from = floor(min(data$massdiff)), to = ceiling(max(data$massdiff)), by = (size * 2))
      
      
      #if 0 is contained in the sequcen of bins generated
      # generate a new sequence with more extreme lower and upper bounds, 
      # and shift that new sequence to contain 0
      #this is done to ensure 0 is contained in the sequence of bins because 0 mass difference
      # will be the most common bin
      if(0 %in% b){
        b <- seq(from = (floor(min(data$massdiff)) - 1), to = (ceiling(max(data$massdiff)) + 1), by = (size * 2))
        b <- b - size
      }
      
      #run .bincode function to get the bins of all mass differences
      bin_results <- .bincode(data$massdiff, b, FALSE, TRUE)
      #return only the sequence of bins that occur in the results
      bins <- b[bin_results]
      
      #get table of mass difference by bins
      cc <- aggregate(data$massdiff, list(bins), length)
      #order the bins from highest counts to lowest
      cc <- cc[order(-cc[ , 2]),]
      #return only the top n bins
      cc <- cc[1:n , ]
      #add the size of the bin (center to one edge) to the bins 
      # so that the center of the bin is displayed rather than the lowest bound
      cc[ , 1] <- cc[ , 1] + size
      
      #add column and row names
      colnames(cc) <- c(paste("Bin Center", "(+/- ", size , ")", sep = ""), "Count")
      rownames(cc) <- 1:n
      return(cc)
      
  }


#################################################
# GENERATE HISTOGRAMS BY MASS DIFFERENCE TOP BINS
#################################################

  getModHists <- function(data, topbins, type = "hist"){
    #make variables a little more readable
    d <- data
    tb <- topbins
    
    #bins in data that appear in the top bins
    nd <- d[d$`Bin Center` %in% tb$`Bin Center`, ]
    colnames(nd)[ncol(nd)] <- "b" 
    
    #the three if statements below just change the ggplot2 code depending on whether the user wants
    # a histogram, density curve, or both
    if(type == "hist"){
      g <- ggplot(nd, aes(x = massdiff)) + geom_histogram(binwidth = .01, colour = "black", fill = "white") + 
        facet_wrap(~ b, scales = "free")
    }
    
    if(type == "freq"){
      g <- ggplot(nd, aes(x = massdiff)) + geom_freqpoly(binwidth = .01, colour = "red") + 
        facet_wrap(~ b, scales = "free")
    }
    if(type == "both"){
      g <- ggplot(nd, aes(x = massdiff)) + geom_histogram(binwidth = .01, colour = "black", fill = "white") +
        geom_freqpoly(binwidth = .01, colour = "red") + 
        facet_wrap(~ b, scales = "free")
    }
    return(g)
  }


#################################################
# GENERATE HISTOGRAMS BY PROTEIN
#################################################

  getProteinHists <- function(data, type = "hist"){
    
    #the three if statements below just change the ggplot2 code depending on whether the user wants
    # a histogram, density curve, or both
    if(type == "hist"){
    g <- ggplot(data, aes(x = massdiff)) + geom_histogram(binwidth = .01, colour = "black", fill = "white") + 
             facet_wrap(~ protein, scales = "free")
    }
    
    if(type == "freq"){
      g <- ggplot(data, aes(x = massdiff)) + geom_freqpoly(binwidth = .01, colour = "red") + 
               facet_wrap(~ protein, scales = "free")
    }
    if(type == "both"){
      g <- ggplot(data, aes(x = massdiff)) + geom_histogram(binwidth = .01, colour = "black", fill = "white") +
               geom_freqpoly(binwidth = .01, colour = "red") + 
               facet_wrap(~ protein, scales = "free")
    }
    return(g)
  }

  
#################################################
# GENERATE SUMMARY OF DATA BY BIN
#################################################
  
  getModSummary <- function(data, topbins){
    
    #explanation of aggregate function arguments:
    # aggreate(variable to summarize, by factor, summary statistic)
    # this function is used once for each summary statistic
    #the output of this function is a data frame with the names of the
    # proteins in the first column, and mass difference summary statistic
    # in the second
    # after the first summary statistic (count/length), the [2] index is
    # used to only add the mass difference summary statistic to the data frame
    #unlike the getProteinSUmmary function, this function only searches the proteins
    # that fall within the massdifference ranges of the top bins
    
    tb <- topbins
    data <- data[data$`Bin Center` %in% tb$`Bin Center`, ]
    
    x <- aggregate(data$massdiff, list(data[ , ncol(data)]), length)
    x <- cbind(x, aggregate(data$massdiff, list(data[ , ncol(data)]), mean)[2])
    x <- cbind(x, aggregate(data$massdiff, list(data[ , ncol(data)]), median)[2])
    x <- cbind(x, aggregate(data$massdiff, list(data[ , ncol(data)]), min)[2])
    x <- cbind(x, aggregate(data$massdiff, list(data[ , ncol(data)]), max)[2])
    x <- cbind(x, aggregate(data$massdiff, list(data[ , ncol(data)]), sd)[2])
    colnames(x) <- c("Bin", "Count", "Mean", "Median", "Min", "Max", "SD")
    return(x)
  }
  
  
