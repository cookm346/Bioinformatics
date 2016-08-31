##########################################################################
#REPRODUCIBLE EXAMPLE
# FOR PTM FUNCTIONS
#
#INSTRUCTIONS:
# 1.) Run all code in PTMFunctions.R
#       This will make all these PTM functions usable to you
# 2.) Use the code below to explore a verifiable, reproducible example
#########################################################################



#########################################################################
#GENEREATE FAKE DATA
#########################################################################
    #simple dataset with 3 unique proteins each occuring 2, 3, and 3 times
    # and associated mass differences between -100 and 100
    data <- data.frame(protein = rep(paste("Protein", letters[1:3]), c(2,3,3)), massdiff = c(0, 15.5, 0, 15.7, 39.3, 0, 15.9, 39.3))
    
    data    #view data
    
    #critically, the protein variable is labelled "protein"
    # and the mass difference variable "massdiff"
    # these are necessary variables names for the functions to work
    
    
    
#########################################################################  
#EXAMPLE USE OF ALL 11 PTM FUNCTIONS
#########################################################################
    
    #DISPLAY NUMBER OF UNIQUE PROTEINS
    nProteins(data)
    
    
    
    #DISYPLAY SUMMARY STATISTICS OF MASS DIFFERENCES FOR EACH PROTEIN
        #notice protein 1 only appears twice (e.g., count = 2)
    getProteinSummary(data)
    
    
    
    #DISPLAY HISTOGRAMS OF MASS DIFFERENCES FOR EACH PROTEIN
    # install.packages("ggplot2")  #install ggplot2 if not installed
    library(ggplot2)               #load the ggplot2 package
    getProteinHists(data)          #not very exciting with this dataset
    
    
    
    #BIN MASS DIFFERENCES
    getBins(data) #notice the new variable added to the dataset
    
    
    
    #SUBSET DATA BY MINIMUM COUNTS
    reduceByCounts(data, thres = 3) #only display proteins with counts >= 3
    reduceByCounts(data, 3)         #identical as above
    
    
    
    #SUBSET DATA BY RANGE OF MASS DIFFERENCE (INCLUSIVE)
    reduceByRange(data, min = 15, max = 16) #proteins with mass differences between 15 and 16
    reduceByRange(data, 15, 16)             #identical
    
    
    
    #SUBSET BY PROTEIN NAME
    reduceByProtein(data, "Protein a")                   #only display protein a data
    reduceByProtein(data, c("Protein a", "Protein c"))   #protein a and c
    
    
    
    #GET TOP 2 BINS (with total bin size of 0.2)
    topBins(data, size = 0.1, n = 2)
    topBins(data, 0.1, 2)            #identical as above
    tb <- topBins(data, 0.1, 2)      #save above as "tb"
    tb                               #look at "tb"
    
    
    
    #DISPLAY HISTOGRAMS FOR TOP 2 BINS
    data <- getBins(data)          #first, add bins to data
    getModHists(data, tb)          #get histograms for top 2 bins
    getModHists(data, tb, "hist")  #identical
    getModHists(data, tb, "freq")  #frequency (density) plot
    getModHists(data, tb, "both")  #histogram with overlaid frequency
    
    
    #ONLY RETAIN DATA FROM TOP BINS
    reduceByBins(data, tb)    #only data from bin center 0, and 39.2
    
    
    
    #GET MODIFICATION (BIN) SUMMARY
    getModSummary(data, tb)   #summary stats for each mod (bin)
    

    
#########################################################################    
#ADVANCED USE
#########################################################################
    
    #these functions can be used in conjuction with one another
    #for example:
    reduceByCounts(data, thres = 3)        #reduce data by counts (>= 3)
    rd <- reduceByCounts(data, thres = 3)  #save as "rd"
    rd                                     #view reduced data
    
    reduceByRange(rd, 15, 40)              #reduce rd (15 >= mass diff <= 40)
    rd <- reduceByRange(rd, 15, 40)        #overwrite "rd
    rd                                     #view reduced data
    
    reduceByProtein(rd, "Protein b")       #only proteins that occur 3 or mores times
                                           # with mass differences between 14 and 40
    
    
    
    #you could also perform the exact same command by nesting those three functions:
    reduceByProtein(reduceByRange(reduceByCounts(data, thres = 3), 15, 40), "Protein b")
    
    