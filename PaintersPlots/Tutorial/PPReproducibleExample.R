##########################################################################
#REPRODUCIBLE EXAMPLE
# FOR PTM FUNCTIONS
#
#INSTRUCTIONS:
# 1.) Run all code in PTMFunctions.R
#       This will load all the functions into your workspace
# 2.) use the code below to explore a verifiable and reproducible example
#########################################################################



#########################################################################
#GENEREATE FAKE DATA
#########################################################################
    data <- matrix(round(runif(100, 1, 1000), 2), 10, 10)
    rownames(data) <- paste("Protein", letters[1:10])
    colnames(data) <- paste("Subject", 1:10)
    
    data   #view data



#########################################################################  
#EXAMPLE USE OF ALL 6 PAINTER'S PLOT FUNCTIONS
#########################################################################

    
    #GET SUMMARY OF DATA
    getSummary(data)                 #default
    getSummary(data, prop = FALSE)   #indetical
    getSummary(data, prop = TRUE)    #if data was already proportional
    
    
    #GENERATE PAINTER'S PLOT
    # install.packages("ggplot2)                           #install ggplot2 package (if not yet installed)
    library(ggplot2)                                       #load required ggplot2 package
    
    PPlot(data)                                            #default
    PPlot(data, 2)                                         #2 branches
    PPlot(data, 2, subject_numbers = FALSE)                #original subjects names (rather than numbers)     
    PPlot(data, 2, omit = "Protein b")                     #omit protein b
    PPlot(data, 2, omit = c("Protein b", "Protein f"))     #omit protein b and f
    
    
    #GENERATE SHANNON'S H BAR PLOT
    # install.packages("vegan")                            #install vegan package (if not yet installed)
    library(vegan)                                         #load vegan package
    
    shanH(data)                                            #default
    shanH(data, dist = "pearson")                          #pearson as opposed to euclidean distance metric
    shanH(data, 2)                                         #2 bracnhes
    
    
    
    
    #THESE FUNCTIONS ARE MAINLY USED AS INTERNAL FUNCTIONS FOR THE PPLOT FUNCTION
    # but may on occasion be used on their own

        #CONVERT DATA TO PROPORTIONS
        getProp(data)
        pd <- getProp(data)  #save proportional data as pd
        colSums(pd)          #verify proportions were computed correctly (columns sum to 1)
        
        
        #CONVERT DATA TO LONG FORMAT
        meltData(data)                             #use default arguments
        meltData(data, n_branches = 2)             #2 branches
        meltData(data, subject_numbers = FALSE)    #use original subject names
        meltData(data, dist = "pearson")           #pearson distance metric
    
        
        #VIEW HIERARCHICAL CLUSTERING DENDROGRAMS
        plotClust(data)                            #use default arguments
        plotClust(data, subject_numbers = FALSE)   #include original subject names (rather than numbers)
    
    
    
    
    
    