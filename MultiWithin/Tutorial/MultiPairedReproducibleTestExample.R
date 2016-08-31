#generate some fake data
    n_subjects <- 7  #change the number of subjects
    n_proteins <- 8  #change the number of proteins
    data <- matrix(runif((n_subjects * 2 * n_proteins), 1, 100), n_proteins, (n_subjects * 2))
    data <- round(data, 4)
    rownames(data) <- paste("Protein", 1:n_proteins)
    colnames(data) <- c(paste("Subject", 1:n_subjects, "pre"), paste("Subject", 1:n_subjects, "post"))
    
    

#get test statistic for each protein (log2fc, t, ANOVA, Friedman, Wilcoxon, or all)
# examples:
    results <- multiPaired(data)                       #paired student's t
    results <- multiPaired(data, "t")                  #identical
    results <- multiPaired(data, "anova")              #identical (Student's t = one-way ANOVA)
    results <- multiPaired(data, "log2fc")             #mean log2 fold change
    results <- multiPaired(data, "friedman")           #Friedman's test
    results <- multiPaired(data, c("t", "wilcoxon"))   #paired student's t and wilcoxon
    results <- multiPaired(data, "all")                #all tests
    