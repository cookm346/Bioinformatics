multiBetween <- function(data){
    
    #vector to store p values
    p <- matrix(0, (nrow(data) - 1))
    
    #conduct all between subject anaovas
    for(i in 2:(nrow(data))){
        #conduct ANOVA
        a <- summary(aov(data[i, ] ~ data[1, ]))
        #extract p value and update p values vector
        p[i - 1] <- as.numeric(unlist(a[[1]])["Pr(>F)1"])
    }
    
    #append p values to original dataset and add column names
    r <- data.frame(rbind(NA , p), data)
    colnames(r) <- c("ANOVA p", colnames(data))
    
    #return results: original dataset with the ANOVA p values appended
    return(r)
}
