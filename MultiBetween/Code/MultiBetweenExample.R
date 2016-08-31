#generate some fake data with condition row variable
    d <- matrix(rnorm(100), 10, 10)
    cond <- rep(c("control", "treatment"), each = 5)
    data <- rbind(cond, d)
    rm(d, cond)
    rownames(data)[2:11] <- paste("Microbe", letters[1:10])
    colnames(data) <- paste("Subject", 1:10)
    
    
#conduct ANOVA for each microbe
    results <- multiBetween(data)
    