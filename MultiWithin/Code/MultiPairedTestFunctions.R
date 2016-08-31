#convert within subject wide format data to long format (required for within-subjects ANOVA)
pairedLong <- function(data, n_groups){
    
    #calculate number of subjects
    n_subjects <- ncol(data) / n_groups
    
    #generate time condition categorical veriable
    cond <- rep(paste("M", 1:n_groups, sep = ""), each = n_subjects)
    
    protein <- rep(rownames(data), each = ncol(data))
    subject <- rep(1:(ncol(data) / n_groups), nrow(data))
    time <- rep(cond, nrow(data))
    scores <- c(t(data))
    
    dl <- data.frame(protein, subject, time, scores)
    
    #make sure subject variable is factor (not numeric)
    dl$subject <- factor(dl$subject)
    
    return(dl)
}
    
    
    
#paired student's t test p value
    pairedStudentsp <- function(a, b){
        t <- t.test(a, b, paired = TRUE)$p.value
        return(t)
    }
    
    
    
#paired ANOVA p value    
    pairedANOVAp <- function(data_long){
        a <- summary(aov(scores ~ time + Error(subject/time), data = data_long))
        a <- as.numeric(unlist(a[[2]])["Pr(>F)1"])
        return(a)
    }

    
    
#paired Wilcoxon signed rank p value
   pairedWilcoxonp <- function(a, b){
        w <- wilcox.test(a, b, paired = TRUE, alternative = "two.sided")$p.value
        return(w)
    }
    
    
   
#friedman test p value
    friedmanp <- function(a, b){
        f <- friedman.test(cbind(a, b))$p.value
        return(f)
    }
    
    

#log2 fold change
    log2fc <- function(a, b){
        l <- mean(log2((b / a)))
        return(l)
    }
    
    
    
#compute requested paired test of significance and return the 
# dataset with the requested test statistics appended to data   
    multiPaired <- function(data, test = "t", n_groups = 2){
        
        #data in long format for within-subject ANOVA
        data_long <- pairedLong(data, n_groups)
        
        if(n_groups == 2){
                
                #vector for storing p values
                p <- matrix(0, nrow(data), 5)
                
                #for each row (protein):
                for(i in 1:nrow(data)){
                    
                    #extract pre and post data
                    pre <- as.numeric(data[i, 1:(ncol(data) / 2)])
                    post <- as.numeric(data[i, ((ncol(data) / 2) + 1):ncol(data)])
                    
                    #students p
                    p[i, 1] <- pairedStudentsp(pre, post)
                    
                    #long data for the ith protein
                    dl <- data_long[data_long$protein %in% as.character(unique(data_long$protein)[i]), ]
                    #paired ANOVA
                    p[i, 2] <- pairedANOVAp(dl)
                    
                    #wilcoxon p
                    p[i, 3] <- pairedWilcoxonp(pre, post)
                    
                    #friedman p
                    p[i, 4] <- friedmanp(pre, post)
                    
                    #average log2 fold change
                    p[i, 5] <- log2fc(pre, post)
                    
                }
                
                #as default only have t test included in p values
                pt <- data.frame(p[ , 1])
                colnames(pt) <- "Student's t p"
                
                #if all tests are chosen, 
                # make this explicit so rest of function can work
                if(any(test == "all")){
                    test <- c("log2fc", "t", "anova", "friedman", "wilcoxon")
                }
                
                #add ANOVA if requested
                if(any(test == "anova")){
                    pt <- cbind(pt, p[ , 2])
                    colnames(pt)[ncol(pt)] <- "ANOVA p"
                }
                
                #add wilcoxon test if requested
                if(any(test == "wilcoxon")){
                    pt <- cbind(pt, p[ , 3])
                    colnames(pt)[ncol(pt)] <- "Wilcoxon p"
                }
                
                #add friedman's test if requested
                if(any(test == "friedman")){
                    pt <- cbind(pt, p[ , 4])
                    colnames(pt)[ncol(pt)] <- "Friedman p"
                }
                
                #add log2 fold change if requested
                if(any(test == "log2fc")){
                    pt <- cbind(pt, p[ , 5])
                    colnames(pt)[ncol(pt)] <- "Mean Log2FC"
                }
                
                #append p values to data
                r <- cbind(pt, data)
                
                #if default t test is not included
                # removet results from output
                if(all((test %in% "t") == FALSE)){
                    r <- r[ , -1]
                }
        
                } else {
            
            p <- matrix(0, nrow(data))
            
            for(i in 1:nrow(data)){
                dl <- data_long[data_long$protein %in% as.character(unique(data_long$protein)[i]), ]
                
                #paired ANOVA
                p[i] <- pairedANOVAp(dl)
            }
            
            r <- cbind(p, data)
            colnames(r)[1] <- "ANOVA p"
        }
        
        return(r)
    }
    
