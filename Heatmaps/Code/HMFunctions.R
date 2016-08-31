####################################################################
####################################################################
#
# Heatmap functions
#
# National Lab for HIV Immunology
# University of Manitoba
#
# If you have questions contact Matt Cook (cookm346@myumanitoba.ca)
#
####################################################################
####################################################################

#####################
#extract subject row
#####################

  subject <- function(data, n_row_headers = 0, n_col_headers = 0){
    return(data[1, 1:(ncol(data) - n_row_headers)])
  }

  
########################
#extract column headers
########################
  
  colHeaders <- function(data, n_row_headers = 0, n_col_headers = 0){
    return(data[1:(n_col_headers + 1), 1:(ncol(data) - n_row_headers)])
  }

  
#####################
#extract row headers
#####################
  
  rowHeaders <- function(data, n_row_headers = 0, n_col_headers = 0){
    return(data[(n_col_headers + 2):nrow(data), (ncol(data) - n_row_headers + 1):ncol(data)])
  }

###################### 
#extract heatmap data
###################### 
  
  heatData <- function(data, n_row_headers = 0, n_col_headers = 0){
    heatdata <- data[(n_col_headers + 2):nrow(data), 1:(ncol(data) - n_row_headers)]
    heatdata <- as.matrix(heatdata)
    colnames(heatdata) <- data[1, (1:(ncol(data) - n_row_headers))]
    return(heatdata)
  }
  
  
################################
#generate column sidebar colors
################################
  
  colColors <- function(data, n_row_headers = 0, n_col_headers = 0){
    
    #use colHeaders function to extract column headers from dataset
    d <- data.matrix(colHeaders(data, n_row_headers, n_col_headers))
    
    #default colors
    col_colors <- c("green", "darkgreen", "yellow", "purple", "darkblue", "lemonchiffon",
                    "maroon", "honeydew", "ivory1", "gold", "cyan",
                    colors())
    
    #color counting variable to keep track of the colors used
    cc <- 1
    
    #this loop starts at 2 to exclude the subject variable row
    for(i in 2:nrow(d)){
      u <- unique(d[i, ])      #extract the unique values from the current row
      ul <- length(u)          #the number of those unique values
      dr <- as.matrix(d[i, ])  #the current row from the col header
      
      for(j in 1:ul){
        dr[dr %in% u[j]] <- col_colors[cc] #change dummy variables in data to colors
        
        #print legend details for column sidebar
        if(j == 1 & i == 2) cat("Column Sidebar legend: \n")
        cat(paste(u[j], ": ", col_colors[cc], "\n", sep = ""))
        
        cc <- cc + 1                       #update count of colors
      }
      d[i, ] <- dr                         #put that colored row back in column header
    }
    
    d <- d[-1, ] #remove top (subject) row
    d <- t(d)    #transpose column header (requirement for heatmap.plus package)
    
    #a requirement for the heatmap.plus package is that the column header columns
    # be a matrix, not a vector
    # if the number of column headers is 1, this doubles up that vector to form a 
    # matrix so the heatmap.plus code works
    if(n_col_headers == 1){
      d <- c(d, d)
      d <- matrix(d, (length(d) / 2), 2)
      colnames(d) <- c(" " ,rownames(data)[2])
    }
    
    return(as.matrix(d))
  }

  
#############################
#generate row sidebar colors
#############################
  
  rowColors <- function(data, n_row_headers = 0, n_col_headers = 0){
    
    #use rowHeaders function to get row headers
    d <- data.matrix(rowHeaders(data, n_row_headers, n_col_headers))
    
    #default row header sidebar colors
    row_colors <- c("magenta", "peru", "aliceblue", "darkorange", "lavenderblush", "lawngreen",
                    "lightblue1", "chocolate", "black",
                    colors())
    
    #color counting variable to keep track of the colors used
    cc <- 1
    
    #a requirement for the heatmap.plus package is that the row header columns
    # be a matrix, not a vector
    # if the number of row headers is 1, this doubles up that vector to form a 
    # matrix so the heatmap.plus code works
    if(n_row_headers == 1){
      d <- cbind(d, d)
    }
    
    #transpose the data so that the same basic code can be used from the colColors function
    d <- t(d)
    
    #identical as in the colColors function
    for(i in 1:nrow(d)){
      u <- unique(d[i, ])
      ul <- length(u)
      dr <- as.matrix(d[i, ])
        for(j in 1:ul){
          
          dr[dr %in% u[j]] <- row_colors[cc]
          
          #print legend details for column sidebar
          if(j == 1 & i == 1) cat("Row Sidebar legend: \n")
          cat(paste(u[j], ": ", row_colors[cc], "\n", sep = ""))
          if(j == ul & i == nrow(d)) cat("\n")
          
          cc <- cc + 1
        }
      d[i, ] <- dr
    }
    
    ##transpose the data back
    d <- t(d)

    #add column names
    # if number of row headers is 1, this is done differently because the 
    # two columns of that matrix are identical
    if(n_row_headers == 1){
      d <- cbind(d[ , 1], d[ , 1])
      colnames(d) <- c(data[(n_col_headers + 1), (ncol(data) - n_row_headers + 1):ncol(data)], "")
    } else {
      colnames(d) <- data[(n_col_headers + 1), (ncol(data) - n_row_headers + 1):ncol(data)]
    }
    return(as.matrix(d))
  }
    

######################
#reorder heatmap data  
######################  
  
  rHeatData <- function(heatdata, heatmap){
    return(heatdata[rev(heatmap$rowInd), heatmap$colInd])
  }
  

#####################
#reorder col headers
#####################
  
  rColHeaders <- function(data, heatmap, n_row_headers = 0, n_col_headers = 0){
    #use colHeaders function to get column headers
    ch <- colHeaders(data, n_row_headers, n_col_headers)
    
    #reorder each row of the data in the same order as the heatmap order
    for(i in 1:nrow(ch)){
      ch[i, ] <- ch[i, ][heatmap$colInd]
    }
    return(ch)
  }
  
##################### 
#reorder row headers
#####################  
  
  rRowHeaders <- function(data, heatmap, n_row_headers = 0, n_col_headers = 0){
    
    #use rowHeaders function to get row headers
    rh <- rowHeaders(data, n_row_headers, n_col_headers)
    
    if(n_row_headers == 1){
      rh <- as.numeric(rh)[rev(heatmap$rowInd)]
    } else {
      for(j in 1:ncol(rh)){
        rh[ , j] <- rh[, j][rev(heatmap$rowInd)]
      }
      colnames(rh) <- data[(n_col_headers + 1), (ncol(data) - n_row_headers + 1):ncol(data)]
    }
    return(rh)
  }

#################################
#reorder entire original dataset 
#################################
  
  rData <- function(data, heatmap, n_row_headers = 0, n_col_headers = 0){
    
    #copy over data to new object to be used as reordered data
    rd <- data
    
    #add reordered heatmap data
    rd[(n_col_headers + 2):nrow(data), 1:(ncol(data) - n_row_headers)] <- rHeatData(heatData(data, n_row_headers, n_col_headers), heatmap)
    
    #if column headers exist add the reorderd column headers
    if(n_col_headers > 0){
      rd[1:(n_col_headers + 1), 1:(ncol(data) - n_row_headers)] <- rColHeaders(data, heatmap, n_row_headers, n_col_headers)
    }
    
    #add reorderd subject IDs
    if(n_col_headers == 0){
        rd[1, 1:(ncol(rd) - n_row_headers)] <- subject(data, n_row_headers, n_col_headers)[heatmap$colInd]
    }
    if(n_row_headers > 0){
      rd[(n_col_headers + 2):nrow(data), (ncol(data) - n_row_headers + 1):ncol(data)] <- rRowHeaders(data, heatmap, n_row_headers, n_col_headers)
    }
    
    #add reorderd rownames
    rn <- rownames(heatData(data, n_row_headers, n_col_headers))
    rn <- rn[rev(heatmap$rowInd)]
    rownames(rd)[(n_col_headers + 2):nrow(rd)] <- rn
    return(rd)
  }
  
  
##################
#generate heatmap
##################
  
  heat <- function(data, n_row_headers = 0, n_col_headers = 0){
    
    #number of row and column headers
    nrh <- n_row_headers
    nch <- n_col_headers
    
    #if row headers exist, run rowColors function
    if(nrh > 0){
      rc <- rowColors(data, nrh, nch)
    }
    
    #if column headers exist, run rowColors function
    if(nch > 0){
      cc <- colColors(data, nrh, nch)
      cc <- apply(cc, 1, rev)
      cc <- t(cc)
      
      if(nch == 1){
        colnames(cc) <- rev(colnames(cc))
      }
    }
    
    #extract heatmap data
    hd <- heatData(data, nrh, nch)
    
    #the code below runs one of four heatmap.plus codes 
    # depending on whether row and column sidebars exist
    
    if(nrh == 0 & nch == 0){
      h <- heatmap.plus(hd, 
                        margins = c(5, 17), Rowv = NULL, Colv = NULL, cexRow = 0.8, cexCol = 0.8, 
                        col = colorRampPalette(c("blue", "white", "red"))(256))
    }
    
    if(nrh == 0 & nch > 0){
      h <- heatmap.plus(hd,
                        ColSideColors = cc, 
                        margins = c(5, 17), Rowv = NULL, Colv = NULL, cexRow = 0.8, cexCol = 0.8, 
                        col = colorRampPalette(c("blue", "white", "red"))(256))
    }
    
    if(nrh > 0 & nch == 0){
      h <- heatmap.plus(hd,
                        RowSideColors = rc, 
                        margins = c(5, 17), Rowv = NULL, Colv = NULL, cexRow = 0.8, cexCol = 0.8, 
                        col = colorRampPalette(c("blue", "white", "red"))(256))
    }
    
    if(nrh > 0 & nch > 0){
      h <- heatmap.plus(hd,
                        RowSideColors = rc, 
                        ColSideColors = cc, 
                        margins = c(5, 17), Rowv = NULL, Colv = NULL, cexRow = 0.8, cexCol = 0.8, 
                        col = colorRampPalette(c("blue", "white", "red"))(256))
    }
    
    return(h)
  }
  
  
  
############################################
#conduct Fisher's exact test for count data
############################################  
  
  exactTest <- function(data, n_row_headers = 0, n_col_headers = 0, n_branches, col_header){
    
    #heatmap data
    hd <- heatData(data, n_row_headers, n_col_headers)
    
    #hierachical clustering using euclidean distance metric
    hc <- hclust(dist(t(hd)))
    
    #this uses the cutree function to determine the branch that each column header belongs to
    ft <- data.frame(as.numeric(cutree(hc, n_branches)), as.numeric(data[(col_header + 1), 1:(ncol(data) - n_row_headers)]))
    colnames(ft) <- c("Branch", rownames(data)[(col_header + 1)])
    
    #this summarizes the above
    m <- table(ft)
    
    #conduct the fisher's exact test
    if(n_branches == 2){
      f <- fisher.test(m, alternative = "two.sided", hybrid = FALSE)
    } else {
      f <- fisher.test(m, alternative = "two.sided", hybrid = TRUE)
    }
    return(f)    
  }

  
############################################
#Verify Fisher's exact test for count data
############################################  
  
  exactVerify <- function(data, n_row_headers = 0, n_col_headers = 0, n_branches, col_header, table = FALSE){
    
    #heatmap data
    hd <- heatData(data, n_row_headers, n_col_headers)
    
    #hierachical clustering using euclidean distance metric
    hc <- hclust(dist(t(hd)))
    
    #this uses the cutree function to determine the branch that each column header belongs to
    ft <- data.frame(as.numeric(cutree(hc, n_branches)), as.numeric(data[(col_header + 1), 1:(ncol(data) - n_row_headers)]))
    colnames(ft) <- c("Branch", rownames(data)[(col_header + 1)])
    rownames(ft) <- 1:nrow(ft)
    
    #this summarizes the above
    m <- table(ft)
    
    if(table == FALSE){
      return(t(ft))
    } else {
      return(table(ft))
    }
    
  }
  
##############################
# Function to plot color bar
##############################
  
  color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    
    #dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
  }

  
