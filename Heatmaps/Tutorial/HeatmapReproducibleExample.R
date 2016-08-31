####################################################################
####################################################################
#
# Example use of HMFunctions.R code for generating heatmaps
#
# National Lab for HIV Immunology
# University of Manitoba
#
# Instructions:
#   Run all code in HMFunctions.R then Import an included test 
#   dataset or your own
#
# If you have questions contact Matt Cook (cookm346@myumanitoba.ca)
#
####################################################################
####################################################################

#install and load requisite package
  # install.packages("heatmap.plus")
  # library(heatmap.plus)
  # install.packages("TeachingDemos") #to choose location of color key
  # library("TeachingDemos")


#import data from .csv file
  d <- read.csv("33.csv", header = FALSE, row.names = 1, skip = 0, stringsAsFactors = FALSE)
  nrh <- 3  #number of row headers
  nch <- 3  #number of column headers
  
  
#generate heatmap using heat function
  # (data, n row headers, n col headers)
  h <- heat(d, nrh, nch)
  
  
  #add legend
  par(xpd=TRUE)        #allows you to add legend outside of plot area
  coords <- locator(1) #allows you to selectlocation for legend
  
  #example legend
  legend(coords, legend=c("Placebo", "Tenofovir", " ", "Case", "Control"),
         fill=c("purple", "yellow", "white", "green", "darkgreen"),
         border=FALSE, bty="n", y.intersp = 1, cex = 0.8)
  

  #add color key
  coords <- locator(1) #allows you to selectlocation for color key
  subplot(color.bar(colorRampPalette(c("blue", "white", "red"))(256), 
          min = -5, nticks = 5), coords$x, coords$y, size=c(.5, 1), vadj=1, hadj=0, pars=NULL)

  
#save reordered data from heatmap
  # (data, heatmap, n row headers, n col headers)
  rd <- rData(d, h, nrh, nch)
  
  
#conduct Fisher's exact test for count data
# (data, n row headers, n col headers, n branches, column sidebar to test)
  exactTest(d, nrh, nch, 2, 1)

#verify fisher's exact test
# (data, n row headers, n col headers, n branches, column sidebar to test, table)
  exactVerify(d, nrh, nch, 2, 1, FALSE) #branches for column sidebar
  exactVerify(d, nrh, nch, 2, 1, TRUE)  #contingency table
  
  
  