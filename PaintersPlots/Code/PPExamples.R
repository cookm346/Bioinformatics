####################################################################
####################################################################
#
# Example use of PPFunctions.R code for generating Painter's Plots
#
# National Lab for HIV Immunology
# University of Manitoba
#
# Instructions:
#   Run all code in PPFunctions.R then Import an included test 
#   dataset or your own
#
# If you have questions contact Matt Cook (cookm346@myumanitoba.ca)
#
####################################################################
####################################################################

#install and load requisite packages
# install.packages("ggplot2")         #for graphics
# library(ggplot2)
# install.packages("vegan")           #for Shannon's H function
# library(vegan)

#load data
    setwd("C:/Users/HIV-Proteomics/Desktop/Matt/HESN")
    d <- read.csv("MSM_curated%20microbiome%20dataRAW.csv", header = TRUE, row.names = 1, skip = 0, stringsAsFactors = FALSE)


#generate painters plot
  PPlot(d, 5)                     #default is euclidean
  PPlot(d, 5, dist = "pearson")   #pearson
  PPlot(d, 5, dist = "spearman")  #spearman
  
  #Change colors of bugs
  PPlot(d, 5, colors = c("red", "green", "blue"))
  
  #omit certain bugs by name
  PPlot(d, 5, omit = "Saccharomyces")

  
#view hierachical clustering of data
  plotClust(d)
  
  
#summary of dataset by bug
  getSummary(d)
  
  
#generate Shannon's H barplot
  shanH(d, 5, omit = "Saccharomyces")

  
  