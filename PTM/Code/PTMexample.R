#################################################
# 
# R Code for Peter McQueen
#
# If you have any questions contact:
# cookm346@myumanitoba.ca
#
# INSTRUCTIONS:
# 1.) run all code in functions.R file
# 2.) protein column should be named 'protein' and 
#     mass differences should be named 'massdiff'
#
#################################################

#install and load requisite graphics package
  # install.packages("ggplot2")
  # library(ggplot2)


#load data from .csv file
  # getwd()   #view current working directory
  # setwd()   #set working directory
  data <- read.csv("3013-1-31_Control_Sample_2.csv", header = TRUE, skip = 0)

#get bins and save over data
  data <- getBins(data)

#get summary of data
  s <- getSummary(data)

#get top bins (default: top 10, 0.1 bin size)
  tb <- topBins(data, 0.1, 10)
#top 25 bins with a size of 0.05  
  # tb <- topBins(data, 0.05, 25)

#get mod histograms
  getModHists(data, tb)
  getModHists(data, tb, type = "hist")   #use "hist", "freq", or "both"
  
#get summary by bins
  getModSummary(data, tb)
  
#subset data by top 10 bins  
  d2 <- reduceByBins(data, tb)
  
#get summary of new data by protein
  s2 <- getProteinSummary(d2)
  
#check the number of proteins in different subsets of data
  nProteins(data)
  nProteins(d2)

#subset data by massdiff range (inclusive upper bound)
  d3 <- reduceByRange(data, 17.996, 18.933)

#specify minimum number of counts for
# protein to be included (e.g., 200)
  d4 <- reduceByCounts(data, 200)

#get summary of new data
  s3 <- getProteinSummary(d4)

#generate massdiff histograms for each protein
  getProteinHists(d4)
  
#subset data by specific protein
  d5 <- reduceByProtein(data, "sp|K2C1_HUMAN|,sp|P04264|K2C1_HUMAN")
  getProteinHists(d5)
  
  
  
  
  
