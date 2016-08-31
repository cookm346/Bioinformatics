setwd("C:/Users/HIV-Proteomics/Desktop/Matt/HPR")
data <- read.csv("Accession Number Report for Run 1.csv",  header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "  ", "   "))
data <- data[complete.cases(data), ]

h <- HPR(data, 2, 4)
