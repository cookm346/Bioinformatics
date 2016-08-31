#set working directory
    setwd("C:/Users/HIV-Proteomics/Desktop/Matt/BurgenerLab/GenusBinner/Code/ExampleData")

#read in data    
    data <- read.csv("ExampleData.csv", header = TRUE, row.names = 1)    

#bin data by genus    
    results <- binByGenus(data)
