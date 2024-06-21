## Get arguments
args <- commandArgs(trailingOnly = TRUE)

library(data.table)

covidData <- read.csv(args[1])

covidData <- as.data.frame(covidData)
d1 <- covidData[,1]
covidData <- covidData[,-1]
rownames(covidData) <- d1

covidData_t <- transpose(covidData)
rownames(covidData_t) <- colnames(covidData)
colnames(covidData_t) <- rownames(covidData)

write.csv(covidData_t, args[2], row.names = TRUE)
