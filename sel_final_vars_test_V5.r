#!/bin/bash/env Rscript

library(readr)

## Get arguments
args <- commandArgs(trailingOnly = TRUE)

## Load data
d <- read_csv(args[1])
d <- as.data.frame(d)
d1 <- d[,1]
d <- d[,-1]
rownames(d) <- d1

col_names <- colnames(d)
write.csv(col_names, args[2], row.names=FALSE)


print("Obtained Variable list")




quit(save="no")
