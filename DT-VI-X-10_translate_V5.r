#!/bin/bash/env Rscript

# arg1 input data
# arg2 variable keys for dt and rf
# arg3 file name for plot pdf
# arg4 file name for filtered dataset
# arg5 file name for translated dataset
# arg6 number of variable to sample at each iteration for intermidiate datasets

## Load analysis packages
library(MachineShop)
library(recipes)
library(dplyr)
library(readr)

## Register cores for parallel processing
library(doSNOW)
registerDoSNOW(makeCluster(56))


## Get arguments
args <- commandArgs(trailingOnly = TRUE)

## Load data
d <- read_csv(args[1])
d <- as.data.frame(d)
d1 <- d[,1]
d <- d[,-1]
rownames(d) <- d1

d$Targets <- factor(d$Targets)

## Uncomment out line below for testing with the allele variable factored (as categorical variables)
#col_names <- which(grepl("Allele", names(d)))
#d[,col_names] <- lapply(d[,col_names] , factor)

h <- length(d)
print("Number of variables in starting dataset")
print(h)



## Start DT-VI-Iter filter
print("Starting DT-VI-Iter Filtering")

d$Targets <- factor(d$Targets)
print("  *************************   DATASET SUMMARY   *************************  ")
str(d) %>% head


## Fix dataset varibales for tree package
varNames <- variable.names(d)
varNamesL <- length(variable.names(d)) - 1
d <- setNames(d, c("Targets", paste0("Var", 1:varNamesL)))
keyNames <- variable.names(d)
varNames_df <- data.frame(varNames)
varNames_df <- setNames(varNames_df, c("variableNames"))
varNames_df$keyNames <- keyNames
write.csv(as.data.frame(varNames_df),args[2], row.names = FALSE)

t_df <- d[1]
varImp_dataframe <- t_df
targets <- t_df$Targets
i <- 2:length(d)
print("Initial index list Complete")

## Majority Class comparison
## This sets threshold of model to perform better than a mjority class model. majCl can be set to any desired threshold.
sampleNumber <- length(t_df$Targets)
sampleNumberhalf = sampleNumber/2
if (sum(t_df==1) >= sampleNumberhalf){
majorityClass <- sum(t_df==1)
} else {
majorityClass <- sum(t_df==0)
}
majCl <- majorityClass/sampleNumber

model_dt <- TunedModel(
  TreeModel,
  metrics = c(accuracy, brier, kappa2, roc_auc, sensitivity, specificity),
  control = CVControl(seed = 123),
  grid = expand_params(
    mincut = c(5,6),
    minsize = c(10,15),
    mindev = c(0.01,0.001),
    split = c("gini", "deviance"),
    best = 10)
)
print(model_dt)

iterations <- c(1,2,3,4,5,6,7,8,9,10)
number_of_variables <- c()

## Resample control (10-fold, defualt)
settings(control = function() CVControl(seed = 123))

whileLoopSampleNumber <- as.numeric(args[6]) - 1

## Change seed for different random intermidiate dataset selections.
sseed <- 122
## Number of iterations corresponds to DT-VI-1000-X
for (iter in 1:10) 
{
print(iter)
sseed <- sseed+iter

while (length(i) > whileLoopSampleNumber)
{
print(length(i))
set.seed(sseed)
variables_selected <- sample(x=i, size=as.numeric(args[6]))
df_2 <- d[variables_selected]
df_2$Targets <- targets


tryCatch(
expr = {
## Model fit
ML_fit <- fit(Targets ~ ., data = df_2, model = model_dt)
x <- as.MLModel(ML_fit)
x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]

perf_Acc <- r[["accuracy"]]

## if state to perform vi and keep variables 
if (perf_Acc > majCl) {
## relative importance of included predictor variables
vi <- varimp(ML_fit)
colnames(vi) <- c('score')
selectVI <- vi[1] > 0
selectedVI <- as.data.frame(selectVI) %>% filter(score == TRUE)
selectedVariables <- row.names(selectedVI)

newSelectedVariables = setdiff(selectedVariables, colnames(varImp_dataframe))

VI_dataset <- select(df_2, all_of(newSelectedVariables))

## Add selected variable to final dataset
varImp_dataframe <- cbind(varImp_dataframe, VI_dataset)
}

## Remove selected variable from index list for next iteration
#i <- i[! i %in% variables_selected]
},
error = function(e){
print("Error in larger dataset, breaking down dataset to proceed")
variables_selected_1 <- variables_selected[1:500]
variables_selected_2 <- variables_selected[(501):1000]
df_2 <- d[variables_selected_1]
df_2$Targets <- targets
ML_fit <- fit(Targets ~ ., data = df_2, model = model_dt)
x <- as.MLModel(ML_fit)
x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]

## selection by accuracy metric, can be set to other metrics here. If changed make sure to change below code to match.
perf_Acc <- r[["accuracy"]]

## if state to perform vi and keep variables 
if (perf_Acc > majCl) {
## relative importance of included predictor variables
vi <- varimp(ML_fit)
colnames(vi) <- c('score')
selectVI <- vi[1] > 0
selectedVI <- as.data.frame(selectVI) %>% filter(score == TRUE)
selectedVariables <- row.names(selectedVI)

newSelectedVariables = setdiff(selectedVariables, colnames(varImp_dataframe))

VI_dataset <- select(df_2, all_of(newSelectedVariables))

## Add selected variable to final dataset
varImp_dataframe <- cbind(varImp_dataframe, VI_dataset)
}
print("Successfully ran breakdown model 1")
df_2 <- d[variables_selected_2]
df_2$Targets <- targets
ML_fit <- fit(Targets ~ ., data = df_2, model = model_dt)
x <- as.MLModel(ML_fit)
x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]

perf_Acc <- r[["accuracy"]]

## if state to perform vi and keep variables 
if (perf_Acc > majCl) {
## relative importance of included predictor variables
vi <- varimp(ML_fit)
colnames(vi) <- c('score')
selectVI <- vi[1] > 0
selectedVI <- as.data.frame(selectVI) %>% filter(score == TRUE)
selectedVariables <- row.names(selectedVI)

newSelectedVariables = setdiff(selectedVariables, colnames(varImp_dataframe))

VI_dataset <- select(df_2, all_of(newSelectedVariables))

## Add selected variable to final dataset
varImp_dataframe <- cbind(varImp_dataframe, VI_dataset)
}
print("Successfully ran breakdown model 2")
#i <- i[! i %in% variables_selected]
print("Error handled successfully, moving on to next iteration")
},
finally = {
i <- i[! i %in% variables_selected]
}
)
}

print(length(i))
if (length(i) > 0){
variables_selected <- i
df_2 <- d[variables_selected]
df_2$Targets <- targets

## Model fit
ML_fit <- fit(Targets ~ ., data = df_2, model = model_dt)
x <- as.MLModel(ML_fit)
x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]

perf_Acc <- r[["accuracy"]]

## if state to perform vi and keep variables 
if (perf_Acc > majCl) {
## relative importance of included predictor variables
vi <- varimp(ML_fit)
colnames(vi) <- c('score')
selectVI <- vi[1] > 0
selectedVI <- as.data.frame(selectVI) %>% filter(score == TRUE)
selectedVariables <- row.names(selectedVI)

newSelectedVariables = setdiff(selectedVariables, colnames(varImp_dataframe))

VI_dataset <- select(df_2, all_of(newSelectedVariables))


## Add selected variable to final dataset
varImp_dataframe <- cbind(varImp_dataframe, VI_dataset)


}
}

number_of_variables <- append(number_of_variables, length(varImp_dataframe))

i <- 2:length(d)

}

# Initial trouble shooting plotting of variables filtered through
print("number_of_variables")
print(number_of_variables)
pdf(args[3])
plot(x = iterations, y = number_of_variables, type = "o")
points(x = iterations, y = number_of_variables, pch = 16)
dev.off()


## Save results
d <- as.data.frame(varImp_dataframe)
write.csv(d, args[4], row.names = TRUE)

h <- length(d)
print("Number of selected variables")
print(h)

print("Variable Importance Filtering Complete")





## Start translating variables
print("Starting to translate variables")

d$Targets <- factor(d$Targets)

keys <- read.csv(args[2])

selected <- as.data.frame(colnames(d))

varSel <- c()

for(i in 1:nrow(selected)) {
  value <- keys[keys$keyNames==selected[i, ],]$variableNames
  varSel <- append(varSel, value)
}

colnames(d) <- c(varSel)

df_final <- as.data.frame(d)
write.csv(df_final, args[5], row.names = TRUE)




quit(save="no")
