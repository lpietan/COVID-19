#!/bin/bash/env Rscript

# Modified due to not having to combine that truncated outputs of CMI 

## Load analysis packages
library(MachineShop)
library(recipes)
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



## Start GLM filter
print("Starting GLM Filtering")
print("  *************************   DATASET SUMMARY   *************************  ")
str(d) %>% head


t_df <- d[1]
UniGLM_dataframe <- t_df
targets <- t_df$Targets
i <- 2:length(d)
print("Initial index list Complete")

## GLM filter for step_sbf(), chi-squared statistic
glm_filter <- function(x, y, step) {
  model_fit <- glm(y ~ ., family = binomial, data = data.frame(y, x))
  p_value <- drop1(model_fit, test = "Chisq")[-1, "Pr(>Chi)"]
  p_value < step$threshold
}

## For testing F test statistic for LR filter. Comment out above glm_filter function and uncomment below function. 
#glm_filter <- function(x, y, step) {
#  model_fit <- glm(y ~ ., family = binomial, data = data.frame(y, x))
#  p_value <- drop1(model_fit, test = "F")[-1, "Pr(>F)"]
#  p_value < step$threshold
#}


while (length(i) > 999)
{
print(length(i))
variables_selected <- sample(x=i, size=1000)
df_2 <- d[variables_selected]
df_2$Targets <- targets


## Recipe with custom univariate SBF step
## options = list(threshold = 0.10)) line below sets pvalue threshold. Any variable with pvalue below will pass filter. This
## variable can be tuned. If tuning this variable make sure to change all instances of below code to match.
sbf_rec <- recipe(Targets ~ ., data = df_2) %>%
  step_sbf(all_predictors(), filter = glm_filter,
           options = list(threshold = 0.10))

## Trained recipe
sbf_prep <- prep(sbf_rec)
#tidy(sbf_prep, number = 1)

## Applied recipe
f <- bake(sbf_prep, df_2)
f <- f[1:length(f) - 1]

## Add CMI selected variable to final dataset
UniGLM_dataframe <- cbind(UniGLM_dataframe, f)

## Remove selected variable from index list for next iteration
i <- i[! i %in% variables_selected]
}

print(length(i))
variables_selected <- i
df_2 <- d[variables_selected]
df_2$Targets <- targets

## Recipe with custom univariate SBF step
sbf_rec <- recipe(Targets ~ ., data = df_2) %>%
  step_sbf(all_predictors(), filter = glm_filter,
           options = list(threshold = 0.10))

## Trained recipe
sbf_prep <- prep(sbf_rec)
#tidy(sbf_prep, number = 1)

## Applied recipe
f <- bake(sbf_prep, df_2)
f <- f[1:length(f) - 1]

## Add CMI selected variable to final dataset
UniGLM_dataframe <- cbind(UniGLM_dataframe, f)

## Save results
d <- as.data.frame(UniGLM_dataframe)
write.csv(d, args[2], row.names = TRUE)

h <- length(d)
print("Number of selected variables")
print(h)

print("Univariate GLM Filtering Complete")





## Start DT-VI filter
print("Starting DT-VI Filtering")

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
write.csv(as.data.frame(varNames_df),args[3], row.names = FALSE)

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


## Resample control (10-fold, defualt)
settings(control = function() CVControl(seed = 123))

## Change seed for different random intermidiate dataset selections.
sseed <- 123

## Change 999 and size=1000 for a different number of randomly selected variables for intermidiate datasets. If tuning this variable make sure to change below
## code to match.
while (length(i) > 999)
{
print(length(i))
variables_selected <- sample(x=i, size=1000)
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

VI_dataset <- select(df_2, all_of(selectedVariables))

## Add selected variable to final dataset
varImp_dataframe <- cbind(varImp_dataframe, VI_dataset)
}

## Remove selected variable from index list for next iteration
i <- i[! i %in% variables_selected]
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

VI_dataset <- select(df_2, all_of(selectedVariables))

## Add selected variable to final dataset
varImp_dataframe <- cbind(varImp_dataframe, VI_dataset)
}
}

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

keys <- read.csv(args[3])

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
