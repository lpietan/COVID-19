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
## options = list(threshold = 0.05)) line below sets pvalue threshold. Any variable with pvalue below will pass filter. This
## variable can be tuned. If tuning this variable make sure to change all instances of below code to match.
sbf_rec <- recipe(Targets ~ ., data = df_2) %>%
  step_sbf(all_predictors(), filter = glm_filter,
           options = list(threshold = 0.05))

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
           options = list(threshold = 0.05))

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



quit(save="no")
