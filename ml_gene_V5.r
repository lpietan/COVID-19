#!/bin/bash/env Rscript

# "TRAIN.RData" "TEST.RData" "variableKeys_train_dt_rf_X-dataset.csv" "variableKeys_test_dt_rf_X-dataset.csv" "VI_list_prefix"


## Load analysis packages
library(MachineShop)
library(recipes)
library(dplyr)
library(readr)
library(doSNOW)
registerDoSNOW(makeCluster(56))

## Get arguments
args <- commandArgs(trailingOnly = TRUE)

settings(control = function() CVControl(seed = 123))

## Load train data
d_train <- read_csv(args[1])
d_train <- as.data.frame(d_train)
d1 <- d_train[,1]
d_train <- d_train[,-1]
rownames(d_train) <- d1
d_train$Targets <- factor(d_train$Targets)
## these next two lines added
#col_names <- which(grepl("Allele", names(d_train)))
#d_train[,col_names] <- lapply(d_train[,col_names] , factor, levels = c(0,1))
d_train_dt_rf <- d_train
str(d_train) %>% head

d_test <- read_csv(args[2])
d_test <- as.data.frame(d_test)
d1 <- d_test[,1]
d_test <- d_test[,-1]
rownames(d_test) <- d1
d_test$Targets <- factor(d_test$Targets)
#col_names <- which(grepl("Allele", names(d_test)))
#d_test[,col_names] <- lapply(d_test[,col_names] , factor, levels = c(0,1))
d_test_dt_rf <- d_test
str(d_test, list.len=Inf) %>% head

## Fix dataset varibales for tree packages
varNames <- variable.names(d_train_dt_rf)
varNamesL <- length(variable.names(d_train_dt_rf)) - 1
d_train_dt_rf <- setNames(d_train_dt_rf, c("Targets", paste0("Var", 1:varNamesL)))
keyNames <- variable.names(d_train_dt_rf)
varNames_df <- data.frame(varNames)
varNames_df <- setNames(varNames_df, c("variableNames"))
varNames_df$keyNames <- keyNames
write.csv(as.data.frame(varNames_df), args[3], row.names = FALSE)

varNames <- variable.names(d_test_dt_rf)
keys <- read.csv(args[3])
# Iterate through each row in the keys dataframe and rename columns
for (i in 1:nrow(keys)) {
  old_name <- keys[i, "variableNames"]
  new_name <- keys[i, "keyNames"]
  if (old_name %in% names(d_test_dt_rf)) {
    names(d_test_dt_rf)[names(d_test_dt_rf) == old_name] <- new_name
  } else {
    warning(paste("Column", old_name, "not found in main dataframe. Skipping."))
  }
}
keyNames <- variable.names(d_test_dt_rf)
varNames_df <- data.frame(varNames)
varNames_df <- setNames(varNames_df, c("variableNames"))
varNames_df$keyNames <- keyNames
write.csv(as.data.frame(varNames_df), args[4], row.names = FALSE)

file_dt_rds <- paste(args[5], "_dt_modelFit.rds", sep = "")
file_rf_rds <- paste(args[5], "_rf_modelFit.rds", sep = "")
file_lasso_rds <- paste(args[5], "_lasso_modelFit.rds", sep = "")
file_nb_rds <- paste(args[5], "_nb_modelFit.rds", sep = "")
file_svmLin_rds <- paste(args[5], "_svmLin_modelFit.rds", sep = "")
file_svmPoly_rds <- paste(args[5], "_svmPoly_modelFit.rds", sep = "")
file_svmRB_rds <- paste(args[5], "_svmRB_modelFit.rds", sep = "")
file_logReg_rds <- paste(args[5], "_logReg_modelFit.rds", sep = "")
file_xgbt_rds <- paste(args[5], "_xgbt_modelFit.rds", sep = "")
file_enet_rds <- paste(args[5], "_enet_modelFit.rds", sep = "")
file_out_train <- paste(args[5], "_train_results.csv", sep = "")
file_out_test <- paste(args[5], "_test_results.csv", sep = "")
file_dt_csv <- paste(args[5], "_dt_vi.csv", sep = "")
file_dt_plot <- paste(args[5], "_dt_vi.pdf", sep = "")
file_rf_csv <- paste(args[5], "_rf_vi.csv", sep = "")
file_rf_plot <- paste(args[5], "_rf_vi.pdf", sep = "")
file_lasso_csv <- paste(args[5], "_lasso_vi.csv", sep = "")
file_lasso_plot <- paste(args[5], "_lasso_vi.pdf", sep = "")
file_nb_csv <- paste(args[5], "_nb_vi.csv", sep = "")
file_nb_plot <- paste(args[5], "_nb_vi.pdf", sep = "")
file_svmLin_csv <- paste(args[5], "_svmLin_vi.csv", sep = "")
file_svmLin_plot <- paste(args[5], "_svmLin_vi.pdf", sep = "")
file_svmPoly_csv <- paste(args[5], "_svmPoly_vi.csv", sep = "")
file_svmPoly_plot <- paste(args[5], "_svmPoly_vi.pdf", sep = "")
file_svmRB_csv <- paste(args[5], "_svmRB_vi.csv", sep = "")
file_svmRB_plot <- paste(args[5], "_svmRB_vi.pdf", sep = "")
file_logReg_csv <- paste(args[5], "_logReg_vi.csv", sep = "")
file_logReg_plot <- paste(args[5], "_logReg_vi.pdf", sep = "")
file_xgbt_csv <- paste(args[5], "_xgbt_vi.csv", sep = "")
file_xgbt_plot <- paste(args[5], "_xgbt_vi.pdf", sep = "")
file_enet_csv <- paste(args[5], "_enet_vi.csv", sep = "")
file_enet_plot <- paste(args[5], "_enet_vi.pdf", sep = "")

print("DT")
model_dt <- TunedModel(
  TreeModel,
  metrics = c(brier, accuracy, kappa2, roc_auc, sensitivity, specificity),
  grid = expand_params(
    mincut = c(5,6),
    minsize = c(10,15),
    mindev = c(0.01,0.001),
    split = c("gini", "deviance"),
    best = 10),
)
print(model_dt)
ML_fit_dt_full_pipeline <- fit(Targets ~ ., data = d_train_dt_rf, model = model_dt)
summary(ML_fit_dt_full_pipeline)
options(dplyr.width = Inf)
saveRDS(ML_fit_dt_full_pipeline, file_dt_rds)
print("DT training performance")
x <- print(as.MLModel(ML_fit_dt_full_pipeline), n = Inf, width = Inf)
print("DT test performance")
obs_test_dt <- response(ML_fit_dt_full_pipeline, newdata = d_test_dt_rf)
pred_test_prob_dt <- predict(ML_fit_dt_full_pipeline, newdata = d_test_dt_rf, type = "prob")
print(performance(obs_test_dt, pred_test_prob_dt))

x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]
r <- as.numeric(as.vector(r))
final_train_metric_df <- as.data.frame(r[2:7])
rownames(final_train_metric_df) <- c("Brier", "Accuracy", "Kappa", "ROC AUC", "Sensitivity", "Specificity")
colnames(final_train_metric_df) <- "DT"

per <- performance(obs_test_dt, pred_test_prob_dt)
final_test_metric_df <- as.data.frame(per)
colnames(final_test_metric_df) <- "DT"

vi <- varimp(ML_fit_dt_full_pipeline, samples = 25)
pdf(file_dt_plot)
plot(vi, n = 40)
dev.off()
vi_df <- as.data.frame(vi)
write.csv(vi_df, file_dt_csv, row.names=TRUE)


print("RF")
modelRF <- TunedModel(RandomForestModel,grid = 5)
print(modelRF)
ML_fitRF <- fit(Targets ~ ., data = d_train_dt_rf, model = modelRF)
summary(ML_fitRF)
options(dplyr.width = Inf)
saveRDS(ML_fitRF, file_rf_rds)
print("RF training performance")
x <- print(as.MLModel(ML_fitRF), n = Inf, width = Inf)
print("RF test performance")
obs_test <- response(ML_fitRF, newdata = d_test_dt_rf)
pred_test_prob <- predict(ML_fitRF, newdata = d_test_dt_rf, type = "prob")
performance(obs_test, pred_test_prob)

x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]
r <- as.numeric(as.vector(r))
r <- r[2:7]
final_train_metric_df["RF"] <- as.data.frame(r)

per <- performance(obs_test, pred_test_prob)
final_test_metric_df["RF"] <- as.numeric(as.vector(per))

vi <- varimp(ML_fitRF, samples = 25)
pdf(file_rf_plot)
plot(vi, n = 40)
dev.off()
vi_df <- as.data.frame(vi)
write.csv(vi_df, file_rf_csv, row.names=TRUE)


print("lasso")
modelL <- TunedModel(GLMNetModel(alpha = 1), grid = c(lambda = 5))
print(modelL)
ML_fitL <- fit(Targets ~ ., data = d_train, model = modelL)
summary(ML_fitL)
options(dplyr.width = Inf)
saveRDS(ML_fitL, file_lasso_rds)
print("lasso training performance")
x <- print(as.MLModel(ML_fitL), n = Inf, width = Inf)
print("lasso test performance")
obs_test <- response(ML_fitL, newdata = d_test)
pred_test_prob <- predict(ML_fitL, newdata = d_test, type = "prob")
print(performance(obs_test, pred_test_prob))

x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]
r <- as.numeric(as.vector(r))
r <- r[2:7]
final_train_metric_df["lasso"] <- as.data.frame(r)

per <- performance(obs_test, pred_test_prob)
final_test_metric_df["lasso"] <- as.numeric(as.vector(per))

vi <- varimp(ML_fitL, samples = 25)
pdf(file_lasso_plot)
plot(vi, n = 40)
dev.off()
vi_df <- as.data.frame(vi)
write.csv(vi_df, file_lasso_csv, row.names=TRUE)


print("NB")
modelNB <- TunedModel(NaiveBayesModel, grid = expand_params(laplace = 0))
print(modelNB)
ML_fitNB <- fit(Targets ~ ., data = d_train, model = modelNB)
summary(ML_fitNB)
options(dplyr.width = Inf)
saveRDS(ML_fitNB, file_nb_rds)
print("NB training performance")
x <- print(as.MLModel(ML_fitNB), n = Inf, width = Inf)
print("NB test performance")
obs_test <- response(ML_fitNB, newdata = d_test)
pred_test_prob <- predict(ML_fitNB, newdata = d_test, type = "prob")
print(performance(obs_test, pred_test_prob))

x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]
r <- as.numeric(as.vector(r))
r <- r[2:7]
final_train_metric_df["NB"] <- as.data.frame(r)

per <- performance(obs_test, pred_test_prob)
final_test_metric_df["NB"] <- as.numeric(as.vector(per))

vi <- varimp(ML_fitNB, samples = 25)
pdf(file_nb_plot)
plot(vi, n = 40)
dev.off()
vi_df <- as.data.frame(vi)
write.csv(vi_df, file_nb_csv, row.names=TRUE)


print("SVM Linear")
modelSVML <- TunedModel(SVMLinearModel, grid = c("C" = 10))
print(modelSVML)
ML_fitSVML <- fit(Targets ~ ., data = d_train, model = modelSVML)
summary(ML_fitSVML)
options(dplyr.width = Inf)
saveRDS(ML_fitSVML, file_svmLin_rds)
print("SVM linear training performance")
x <- print(as.MLModel(ML_fitSVML), n = Inf, width = Inf)
print("SVM linear test performance")
obs_test <- response(ML_fitSVML, newdata = d_test)
pred_test_prob <- predict(ML_fitSVML, newdata = d_test, type = "prob")
print(performance(obs_test, pred_test_prob))

x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]
r <- as.numeric(as.vector(r))
r <- r[2:7]
final_train_metric_df["SVM Linear"] <- as.data.frame(r)

per <- performance(obs_test, pred_test_prob)
final_test_metric_df["SVM Linear"] <- as.numeric(as.vector(per))

vi <- varimp(ML_fitSVML, samples = 25)
pdf(file_svmLin_plot)
plot(vi, n = 40)
dev.off()
vi_df <- as.data.frame(vi)
write.csv(vi_df, file_svmLin_csv, row.names=TRUE)


print("SVM Poly")
modelSVMP <- TunedModel(SVMPolyModel, grid = expand_params(C = c(0.05,0.25,1,4,16,20), degree = as.integer(c(1,2,3,4,5)), scale = c(0.001,0.0012,0.015,0.17,1,2)))
print(modelSVMP)
ML_fitSVMP <- fit(Targets ~ ., data = d_train, model = modelSVMP)
summary(ML_fitSVMP)
saveRDS(ML_fitSVMP, file_svmPoly_rds)
options(dplyr.width = Inf)
print("SVM poly training performance")
x <- print(as.MLModel(ML_fitSVMP), n = Inf, width = Inf)
print("SVM poly test performance")
obs_test <- response(ML_fitSVMP, newdata = d_test)
pred_test_prob <- predict(ML_fitSVMP, newdata = d_test, type = "prob")
print(performance(obs_test, pred_test_prob))

x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]
r <- as.numeric(as.vector(r))
r <- r[2:7]
final_train_metric_df["SVM Poly"] <- as.data.frame(r)

per <- performance(obs_test, pred_test_prob)
final_test_metric_df["SVM Poly"] <- as.numeric(as.vector(per))

vi <- varimp(ML_fitSVMP, samples = 25)
pdf(file_svmPoly_plot)
plot(vi, n = 40)
dev.off()
vi_df <- as.data.frame(vi)
write.csv(vi_df, file_svmPoly_csv, row.names=TRUE)


print("SVM radial basis")
model_SVM_rad <- TunedModel(SVMRadialModel, grid = expand_params(C = c(0.05,0.25,1,4,16,20), sigma = c(0.01,0.03,0.05,0.07,0.09,0.1)))
print(model_SVM_rad)
ML_fit_SVM_rad <- fit(Targets ~ ., data = d_train, model = model_SVM_rad)
saveRDS(ML_fit_SVM_rad, file_svmRB_rds)
options(dplyr.width = Inf)
print("SVM rad basis training performance")
x <- print(as.MLModel(ML_fit_SVM_rad), n = Inf, width = Inf)
print("SVM rad basis test performance")
obs_test <- response(ML_fit_SVM_rad, newdata = d_test)
pred_test_prob <- predict(ML_fit_SVM_rad, newdata = d_test, type = "prob")
print(performance(obs_test, pred_test_prob))

x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]
r <- as.numeric(as.vector(r))
r <- r[2:7]
final_train_metric_df["SVM radial basis"] <- as.data.frame(r)

per <- performance(obs_test, pred_test_prob)
final_test_metric_df["SVM radial basis"] <- as.numeric(as.vector(per))

vi <- varimp(ML_fit_SVM_rad, samples = 25)
pdf(file_svmRB_plot)
plot(vi, n = 40)
dev.off()
vi_df <- as.data.frame(vi)
write.csv(vi_df, file_svmRB_csv, row.names=TRUE)


print("XGBTreeModel")
model_XGBTree <- TunedModel(XGBTreeModel, grid = 5)
print(model_XGBTree)
ML_fit_XGBTree <- fit(Targets ~ ., data = d_train, model = model_XGBTree)
saveRDS(ML_fit_XGBTree, file_xgbt_rds)
options(dplyr.width = Inf)
print("XGBTree training performance")
x <- print(as.MLModel(ML_fit_XGBTree), n = Inf, width = Inf)
print("XGBTree test performance")
obs_test <- response(ML_fit_XGBTree, newdata = d_test)
pred_test_prob <- predict(ML_fit_XGBTree, newdata = d_test, type = "prob")
print(performance(obs_test, pred_test_prob))

x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]
r <- as.numeric(as.vector(r))
r <- r[2:7]
final_train_metric_df["XGBTreeModel"] <- as.data.frame(r)

per <- performance(obs_test, pred_test_prob)
final_test_metric_df["XGBTreeModel"] <- as.numeric(as.vector(per))

vi <- varimp(ML_fit_XGBTree, samples = 25)
pdf(file_xgbt_plot)
plot(vi, n = 40)
dev.off()
vi_df <- as.data.frame(vi)
write.csv(vi_df, file_xgbt_csv, row.names=TRUE)


print("Elasticnet")
model_net <- TunedModel(GLMNetModel, grid = 5)
print(model_net)
ML_fit_net <- fit(Targets ~ ., data = d_train, model = model_net)
saveRDS(ML_fit_net, file_enet_rds)
options(dplyr.width = Inf)
print("Elasticnet training performance")
x <- print(as.MLModel(ML_fit_net), n = Inf, width = Inf)
print("Elasticnet test performance")
obs_test <- response(ML_fit_net, newdata = d_test)
pred_test_prob <- predict(ML_fit_net, newdata = d_test, type = "prob")
print(performance(obs_test, pred_test_prob))

x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]
r <- as.numeric(as.vector(r))
r <- r[2:7]
final_train_metric_df["Elasticnet"] <- as.data.frame(r)

per <- performance(obs_test, pred_test_prob)
final_test_metric_df["Elasticnet"] <- as.numeric(as.vector(per))

vi <- varimp(ML_fit_net, samples = 25)
pdf(file_enet_plot)
plot(vi, n = 40)
dev.off()
vi_df <- as.data.frame(vi)
write.csv(vi_df, file_enet_csv, row.names=TRUE)


print("Training Metrics (without Log Reg)")
print(final_train_metric_df)
print("Testing Metrics (without Log Reg)")
print(final_test_metric_df)


print("Log Reg")
model_log_reg <- TunedModel(GLMModel)
print(model_log_reg)
ML_fit_log_reg <- fit(Targets ~ ., data = d_train, model = model_log_reg)
saveRDS(ML_fit_log_reg, file_logReg_rds)
options(dplyr.width = Inf)
print("Log reg training performance")
x <- print(as.MLModel(ML_fit_log_reg), n = Inf, width = Inf)
print("Log reg test performance")
obs_test <- response(ML_fit_SVM_rad, newdata = d_test)
pred_test_prob <- predict(ML_fit_SVM_rad, newdata = d_test, type = "prob")
print(performance(obs_test, pred_test_prob))

x <- x@steps
x_ts <- x$TrainingStep1
x_ts_log <- x_ts@log
x_ts_log_sel <- x_ts_log$selected
metrics_df <- as.data.frame(x_ts_log_sel)
x_ts_log_met <- x_ts_log$metrics
metrics_df <- cbind(metrics_df,as.data.frame(x_ts_log_met))
r <- metrics_df[metrics_df$x_ts_log_sel == TRUE,]
r <- as.numeric(as.vector(r))
r <- r[2:7]
final_train_metric_df["Log Reg"] <- as.data.frame(r)

per <- performance(obs_test, pred_test_prob)
final_test_metric_df["Log Reg"] <- as.numeric(as.vector(per))

vi <- varimp(ML_fit_log_reg, samples = 25)
pdf(file_logReg_plot)
plot(vi, n = 40)
dev.off()
vi_df <- as.data.frame(vi)
write.csv(vi_df, file_logReg_csv, row.names=TRUE)


print("Training Metrics")
print(final_train_metric_df)
print("Testing Metrics")
print(final_test_metric_df)


write.csv(final_train_metric_df, file_out_train, row.names=TRUE, quote = FALSE)
write.csv(final_test_metric_df, file_out_test, row.names=TRUE, quote = FALSE)


quit(save="no")
