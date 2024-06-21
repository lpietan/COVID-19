#!/bin/bash

## How this file will be called from the command line.
## bash CMI-5_LR_DT-VI-1000-1.sh "example_ML-dataset_V4.2.csv"

dataset=$1
file=${dataset::-4}


echo "STARTING SPLIT NC80 ZV"
python split80-20_NC80_ZV_V5.py $1 $file"_TEST.csv" $file"_TRAIN.csv" $file"_TRAIN_NC80.csv" $file"_TRAIN_NC80_zv.csv"

echo "STARTING TRANSPOSE"
python transpose.py $file"_TRAIN_NC80_zv.csv" $file"_TRAIN_NC80_zv_T.csv"

echo "STARTING CMI"
Rscript CMI-5_V5.r $file"_TRAIN_NC80_zv_T.csv" $file"_TRAIN_NC80_zv_T_CMI-5.csv"

echo "STARTING LR DT-VI translate"
Rscript LR_DT-VI-1000-1_translate_V5.r $file"_TRAIN_NC80_zv_T_CMI-5.csv" $file"_TRAIN_NC80_zv_T_CMI-5_LR.csv" "variableKeys_DT-VI.csv" $file"_TRAIN_NC80_zv_T_CMI-5_LR_DT-VI-1000-1.csv" $file"_TRAIN_NC80_zv_T_CMI-5_LR_DT-VI-1000-1_Translated.csv"

echo "STARTING VARIABLE SELECTION FOR TEST SET"
Rscript sel_final_vars_test_V5.r $file"_TRAIN_NC80_zv_T_CMI-5_LR_DT-VI-1000-1_Translated.csv" "variable_names_for_test_selection.csv"

python sel_final_vars_test_NCfix.py $file"_TEST.csv" "variable_names_for_test_selection.csv" $file"_TEST_selectedVars.csv" $file"_TEST_selectedVars_NCfix.csv"

Rscript test_transpose_V5.r $file"_TEST_selectedVars_NCfix.csv" $file"_TEST_selectedVars_NCfix_T.csv"

echo "STARTING ML ANALYSIS"
## The last (5th) argument for the script (CMI-5_LR_DT-VI-1000-1) sets the prefix for the output file names.
Rscript ml_V5.r $file"_TRAIN_NC80_zv_T_CMI-5_LR_DT-VI-1000-1_Translated.csv" $file"_TEST_selectedVars_NCfix_T.csv" "variableKeys_ML_train_dt_rf.csv" "variableKeys_ML_test_dt_rf.csv" "CMI-5_LR_DT-VI-1000-1"

