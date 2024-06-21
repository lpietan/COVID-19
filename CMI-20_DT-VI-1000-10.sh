#!/bin/bash


## How this file will be called from the command line.
## bash CMI-20_DT-VI-1000-10.sh "example_ML-dataset_V4.2.csv"


dataset=$1
file=${dataset::-4}



echo "STARTING SPLIT NC80 ZV"
python split80-20_NC80_ZV_V5.py $1 $file"_TEST.csv" $file"_TRAIN.csv" $file"_TRAIN_NC80.csv" $file"_TRAIN_NC80_zv.csv"

echo "STARTING TRANSPOSE"
python transpose.py $file"_TRAIN_NC80_zv.csv" $file"_TRAIN_NC80_zv_T.csv"

echo "STARTING CMI"
Rscript CMI-20_V5.r $file"_TRAIN_NC80_zv_T.csv" $file"_TRAIN_NC80_zv_T_CMI-20.csv"

echo "STARTING DT-VI translate"
## The last (6th) argument for the script (1000) determines the number of random variables to select for each intermidiate
## dataset. If the number is changed the user would also change the naming of the files (e.g. $file"_CMI-20_DT-VI-1000-10.csv" to
## $file"_CMI-20_DT-VI-500-10.csv").
Rscript DT-VI-X-10_translate_V5.r $file"_TRAIN_NC80_zv_T_CMI-20.csv" "variableKeys_DT-VI-1000-10.csv" "DT-VI-1000-10_plot.pdf" $file"_TRAIN_NC80_zv_T_CMI-20_DT-VI-1000-10.csv" $file"_TRAIN_NC80_zv_T_CMI-20_DT-VI-1000-10_Translated.csv" 1000

echo "STARTING VARIABLE SELECTION FOR TEST SET"
Rscript sel_final_vars_test_V5.r $file"_TRAIN_NC80_zv_T_CMI-20_DT-VI-1000-10_Translated.csv" "var_names_for_test_selection.csv"

python sel_final_vars_test_NCfix.py $file"_TEST.csv" "var_names_for_test_selection.csv" $file"_TEST_selectedVars.csv" $file"_TEST_selectedVars_NCfix.csv"

Rscript test_transpose_V5.r $file"_TEST_selectedVars_NCfix.csv" $file"_TEST_selectedVars_NCfix_T.csv"

echo "STARTING ML ANALYSIS"
## The last (5th) argument for the script (CMI-20_DT-VI-1000-10) sets the prefix for the output file names.
Rscript ml_V5.r $file"_TRAIN_NC80_zv_T_CMI-20_DT-VI-1000-10_Translated.csv" $file"_TEST_selectedVars_NCfix_T.csv" "variableKeys_ML_train_dt_rf.csv" "variableKeys_ML_test_dt_rf.csv" "CMI-20_DT-VI-1000-10"

