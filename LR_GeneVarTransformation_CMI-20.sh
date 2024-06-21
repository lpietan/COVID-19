#!/bin/bash

## How this file will be called from the command line.
## bash LR_GeneVarTransformation_CMI-20.sh "example_ML-dataset_V4.2.csv"

dataset=$1
file=${dataset::-4}


echo "STARTING SPLIT NC80 ZV"
python split80-20_NC80_ZV_V5.py $1 $file"_TEST.csv" $file"_TRAIN.csv" $file"_TRAIN_NC80.csv" $file"_TRAIN_NC80_zv.csv"

echo "STARTING TRANSPOSE"
python transpose.py $file"_TRAIN_NC80_zv.csv" $file"_TRAIN_NC80_zv_T.csv"

echo "STARTING LR"
Rscript LR_gene_V5.r $file"_TRAIN_NC80_zv_T.csv" $file"_TRAIN_NC80_zv_T_LR.csv"

echo "STARTING TRANSPOSE"
python transpose.py $file"_TRAIN_NC80_zv_T_LR.csv" $file"_TRAIN_NC80_zv_T_LR_T.csv"

## Preprocessing for gene variable transformation
tail -n +3 $file"_TRAIN_NC80_zv_T_LR_T.csv" > $file"_TRAIN_NC80_zv_T_LR_T_VARS.csv"
head -n 2 $file"_TRAIN_NC80_zv_T_LR_T.csv" > $file"_TRAIN_NC80_zv_T_LR_T_gene.csv"
sort --field-separator=':' -V -k 3,3 -k 1,1 -k 2,2 $file"_TRAIN_NC80_zv_T_LR_T_VARS.csv" > $file"_TRAIN_NC80_zv_T_LR_T_VARS_SORTED.csv"

echo "STARTING TO CREATE GENE VARIABLES"
## The last (3rd) argument for the script (25000) specifies region bin size (25kb) for variant variables in the same genomic
## region that do not have a gene symbol annotation. 
python gene_regions_transformation_V5.py $file"_TRAIN_NC80_zv_T_LR_T_VARS_SORTED.csv" $file"_TRAIN_NC80_zv_T_LR_T_VARS_SORTED_geneVARS.csv" 25000

## Sample allele frequency correction gene transformation. Comment out code line above and uncomment code line below. 
## The 4th argument specifies the sample allele frequency threshold. Alternative allele count for a variable greater than the
## threshold the variable will be corrected.
## python gene_regions_transformation_sampleFreqCorr_V5.py $file"_TRAIN_NC80_zv_T_LR_T_VARS_SORTED.csv" $file"_TRAIN_NC80_zv_T_LR_T_VARS_SORTED_geneVARS.csv" 25000 0.5 $file"_TRAIN_NC80_zv_T_LR_T_gene.csv" "samCorr_index_corrected.csv"

## Directional correction gene transformation. Comment out code line above and uncomment code line below.
## The 4th argument specifies the ratio of alternative alleles for control samples threshold. If the ratio is less than the
## threshold the variable will be corrected. 
## python gene_regions_transformation_directionalCorr_V5.py $file"_TRAIN_NC80_zv_T_LR_T_VARS_SORTED.csv" $file"_TRAIN_NC80_zv_T_LR_T_VARS_SORTED_geneVARS.csv" 25000 0.5 $file"_TRAIN_NC80_zv_T_LR_T_gene.csv" "dirCorr_index_corrected.csv"

cat $file"_TRAIN_NC80_zv_T_LR_T_VARS_SORTED_geneVARS.csv" >> $file"_TRAIN_NC80_zv_T_LR_T_gene.csv"

echo "STARTING TRANSPOSE"
python transpose.py $file"_TRAIN_NC80_zv_T_LR_T_gene.csv" $file"_TRAIN_NC80_zv_T_LR_T_gene_T.csv"

## Removal of temp files. 
rm $file"_TRAIN_NC80_zv_T_LR_T_VARS.csv"
rm $file"_TRAIN_NC80_zv_T_LR_T_VARS_SORTED.csv"
rm $file"_TRAIN_NC80_zv_T_LR_T_VARS_SORTED_geneVARS.csv"

echo "STARTING GENE CMI-20 FILTERING"
Rscript CMI-20_V5.r $file"_TRAIN_NC80_zv_T_LR_T_gene_T.csv" $file"_TRAIN_NC80_zv_T_LR_T_gene_T_CMI-20.csv"

echo "STARTING VARIABLE SELECTION FOR TEST SET"
Rscript sel_final_vars_test_V5.r $file"_TRAIN_NC80_zv_T_LR.csv" "variable_names_for_test_selection.csv"

python sel_final_vars_test_NCfix.py $file"_TEST.csv" "variable_names_for_test_selection.csv" $file"_TEST_selectedVars.csv" $file"_TEST_selectedVars_NCfix.csv"

## Preprocessing for test set gene variable transformation. 
tail -n +3 $file"_TEST_selectedVars_NCfix.csv" > $file"_TEST_selectedVars_NCfix_VARS.csv"
head -n 2 $file"_TEST_selectedVars_NCfix.csv" > $file"_TEST_selectedVars_NCfix_gene.csv"
sort --field-separator=':' -V -k 3,3 -k 1,1 -k 2,2 $file"_TEST_selectedVars_NCfix_VARS.csv" > $file"_TEST_selectedVars_NCfix_VARS_SORTED.csv"

echo "STARTING TO CREATE GENE VARIABLES"
python gene_regions_transformation_V5.py $file"_TEST_selectedVars_NCfix_VARS_SORTED.csv" $file"_TEST_selectedVars_NCfix_VARS_SORTED_geneVARS.csv" 25000

## Use the two lines of code for using either the sample allele frequency correction or the directional correction. 
## python gene_regions_transformation_testSetCorrection_V5.py $file"_TEST_selectedVars_NCfix_VARS_SORTED.csv" $file"_TEST_selectedVars_NCfix_VARS_SORTED_geneVARS.csv" 25000 "samCorr_index_corrected.csv"
## python gene_regions_transformation_testSetCorrection_V5.py $file"_TEST_selectedVars_NCfix_VARS_SORTED.csv" $file"_TEST_selectedVars_NCfix_VARS_SORTED_geneVARS.csv" 25000 "dirCorr_index_corrected.csv"

cat $file"_TEST_selectedVars_NCfix_VARS_SORTED_geneVARS.csv" >> $file"_TEST_selectedVars_NCfix_gene.csv"

echo "STARTING TRANSPOSE"
python transpose.py $file"_TEST_selectedVars_NCfix_gene.csv" $file"_TEST_selectedVars_NCfix_gene_T.csv"

## Removal of temp files.
rm $file"_TEST_selectedVars_NCfix_VARS.csv"
rm $file"_TEST_selectedVars_NCfix_VARS_SORTED.csv"
rm $file"_TEST_selectedVars_NCfix_VARS_SORTED_geneVARS.csv"

echo "STARTING ML ANALYSIS"
Rscript ml_gene_V5.r $file"_TRAIN_NC80_zv_T_LR_T_gene_T_CMI-20.csv" $file"_TEST_selectedVars_NCfix_gene_T.csv" "variableKeys_ML_train_dt_rf.csv" "variableKeys_ML_test_dt_rf.csv" "LR_geneNoCorrection_CMI-20"


