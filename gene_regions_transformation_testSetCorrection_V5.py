#!/usr/bin/env python3

## arg 1 IN FILE
## arg 2 OUT FILE
## arg 3 bin size
## arg 4 index csv from train

import sys
import pandas as pd
import math
import csv


binVar = int(sys.argv[3]) ##25000
binThreshold = binVar + 1
lastGeneVar = ''
chromosomes = 'chr1'
startPos = 1
alleleBin = []
CADDBin = []
lineIndex = 0

j = ','
# dataset file
with open(sys.argv[4], "r") as f:
	index_alt_lines = f.readlines()
	index_line = [eval(i) for i in index_alt_lines[0].split(",")]
	alt_line = [eval(i) for i in index_alt_lines[1].split(",")]
with open(sys.argv[1], "r") as f:
	lines = f.readlines()
with open(sys.argv[2], "w") as f:
	for line in lines:
		lineN = line.rstrip('\n')
		lineList = lineN.split(',')
		lineVarName = lineList[0]
		lineVarNameList = lineVarName.split(':')
		lineVarNameListGene = lineVarNameList[2]
		lineVarNameListChr = lineVarNameList[0]
		lineVarNameListPos = lineVarNameList[1]
		lineVarNameListVar = lineVarNameList[3]
		lineListNoName = lineList[1:]
		num_var_list = [eval(str(i)) for i in lineListNoName]
		## perform correction here
		if lineIndex in index_line:
			correctedVar = []
			if sum(num_var_list) > 0:
				altValue = [i for i in list(set(num_var_list)) if i > 0][0]
			else:
				altValue = alt_line[index_line.index(lineIndex)]
			for i in num_var_list:
				if i > 0:
					correctedVar.append(0)
				else:
					correctedVar.append(altValue)
			num_var_list = correctedVar
		lineIndex = lineIndex + 1
		if lineVarNameListGene == '':
			if lineVarNameListChr == chromosomes:
				if int(lineVarNameListPos) < binThreshold:
					if 'Allele' in lineVarNameListVar:
						alleleBin.append(num_var_list)
					else:
						CADDBin.append(num_var_list)
				elif len(alleleBin) == 0 and len(CADDBin) == 0:
					if 'Allele' in lineVarNameListVar:
						alleleBin.append(num_var_list)
					else:
						CADDBin.append(num_var_list)
					binThreshold = math.ceil(int(lineVarNameListPos)/binVar) * binVar + 1
					startPos = binThreshold - binVar
				else:
					if len(alleleBin) > 0:
						varNameAllele = chromosomes + ":" + str(startPos) + "_" + str(binThreshold - 1) + "_Allele,"
						alleleDF = pd.DataFrame(alleleBin)
						finAlleleVar = alleleDF.sum(axis = 0).tolist()
						lineFix = varNameAllele + j.join(map(str, finAlleleVar)) + '\n'
						f.write(lineFix)
						alleleBin = []
					if len(CADDBin) > 0:
						varNameCADD = chromosomes + ":" + str(startPos) + "_" + str(binThreshold - 1) + "_CADD,"
						CADDDF = pd.DataFrame(CADDBin)
						finCADDVar = CADDDF.sum(axis = 0).tolist()
						lineFix = varNameCADD + j.join(map(str, finCADDVar)) + '\n'
						f.write(lineFix)
						CADDBin = []
					if 'Allele' in lineVarNameListVar:
						alleleBin.append(num_var_list)
					else:
						CADDBin.append(num_var_list)
					binThreshold = math.ceil(int(lineVarNameListPos)/binVar) * binVar + 1
					startPos = binThreshold - binVar
			else:
				if len(alleleBin) > 0:
					varNameAllele = chromosomes + ":" + str(startPos) + "_" + str(binThreshold - 1) + "_Allele,"
					alleleDF = pd.DataFrame(alleleBin)
					finAlleleVar = alleleDF.sum(axis = 0).tolist()
					lineFix = varNameAllele + j.join(map(str, finAlleleVar)) + '\n'
					f.write(lineFix)
					alleleBin = []
				if len(CADDBin) > 0:
					varNameCADD = chromosomes + ":" + str(startPos) + "_" + str(binThreshold - 1) + "_CADD,"
					CADDDF = pd.DataFrame(CADDBin)
					finCADDVar = CADDDF.sum(axis = 0).tolist()
					lineFix = varNameCADD + j.join(map(str, finCADDVar)) + '\n'
					f.write(lineFix)
					CADDBin = []
				if 'Allele' in lineVarNameListVar:
					alleleBin.append(num_var_list)
				else:
					CADDBin.append(num_var_list)
				chromosomes = lineVarNameListChr
				binThreshold = math.ceil(int(lineVarNameListPos)/binVar) * binVar + 1
				startPos = binThreshold - binVar
		else:
			if lastGeneVar != '':
				if lineVarNameListGene == lastGeneVar:
					if 'Allele' in lineVarNameListVar:
						alleleBin.append(num_var_list)
					else:
						CADDBin.append(num_var_list)
				else:
					if len(alleleBin) > 0:
						varNameAllele = lastGeneVar + "_Allele,"
						alleleDF = pd.DataFrame(alleleBin)
						finAlleleVar = alleleDF.sum(axis = 0).tolist()
						lineFix = varNameAllele + j.join(map(str, finAlleleVar)) + '\n'
						f.write(lineFix)
						alleleBin = []
					if len(CADDBin) > 0:
						varNameCADD = lastGeneVar + "_CADD,"
						CADDDF = pd.DataFrame(CADDBin)
						finCADDVar = CADDDF.sum(axis = 0).tolist()
						lineFix = varNameCADD + j.join(map(str, finCADDVar)) + '\n'
						f.write(lineFix)
						CADDBin = []
					if 'Allele' in lineVarNameListVar:
						alleleBin.append(num_var_list)
					else:
						CADDBin.append(num_var_list)
					lastGeneVar = lineVarNameListGene
			else:
				if len(alleleBin) > 0:
					varNameAllele = chromosomes + ":" + str(startPos) + "_" + str(binThreshold - 1) + "_Allele,"
					alleleDF = pd.DataFrame(alleleBin)
					finAlleleVar = alleleDF.sum(axis = 0).tolist()
					lineFix = varNameAllele + j.join(map(str, finAlleleVar)) + '\n'
					f.write(lineFix)
					alleleBin = []
				if len(CADDBin) > 0:
					varNameCADD = chromosomes + ":" + str(startPos) + "_" + str(binThreshold - 1) + "_CADD,"
					CADDDF = pd.DataFrame(CADDBin)
					finCADDVar = CADDDF.sum(axis = 0).tolist()
					lineFix = varNameCADD + j.join(map(str, finCADDVar)) + '\n'
					f.write(lineFix)
					CADDBin = []
				if 'Allele' in lineVarNameListVar:
					alleleBin.append(num_var_list)
				else:
					CADDBin.append(num_var_list)
				lastGeneVar = lineVarNameListGene
	if len(alleleBin) > 0:
		varNameAllele = lineVarNameListGene + "_Allele,"
		alleleDF = pd.DataFrame(alleleBin)
		finAlleleVar = alleleDF.sum(axis = 0).tolist()
		lineFix = varNameAllele + j.join(map(str, finAlleleVar)) + '\n'
		f.write(lineFix)
	if len(CADDBin) > 0:
		varNameCADD = lineVarNameListGene + "_CADD,"
		CADDDF = pd.DataFrame(CADDBin)
		finCADDVar = CADDDF.sum(axis = 0).tolist()
		lineFix = varNameCADD + j.join(map(str, finCADDVar)) + '\n'
		f.write(lineFix)





