#!/usr/bin/env python3

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

j = ','
# dataset file
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
		## perform correction here
		if lineVarNameListGene == '':
			if lineVarNameListChr == chromosomes:
				if int(lineVarNameListPos) < binThreshold:
					if 'Allele' in lineVarNameListVar:
						num_var_list = [eval(str(i)) for i in lineListNoName]
						alleleBin.append(num_var_list)
					else:
						num_var_list = [eval(str(i)) for i in lineListNoName]
						CADDBin.append(num_var_list)
				elif len(alleleBin) == 0 and len(CADDBin) == 0:
					if 'Allele' in lineVarNameListVar:
						num_var_list = [eval(str(i)) for i in lineListNoName]
						alleleBin.append(num_var_list)
					else:
						num_var_list = [eval(str(i)) for i in lineListNoName]
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
						num_var_list = [eval(str(i)) for i in lineListNoName]
						alleleBin.append(num_var_list)
					else:
						num_var_list = [eval(str(i)) for i in lineListNoName]
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
					num_var_list = [eval(str(i)) for i in lineListNoName]
					alleleBin.append(num_var_list)
				else:
					num_var_list = [eval(str(i)) for i in lineListNoName]
					CADDBin.append(num_var_list)
				chromosomes = lineVarNameListChr
				binThreshold = math.ceil(int(lineVarNameListPos)/binVar) * binVar + 1
				startPos = binThreshold - binVar
		else:
			if lastGeneVar != '':
				if lineVarNameListGene == lastGeneVar:
					if 'Allele' in lineVarNameListVar:
						num_var_list = [eval(str(i)) for i in lineListNoName]
						alleleBin.append(num_var_list)
					else:
						num_var_list = [eval(str(i)) for i in lineListNoName]
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
						num_var_list = [eval(str(i)) for i in lineListNoName]
						alleleBin.append(num_var_list)
					else:
						num_var_list = [eval(str(i)) for i in lineListNoName]
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
					num_var_list = [eval(str(i)) for i in lineListNoName]
					alleleBin.append(num_var_list)
				else:
					num_var_list = [eval(str(i)) for i in lineListNoName]
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





