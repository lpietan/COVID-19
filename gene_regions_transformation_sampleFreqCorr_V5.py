#!/usr/bin/env python3

## arg1 IN FILE
## arg2 OUT FILE
## arg3 bin size
## arg4 correction threshold 0-1 0.5
## arg5 head lines  
## arg6 training index file

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
correction_threshold = float(sys.argv[4])
lineIndex = 0
testSetIndex = []
testSetAlt = []


j = ','
# dataset file
with open(sys.argv[5], "r") as f:
	head_lines = f.readlines()
samSize = len(head_lines[1].split(',')) - 1
samThreshold = samSize * correction_threshold
with open(sys.argv[1], "r") as f:
	lines = f.readlines()
with open(sys.argv[2], "w") as f:
	with open(sys.argv[6], "w") as f2:
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
			count = len([i for i in num_var_list if i > 0])
			if count > samThreshold:
				testSetIndex.append(lineIndex)
				correctedVar = []
				altValue = [i for i in list(set(num_var_list)) if i > 0][0]
				testSetAlt.append(altValue)
				for i in num_var_list:
					if i > 0:
						correctedVar.append(0)
					else:
						correctedVar.append(altValue)
				num_var_list = correctedVar
			lineIndex = lineIndex + 1
			## Start binning
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
		lineIndexListFix = j.join(map(str, testSetIndex)) + "\n"
		lineAltListFix = j.join(map(str, testSetAlt))
		f2.write(lineIndexListFix)
		f2.write(lineAltListFix)




