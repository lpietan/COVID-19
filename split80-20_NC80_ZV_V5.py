#!/usr/bin/env python3

import sys
import pandas as pd


## Split
filename = sys.argv[1]
df = pd.read_csv(filename, dtype=object)

df_test = df.sample(frac=0.2, replace=False, axis=1, random_state=123)
df_train = df.drop(list(df_test.columns), axis = 1)

df_test.to_csv(sys.argv[2])
df_train.to_csv(sys.argv[3])


# NC80
fail = 0.8*len(df_train.columns)
remove = False
j = ','
# dataset file
with open(sys.argv[3], "r") as f:
        lines = f.readlines()
with open(sys.argv[4], "w") as f:
        for line in lines:
                lineN = line.rstrip('\n')
                lineList = lineN.split(',')
                lineCount = lineList.count('.')
                if remove == True:
                        remove = False
                elif lineCount <= fail:
                        if lineCount > 0:
                                ## Change '.' to 0
                                while lineCount > 0:
                                        index = lineList.index('.')
                                        lineList[index] = '0'
                                        lineCount -= 1
                                lineFix = j.join(lineList) + '\n'
                                f.write(lineFix)
                        else:
                                f.write(line)
                else:
                        remove = True


## Zero Variance
# dataset file
with open(sys.argv[4], "r") as f:
        lines = f.readlines()
with open(sys.argv[5], "w") as f:
        for line in lines:
                lineN = line.rstrip('\n')
                lineList = lineN.split(',')
                line_set = set(lineList)
                lineSetLen = len(line_set)
                if lineSetLen > 2:
                        f.write(line)






