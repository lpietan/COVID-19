#!/usr/bin/env python3

import sys
import pandas as pd

d = pd.read_csv(sys.argv[1],index_col=0)
d_var = pd.read_csv(sys.argv[2],index_col=0)
d_var = d_var.index.tolist()

d = d.loc[d_var]
d.to_csv(sys.argv[3])

j = ','
# dataset file
with open(sys.argv[3], "r") as f:
        lines = f.readlines()
with open(sys.argv[4], "w") as f:
        for line in lines:
                lineN = line.rstrip('\n')
                lineList = lineN.split(',')
                lineCount = lineList.count('.')
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



