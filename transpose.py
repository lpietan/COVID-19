#!/usr/bin/env python3

import sys
import pandas as pd

d = pd.read_csv(sys.argv[1],index_col=0)
df1_transposed = d.T
df1_transposed.to_csv(sys.argv[2])
