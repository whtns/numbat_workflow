#!/usr/bin/python 

import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns
#~ import ipdb
#~ ipdb.set_trace()

fpkm_file = sys.argv[1]
cutoffs = map(float, sys.argv[2].split(","))
cutoffs.sort(reverse=True)
fpkm = pd.read_csv(fpkm_file, sep="\t", index_col=0)

genes_detected = pd.DataFrame(0, columns = fpkm.columns, index=cutoffs)
upper_bound = sys.float_info.max
for c in cutoffs:
	genes_detected.loc[c] = fpkm[(fpkm >= c) & (fpkm < upper_bound)].count()
	upper_bound = c

cell_order = fpkm.columns.values

ax = genes_detected.transpose().plot.bar(stacked=True, width=0.8)
plt.ylabel('Genes detected')
plt.show()
