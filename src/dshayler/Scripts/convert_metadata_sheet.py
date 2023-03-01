#!/usr/bin/python

import pandas as pd
import sys

# argv[1] = metadata sheet in csv (comma delimited)
# Sample_ID,Kit_ID,Seq_Number,Fetal_Age,,Collection_Group,Poor_Read_Number,Moderate_Alignment,Rod_Cells,Non_Photoreceptors,Collection_Method,Outliers
# X313,,1,16,,1,,,,,,
# X314,,1,16,,1,,,,,,

a = pd.read_csv(sys.argv[1], sep=",", index_col=0, dtype=object)
a[a!=a] = "Empty"

for c in a.columns:
	u = a[c].unique()
	for i in u:
		name = c+"_"+str(i)
		print name+"\t"+"\t".join(a.loc[a[c]==i,c].index)
