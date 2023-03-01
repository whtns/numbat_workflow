#!/usr/bin/python

import pandas as pd
import os
import sys
import IPython

in_csv = "/home/pwiner/bulk_single.txt"
output_table = "/home/pwiner/output_table"
p_qs = "/home/pwiner/chromosome_arm"
centrof = "/dataVolume/storage/single_cell_pipeline/output/GT_20180109_SHL_H_sapiens_RB_31_output/copywriter/cen.txt"

smeans = pd.read_csv(in_csv, sep = "\t")

carm = pd.read_csv(p_qs, sep = "\t")
carm['band'] = carm['arm']
carm['arm'] = carm.band.str[:1]

centro = pd.read_csv(centrof, sep="\t", header = None, names = ["chrom", "start", "end", "band", "cen"])
centro = centro[centro['band'].str.contains("p")]
centro.chrom=centro.chrom.str.replace("chr", "")
centro.chrom=centro.chrom.str.replace("X", "23")
centro.chrom=centro.chrom.str.replace("Y", "24")

gb = smeans.groupby('ID')
lsmeans = [gb.get_group(x) for x in gb.groups]

names_means = smeans.columns.tolist()
names_means[names_means.index('loc.end')] = 'end'
names_means[names_means.index('loc.start')] = 'start'
names_means[names_means.index('seg.mean')] = 'seg_mean'
names_means[names_means.index('num.mark')] = 'num_mark'
smeans.columns = names_means

arm={}
for ind, row in centro.iterrows():
	for i, r in smeans.iterrows():
		if (r.chrom == int(row.chrom) and r.end<row.start):
			arm[i]='p'
		elif (r.chrom == int(row.chrom) and r.start > row.end):
			arm[i]='q'
		elif (r.chrom == int(row.chrom) and r.start < row.end and r.end > row.start):
			arm[i]='pq'

smeans['arm'] = smeans.index.to_series().map(arm)
smeans.ID=smeans.ID.str.replace("_R1.*", "").str.replace("log2.DR_Seq_", "")
smeans.ID=smeans.ID.str.replace(".CL_1_.*", "")
pq_region=smeans[smeans.arm=="pq"]

p_region=pq_region.copy()
p_region.arm=p_region.arm.str.replace("pq", "p")

q_region = pq_region.copy()
q_region.arm=q_region.arm.str.replace("pq", "q")

joined_table=pd.concat([smeans,p_region,q_region])
final_table=joined_table[joined_table.arm!='pq']

final_table.to_csv("final_table.csv", sep="\t")

#duplicated=[]
#for ind,row in smeans.iterrows():
#	if (row.arm=='pq'):
#		copied_rows=row.copy(ind)
#		print copied_rows
#		duplicated.append(row)

IPython.embed()

                        
