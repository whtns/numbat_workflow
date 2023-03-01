#!/usr/bin/python

# import system library - takes care of reading arguments from command line
import sys
# pandas library helps with reading csv file and with statistics
import pandas as pd
# library matplotlib takes care of plotting the graphs
import matplotlib.pyplot as plt
# library saborn provides more pleasant color palletes for the graphs
import seaborn as sns

import subprocess
from StringIO import StringIO

import os
import fnmatch

cellset_file = sys.argv[1]
path_to_cells= sys.argv[2]
gene_exone_file = sys.argv[3]
sets = ["shCtrl", "sh733", "sh737"]

# cell sets
cell_sets = {}
with open(cellset_file) as f:
	for line in f:
		x = line.rstrip().split("\t")
		cell_sets[x[0]] = x[1:]

def find_file(path,cell,pattern):
	for root, dirs, files in os.walk(path+cell):
		for f in files:
			if fnmatch.fnmatch(f, pattern):
				return os.path.join(root, f)

set_files = {}
for s in sets:
	set_files[s] = []
	cells = cell_sets[s]
	filenames = []
	for c in cells:
		set_files[s].append(find_file(path_to_cells,c[1:],"*_removed_duplicates.bam"))

gene_exones = open(gene_exone_file).read().splitlines()

min_coverage = 2
coverage = pd.DataFrame(columns=["cell", "exon", "mean_expression", "median_expression", "covered_bases", "treatment"])
for j in range(0,len(sets)):
	s = sets[j]
	cell_names = [x.split("/")[2] for x in set_files[s]]
	for i in range(0,len(gene_exones)):
		e=gene_exones[i]
		cmd = ["samtools","depth","-aa","-r",e]
		cmd += set_files[s]
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
		out, err = p.communicate()
		outf = StringIO(out)
		cov = pd.read_csv(outf, sep="\t", header=None).iloc[:,2:]
		x = pd.DataFrame([
			cell_names,
			[i+0.25*j]*cov.shape[1], #exon coordinate (each set moved 0.25 to display more clearly)
			list(cov.mean()),
			list(cov.median()),
			list(cov[cov>min_coverage].count()),
			[s]*cov.shape[1]
			]).transpose()
		x.columns = ["cell", "exon", "mean_expression", "median_expression", "covered_bases", "treatment"]
		print x
		coverage = coverage.append(x)

# let's be conventional with exon numbers starting from 1 (not 0)
coverage["exon"] = coverage["exon"]+1

coverage["color"]="blue"
coverage.loc[coverage["treatment"]=="sh733","color"]="orange"
coverage.loc[coverage["treatment"]=="sh737","color"]="red"

# annotate point on axis if it's far enough
def annotate_df(df,min_dist,ax, x,y,label_column):
	for index,row in df.iterrows():
		rowx=row[x]
		dist = (df.loc[df[x]==rowx,y] - row[y]).abs().sort_values().iloc[1]
		if(dist > min_dist):
			ax.annotate(row[label_column], [row[x],row[y]],
				xytext=(5,-3), 
				textcoords='offset points',
				size=10, 
				color='darkslategrey')


for plot in ["mean_expression", "median_expression", "covered_bases"]:
	print "plotting:",plot
	ax = coverage.plot.scatter(x="exon", y=plot, c=coverage["color"])
	annotate_df(coverage, 4, ax, "exon", plot, "cell")
	ymin = 0
	ymax = coverage[plot].max()
	xmin = 0.75
	xmax = coverage["exon"].max()
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([ymin, ymax])
	plt.xticks(range(1,28))

plt.show()
