#!/usr/bin/python

# argv[1] = cell_to_index_and_lane.csv
#cell	lanes	i7	i5
#1	1,2,3,4	TAAGGCGA	AGAGGATA
#2	1,2,3,4	TAAGGCGA	CTCCTTAC
#3	1,2,3,4	TAAGGCGA	TATGCAGT

# argv[2] = format of file with read counts (output of index_hopping_rate.py). Replace lane number with the variable "%d"
# example: "index_hopping_read_counts_lane-%d.csv"

# argv[3] = output format, subtitute lane number for "%d"
# example: "index_hopping_read_counts_lane-%d_filled.csv"

# argv[4] = number of reads in each cell
# Example:
# Cell	Total_Sequences
# 1	55354
# 2	8194066
# 3	8191208

from __future__ import print_function
import sys
import pandas as pd
import gzip
import itertools
import ipdb
from IPython import embed

table_of_valid_indeces_file = sys.argv[1]
read_count_table_filename_format = sys.argv[2]
output_format = sys.argv[3]
cell_read_cout_file = sys.argv[4]

# print to stderr
def errprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

if(len(sys.argv)< 3):
	errprint("Usage: ./check_for_accepted_errors_in_indices.py FILE.fastq[.gz]")
	errprint(sys.argv)
	exit(1)

# read tables with hopped indices
tables = {}
def read_table(lane):
	tables[lane] = pd.read_csv(read_count_table_filename_format % lane, sep="\t", index_col=0)


valid_i7_on_lane = {}
valid_i5_on_lane = {}
valid_pairs_on_lane = {}

lanes = []
with open(table_of_valid_indeces_file) as t:
	t.readline() # skip header
	for line in t:
		x=line.rstrip().split("\t")
		lanes_here = map(int,x[1].split(","))
		for l in lanes_here:
			if l not in lanes:
				valid_i7_on_lane[l] = set()
				valid_i5_on_lane[l] = set()
				valid_pairs_on_lane[l] = set()
				lanes.append(l)
			valid_i7_on_lane[l].add(x[2])
			valid_i5_on_lane[l].add(x[3])
			if((x[2],x[3]) in valid_pairs_on_lane[l]):
				errprint("Index pair used twice on the same lane! ",line)
				exit(1)
			valid_pairs_on_lane[l].add((x[2],x[3]))

lanes = valid_pairs_on_lane.keys()
invalid_pairs_on_lane = {}

for lane in lanes:
	invalid_pairs_on_lane[lane] = [i for i in itertools.product(valid_i7_on_lane[lane], valid_i5_on_lane[lane]) if i not in valid_pairs_on_lane[lane] ]

for lane in lanes:
	read_table(lane)
	valid_table = tables[lane].copy()
	for i in invalid_pairs_on_lane[lane]:
		valid_table.loc[i[0],i[1]] = 0
		#print(i)
	hopping_ratio = pd.DataFrame(0, index=valid_table.index, columns=valid_table.columns)
	for i in invalid_pairs_on_lane[lane]:
		hopping_ratio.loc[i[0],i[1]] = tables[lane].loc[i[0],i[1]] / float(valid_table.loc[i[0],].sum()+valid_table[i[1]].sum())

	for i in valid_pairs_on_lane[lane]:
		rowcol_sum = float(hopping_ratio.loc[i[0],].sum()+hopping_ratio[i[1]].sum())
		hopping_ratio.loc[i[0],i[1]] = rowcol_sum/len(hopping_ratio.loc[i[0],]+hopping_ratio[i[1]])
	print(valid_table)
	print(hopping_ratio)
	hopping_ratio.to_csv(output_format % lane, sep="\t")

cell_read_cout = pd.read_csv(cell_read_cout_file, sep="\t", index_col=0).iloc[:,0]
indices = pd.read_csv(table_of_valid_indeces_file, sep = "\t")
tuples = [tuple(x) for x in indices[['i7', 'i5']].values]

ipdb.set_trace()

with open(cell_read_cout_file) as t:
	t.readline() # skip header 
	for line in t:    
		x=line.rstrip().split("\t") 
		cell_id = int(x[0])
		test = indices[(indices['cell'] == cell_id)]
		for i in tuples:
			p = hopping_ratio.loc[i]
			indices.loc[(indices['i7'] == tuples[i][0]) & (indices['i5'] == tuples[i][1])][0]
			invalid_reads = p * cell_read_cout[i]
total_reads = reads_mat.values.sum() 




'''
with open(table_of_valid_indeces_file) as t:
	t.readline() # skip header
	for line in t:
		x=line.rstrip().split("\t")
		print(x)
		cell = int(x[0])
		lanes = map(int, x[1].split(","))
		i7 = x[2]
		i5 = x[3]
		reads_per_lane = cell_read_cout[cell] / len(lanes)
		for l in lanes:
			if l not in tables:
				read_table(l)
			tables[l].loc[i7,i5] = reads_per_lane

for lane in tables:
	tables[lane].to_csv(output_format % lane, sep="\t")
'''
