#!/usr/bin/python

# argv[1] = cell_to_index_and_lane.csv
#cell	lanes	i7	i5
#1	1,2,3,4	TAAGGCGA	AGAGGATA
#2	1,2,3,4	TAAGGCGA	CTCCTTAC
#3	1,2,3,4	TAAGGCGA	TATGCAGT

# argv[2] = number of reads in each cell
# Cell	Total_Sequences
# 1	27677
# 1	27677

from __future__ import print_function
import sys
import pandas as pd
import gzip

## change following function to correctly extract indices from fastq header lines
# the fastq header lines this function was written for loos as follows:
# @K00187:32:HKCGYBBXX:7:1101:5741:1103 1:N:0:NAAGGCGA+NGAGGATA
# returns list [lane, i7, i5]
def extract_lane_and_indices(line):
	return [line.split(":")[3]]+line.rstrip().split(":")[-1].split("+")

# print to stderr
def errprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

if(len(sys.argv)< 3):
	errprint("Usage: ./check_for_accepted_errors_in_indices.py FILE.fastq[.gz]")
	errprint(sys.argv)
	exit(1)

with open(table_of_valid_indeces_file) as t:
	t.readline() # skip header
	for line in t:
		x=line.rstrip().split("\t")
		lanes = x[1].split(",")
		for l in lanes:
			valid_i7_on_lane[l].add(x[2])
			valid_i5_on_lane[l].add(x[3])
			if((x[2],x[3]) in valid_pairs_on_lane[l]):
				errprint("Index pair used twice on the same lane! ",line)
				exit(1)
			valid_pairs_on_lane[l].add((x[2],x[3]))

lanes = valid_pairs_on_lane.keys()

# print how many indices were used on each line
for l in valid_pairs_on_lane:
	print("Lane: ",l)
	print("# valid i7:",len(valid_i7_on_lane[l]))
	print("# valid i5:",len(valid_i5_on_lane[l]))
	print("# valid pairs:",len(valid_pairs_on_lane[l]))
	print("="*100)


lane_readcount_table = {}
for l in lanes:
	lane_readcount_table[l] = pd.DataFrame(0,index=valid_i7_on_lane[l], columns=valid_i5_on_lane[l])

files_processed = 0
for fastq in fastqs:
	print("processing file:",fastq)
	print("="*100)
	if(fastq.endswith(".gz")):
		f = gzip.open(fastq)
	else:
		f = open(fastq)
	line_num = 0
	for line in f:
		if line.startswith("@"):
			line_num += 1
			if(line_num % 1000000 == 0):
				errprint("Processed ",line_num,"reads")
			lane, i7, i5 = extract_lane_and_indices(line)
			m7 = matches_index(i7,valid_i7_on_lane[lane], err_allowed)
			m5 = matches_index(i5,valid_i5_on_lane[lane], err_allowed)
			if(m7 and m5): # if both indices matched valid index on the lane
				lane_readcount_table[lane].loc[m7,m5] += 1 # add read to the table
	f.close()
	files_processed += 1
	print("Momentary read counts after ", files_processed," files:")
	for l in lanes:
		print("-"*70)
		print("LANE ",l)
		print(lane_readcount_table[l])
		lane_readcount_table[l].to_csv("index_hopping_read_counts_"+str(l)+"_files-processed-"+str(files_processed)+".csv", sep="\t")




