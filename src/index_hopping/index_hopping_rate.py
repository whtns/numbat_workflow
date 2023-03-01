#!/usr/bin/python


# Copyright (c) 2017 Martin Triska
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# Following software helps you to estimate index hopping rate from single cell
# Illumina sequencing run (using i7 and i5 indices to demultiplex reads).
# For details, please see publication <citation>

# =============================================================================
# usage:
# python index_hopping_rate.py cell_to_index_and_lane.csv per_cell_read_counts.csv Undetermined_reads.fastq.gz [Undetermined_reads.fastq.gz]*

# argv[1] = cell_to_index_and_lane.csv
# Example:
# cell	lanes	i7	i5
# 1	1,2,3,4	TAAGGCGA	AGAGGATA
# 2	1,2,3,4	TAAGGCGA	CTCCTTAC
# 3	1,2,3,4	TAAGGCGA	TATGCAGT

# argv[2] = number of reads in each cell
# Example:
# Cell	Total_Sequences
# 1	55354
# 2	8194066
# 3	8191208

# argv[3:] = Undetermined.fastq[.gz] file(s)


from __future__ import print_function
import sys
import pandas as pd
import gzip
 

## Function extracts lane and indices from read header in fastq file.
# It is supposed to return list with 3 elements: [lane, i7_index, i5_index]
# The fastq header lines this function was written for look as follows:
# @K00187:32:HKCGYBBXX:7:1101:5741:1103 1:N:0:NAAGGCGA+NGAGGATA
# If read headers in your fastq look different, please change this line
# to extract the information properly. If you don't know how, you can always 
# email me at martin.triska@gmail.com with sample of your fastq file, and
# I will write a replacement function for you (it takes few seconds)
def extract_lane_and_indices(line):
	return [int(line.split(":")[3])] + line.rstrip().split(":")[-1].split("+")

## print to stderr

def errprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

if(len(sys.argv)< 3):

	errprint("Usage: python index_hopping_rate.py cell_to_index_and_lane.csv per_cell_read_counts.csv Undetermined_reads.fastq.gz [Undetermined_reads.fastq.gz]*")

	errprint(sys.argv)
	exit(1)

table_of_valid_indeces_file = sys.argv[1]

cell_read_cout_file = sys.argv[2]
fastqs = sys.argv[3:]

err_allowed = 1

## return distance of two strings of equal length
def str_distance(s1,s2):
	return sum(1 for a, b in zip(s1, s2) if a != b)

## return which index does this index match (possibly with error(s))

def matches_index(s,t,err=1):
	if s in t:
		return s
	for i in t:
		if(str_distance(s,i) <= err):
			return i
	return False


## read table of used indices on each line

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

# print how many indices were used on each line
for l in valid_pairs_on_lane:
	print("Lane: ",l)
	print("# valid i7:",len(valid_i7_on_lane[l]))
	print("# valid i5:",len(valid_i5_on_lane[l]))
	print("# valid pairs:",len(valid_pairs_on_lane[l]))
	print("="*100)


lane_to_readcount_table = {}
for l in lanes:
	lane_to_readcount_table[l] = pd.DataFrame(0,index=valid_i7_on_lane[l], columns=valid_i5_on_lane[l])

# create an empty table of read counts for each lane
cell_read_cout = pd.read_csv(cell_read_cout_file, sep="\t", index_col=0).iloc[:,0]

# pre-fill the table with read counts of cells
with open(table_of_valid_indeces_file) as t:
	t.readline() # skip header
	for line in t:
		x=line.rstrip().split("\t")
		cell = int(x[0])
		lanes_here = map(int, x[1].split(","))
		i7 = x[2]
		i5 = x[3]
		reads_per_lane = cell_read_cout[cell] / len(lanes_here)
		for l in lanes_here:
			lane_to_readcount_table[l].loc[i7,i5] = reads_per_lane


# now we process the "Undetermined.fastq[.gz]" reads and find out how many reads have index combinations,
# that should not exist in the dataset

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

			if lane not in lanes: # this read does not belong to our lanes
				continue
			m7 = matches_index(i7,valid_i7_on_lane[lane], err_allowed)
			m5 = matches_index(i5,valid_i5_on_lane[lane], err_allowed)
			if(m7 and m5): # if both indices matched valid index on the lane
				lane_to_readcount_table[lane].loc[m7,m5] += 1 # add read to the table

	f.close()
	files_processed += 1
	print("Momentary read counts after ", files_processed," files:")
	for l in lanes:
		print("-"*70)
		print("LANE ",l)

		print(lane_to_readcount_table[l])

 



