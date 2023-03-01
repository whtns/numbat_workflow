#!/usr/bin/python

from __future__ import print_function
import sys
import gzip
import os

# this script checks for how many sequencing errors were allowed when
# demultiplexing sequencing data. This information is needed later
# when trying to estimate the amplitude of index hopping
# Usage: ./check_for_accepted_errors_in_indices.py FILE.fastq[.gz]

# argv[1] = .fastq[.gz] file with demultiplexed data

## change following function to correctly extract indices from fastq header lines
# the fastq header lines this function was written for loos as follows:
# @K00187:32:HKCGYBBXX:7:1101:5741:1103 1:N:0:NAAGGCGA+NGAGGATA
def extract_indices(line):
	return line.rstrip().split(":")[-1].split("+")


def errprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

if(len(sys.argv)!= 2):
	errprint("Usage: ./check_for_accepted_errors_in_indices.py FILE.fastq[.gz]")
	errprint(sys.argv)
	exit(1)

filename = sys.argv[1]
if(filename.endswith(".gz")):
	f = gzip.open(filename)
else:
	f = open(filename)

i5_count = {} # index => count
i7_count = {}
for line in f:
	if line.startswith("@"):
		i7, i5 = extract_indices(line)
		if(i5 not in i5_count):
			i5_count[i5] = 0
		if(i7 not in i7_count):
			i7_count[i7] = 0
		i5_count[i5] += 1
		i7_count[i7] += 1

most_common_i5 = i5
most_common_i7 = i7

for i5 in i5_count:
	if i5_count[i5] > i5_count[most_common_i5]:
		most_common_i5 = i5

for i7 in i7_count:
	if i7_count[i7] > i7_count[most_common_i7]:
		most_common_i7 = i7

def str_diff(seq1, seq2):
	return sum(1 for a, b in zip(seq1, seq2) if a != b)

max_errors_i5 = max( [str_diff(most_common_i5, i5) for i5 in i5_count] )
max_errors_i7 = max( [str_diff(most_common_i7, i7) for i7 in i7_count] )

print("Max i5 errors allowed: ",max_errors_i5)
print("Max i7 errors allowed: ",max_errors_i7)
