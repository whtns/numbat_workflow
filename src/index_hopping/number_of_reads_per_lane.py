#!/usr/bin/python
from __future__ import print_function
import sys
import pandas as pd
import gzip

# argv[1:] = R1.fastq.gz filenames. Each file needs to start with cell number followed by underscore

# print to stderr
def errprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

reads_per_lane = [0]*8

for fastq in sys.argv[1:]:
	errprint("Processing file:",fastq)
	if fastq.endswith(".gz"):
		f = gzip.open(fastq)
	else:
		f = open(fastq)
	
	for line in f:
		if line.startswith("@"):
			reads_per_lane[int(line.split(":")[3])-1] += 1
	
	f.close()

for i,l in enumerate(reads_per_lane):
	print("Lane ",l," has ",2*l," reads")
