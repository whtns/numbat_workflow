#!/usr/bin/python -i

# script generates bed file that scpecifies regions that overlap coverage above defined level
# usage: samtools depth <file.bam> | bed_from_areas_covered_above_N_v2.py N

import sys
import re
import pandas as pd
import subprocess
import fileinput
#import pdb

if(len(sys.argv) < 2):
	print "usage: samtools depth <file.bam> | bed_from_areas_covered_above_N.py N"
	exit(0)

#pdb.set_trace()
N = int(sys.argv[1])

# read only first line
line = sys.stdin.readline()
x = line.rstrip().split("\t")
chrom=x[0]
start_pos = int(x[1])
last_pos = start_pos
num = int(x[2])
if(num > N):
	isPileup = True
else:
	isPileup = False

while (1):
	line = sys.stdin.readline()
	if not line:
		break
	x = line.rstrip().split("\t")
	this_pos = int(x[1])
	if(last_pos != this_pos-1):
		if(isPileup):
			print chrom+"\t"+str(start_pos)+"\t"+str(last_pos)
		chrom = x[0]
		start_pos = int(x[1])
		last_pos = int(x[1])
		isPileup = False
	
	if(int(x[2]) > N): # depth is more than we want
		isPileup = True
	
	last_pos = this_pos

