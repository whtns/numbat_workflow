#!/usr/bin/python

import sys
import re
import pandas as pd

header = ["IN:reads", "IN:paired", "%", "PAIR-AL-C:0", "%", "PAIR-AL-C:1", "%", "PAIR-AL-C:>1", "%", "PAIR-AL-D:1", "%"]
header +=["MATE-AL:0", "%", "MATE-AL:1", "%", "MATE-AL:>1", "%"]
header +=[ "IN:unpaired", "%", "AL:0", "%", "AL:1", "%", "AL:>1", "%", "overall alignment rate"]
data = pd.DataFrame(columns=header)

for i in sys.argv[1:]:
	cell=i.split("/")[-1].split(".")[0]
	row=[]
	with open(i) as f:
		l = f.read().splitlines()
	m = re.search("([0-9]*) reads; of these:",l[0])
	row += [m.group(1)]
	m = re.search(" ([0-9]*) \((.*)%\) were paired; of these:",l[1])
	row += [m.group(1), m.group(2)]
	m = re.search(" ([0-9]*) \((.*)%\) aligned concordantly 0 times",l[2])
	row += [m.group(1), m.group(2)]
	m = re.search(" ([0-9]*) \((.*)%\) aligned concordantly exactly 1 time",l[3])
	row += [m.group(1), m.group(2)]
	m = re.search(" ([0-9]*) \((.*)%\) aligned concordantly >1 times",l[4])
	row += [m.group(1), m.group(2)]
	m = re.search(" ([0-9]*) \((.*)%\) aligned discordantly 1 time",l[7])
	row += [m.group(1), m.group(2)]
	
	m = re.search(" ([0-9]*) \((.*)%\) aligned 0 times",l[11])
	row += [m.group(1), m.group(2)]
	m = re.search(" ([0-9]*) \((.*)%\) aligned exactly 1 time",l[12])
	row += [m.group(1), m.group(2)]
	m = re.search(" ([0-9]*) \((.*)%\) aligned >1 times",l[13])
	row += [m.group(1), m.group(2)]
	
	
	m = re.search(" ([0-9]*) \((.*)%\) were unpaired; of these:",l[14])
	row += [m.group(1), m.group(2)]
	m = re.search(" ([0-9]*) \((.*)%\) aligned 0 times",l[15])
	row += [m.group(1), m.group(2)]
	m = re.search(" ([0-9]*) \((.*)%\) aligned exactly 1 time",l[16])
	row += [m.group(1), m.group(2)]
	m = re.search(" ([0-9]*) \((.*)%\) aligned >1 times",l[17])
	row += [m.group(1), m.group(2)]
	m = re.search("(.*)% overall alignment rate",l[18])
	row += [m.group(1)]
	#row = map(float,row)
	row = pd.Series(row, index=header, name=cell)
	data = data.append(row)
	print cell
	
data.to_csv("hisat2_tweaked_alignment_stats.csv", sep="\t")
