#!/usr/bin/python -i


import sys
import re
import pandas as pd
import IPython

if len(sys.argv) == 1:
    print "provide: 1) output filename 2) multiqc_general_stats 3) hisat2.log files using xargs to run script!"

header = ["IN:reads", "IN:paired", "%", "PAIR-AL-C:0", "%", "PAIR-AL-C:1", "%", "PAIR-AL-C:>1", "%", "PAIR-AL-D:1", "%"]
header +=["MATE-AL:0", "%", "MATE-AL:1", "%", "MATE-AL:>1", "%"]
header +=[ "IN:unpaired", "%", "AL:0", "%", "AL:1", "%", "AL:>1", "%", "overall alignment rate"]
data = pd.DataFrame(columns=header)

outfile = sys.argv[1]
general_stats = sys.argv[2]

for i in sys.argv[3:]:
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
	#print cell

general_stats = pd.read_csv(general_stats, sep = "\t")

def return_cell_id(x):
    return int(re.split("[_-]", x)[0])
    
general_stats["sortby"] = general_stats["Sample"].apply(return_cell_id)
general_stats.sort_values(by="sortby", inplace=True)
general_stats = general_stats[general_stats["Sample"].str.contains("R1")]
general_stats.to_csv("fastqc/general_stats.csv", sep="\t")

#IPython.embed()

#~ print data.index
#~ print general_stats.sortby
#~ print general_stats.shape
#~ for i in range(1,573):
	#~ if i not in general_stats.sortby.values:
		#~ print ">>>", i

data.insert(0, "Vanilla from fastqc", general_stats.loc[:,"Total Sequences"].values)
print data.head()
data.sort_values(by="Vanilla from fastqc", inplace=True)


data.to_csv(outfile, sep="\t")


