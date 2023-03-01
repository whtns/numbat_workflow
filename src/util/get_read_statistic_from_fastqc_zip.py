#!/usr/bin/python
# args = fastqc zip files

import sys
import re
from zipfile import *

print "cell\tread_count"

for zipfile in sys.argv[1:]:
	try:
		z = ZipFile(zipfile)
	except:
		print "Problem with:",zipfile
		raise
	nl = z.namelist()
	f = [fn for fn in nl if "fastqc_data.txt" in fn]
	
	if(len(f) != 1):
		print "Error: did not find fastqc_data.txt"
		exit(1)
		
	su = z.read(f[0]).splitlines()
	z.close()
	
	name = [n for n in su if "Filename\t" in n]
	if(len(name) != 1):
		print "Error: did not find \"Filename\" in the file"
		exit(1)
	name = name[0].split("\t")[1].split("_")[0]

#	name = name[0]+"_"+name[1]
	
	count = [c for c in su if "Total Sequences\t" in c]
	if(len(count) != 1):
		print "Error: did not find \"Total Sequences\" in the file"
		exit(1)
	count = count[0].split("\t")[1]
	count = str(int(count)*2)
	print name+"\t"+count
	
