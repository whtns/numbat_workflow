#!/usr/bin/python
# args = fastqc zip files

import sys
from zipfile import *

z = ZipFile(sys.argv[1])
nl = z.namelist()
f = [fn for fn in nl if "summary.txt" in fn]

if(len(f) != 1):
	print "Error: did not find summary.txt"
	exit(1)
	
su = z.read(f[0]).splitlines()
z.close()

if(len(sys.argv) > 1):
	header = "File"
	for i in su:
		header += "\t" + i.split("\t")[1]
	print header



for zipfile in sys.argv[1:]:
	try:
		z = ZipFile(zipfile)
	except:
		print "Problem with:",zipfile
		raise
	nl = z.namelist()
	f = [fn for fn in nl if "summary.txt" in fn]
	
	if(len(f) != 1):
		print "Error: did not find summary.txt"
		exit(1)
		
	su = z.read(f[0]).splitlines()
	z.close()
	out = su[0].split("\t")[2]
	for i in su:
		out += "\t" + i.split("\t")[0]
	print out
	
