#!/usr/bin/python

import sys
import subprocess

for f in sys.argv[1:]:
	#print f
	basename = f.split("/")[-1]
	cell = basename.split("_")[0]
	
	grep = subprocess.Popen(["zgrep","^@",f], stdout=subprocess.PIPE)
	cut  = subprocess.Popen(["cut", "-f", "4", "-d", ":"], stdin=grep.stdout, stdout=subprocess.PIPE)
	uniq1= subprocess.Popen(["uniq"], stdin=cut.stdout, stdout=subprocess.PIPE )
	sort = subprocess.Popen(["sort"], stdin=uniq1.stdout, stdout=subprocess.PIPE )
	uniq2= subprocess.Popen(["uniq"], stdin=sort.stdout, stdout=subprocess.PIPE )
	
#	grep.stdout.close()
#	cut.stdout.close()
#	uniq1.stdout.close()
#	sort.stdout.close()
#	uniq2.stdout.close()
	
	out = uniq2.communicate()[0]
	out = out.rstrip().replace("\n",",")
	print cell+"\t"+out
