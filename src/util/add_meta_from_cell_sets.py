#!/usr/bin/python

import pandas as pd
import sys
import IPython

new_meta = sys.argv[1]	# tab delimited cell_info shCtrl	X1	X2 
old_meta = sys.argv[2]	#sample sheet with cell metadata
 
with open(new_meta) as nm:          
    content = nm.readlines()
    content = [x.rstrip().split("\t") for x in content]
    new_dict = {}
    for i in content:
		for m in i[1:]:
			new_dict[m] = i[0]

test = pd.read_csv(old_meta, sep = ",")
     
test['branch'] = test['sample_id'].map(new_dict)
print old_meta
# ~ IPython.embed()

test.to_csv(old_meta, sep=",", index=False)
    
