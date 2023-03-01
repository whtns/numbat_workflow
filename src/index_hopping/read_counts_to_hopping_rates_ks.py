#!/usr/bin/python

import sys
import pandas as pd
import ipdb

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

# argv[3:] = lane file(s)

table_of_valid_indeces_file = sys.argv[1]

cell_read_cout_file = sys.argv[2]

lane_file = sys.argv[3]

indices = pd.read_csv(table_of_valid_indeces_file, sep = "\t")
reads_mat = pd.read_csv(lane_file, sep = "\t", index_col=0)

subset = indices[['i7', 'i5']]
tuples = [tuple(x) for x in subset.values]

ipdb.set_trace()

cell_read_cout = pd.read_csv(cell_read_cout_file, sep="\t", index_col=0).iloc[:,0]
    
invalid_reads = int(0)       
for i in tuples:                          
     p = reads_mat.loc[i]
     invalid_reads = p - cell_read_cout[i]
total_reads = reads_mat.values.sum() 

hopping_rate = ((total_reads - valid_reads)/float(total_reads))


misattr_reads = (hopping_rate*cell_read_cout)

hopped_reads = pd.concat([cell_read_cout, misattr_reads.astype(int)], axis=1)
hopped_reads.columns = ["Total_Sequences", "Hopped_Reads"]
print("wait")

#indices.loc[(indices['i7'] == tuples[1][0]) & (indices['i5'] == tuples[1][1])]

		lanes_here = map(int,x[1].split(","))
		for l in lanes_here:
			if l not in test:
				test[l] = set()
			test[l].add(hopping_ratio.loc[i[0],i[1]],)

