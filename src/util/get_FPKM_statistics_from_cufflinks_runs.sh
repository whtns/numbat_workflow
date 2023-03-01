#!/bin/bash
# $1 = file name like "genes.fpkm_tracking"
# $2 = which column should be extracted
# $3 = output file
# $4 = input folder
echo "" > $3
for i in `find $4 -regex ".*/cufflinks/$1" | sort -t "/" -k2g | head -n 1`;do # | sort -t "_" -k2g |
	tail -n +2 $i | cut -f 1 | sort >> $3
done
for i in `find $4 -regex ".*/cufflinks/$1" | sort -t "/" -k2g`;do
	echo $i
	cell_name=`echo $i | cut -f 2 -d'/' | cut -f 1,2 -d'_'`
	echo "${cell_name}" > tmp$2
	tail -n +2 $i | sort | cut -f $2 >> tmp$2
	paste $3 tmp$2 > tmp_res
	#read
	mv tmp_res $3
	rm tmp$2
done
