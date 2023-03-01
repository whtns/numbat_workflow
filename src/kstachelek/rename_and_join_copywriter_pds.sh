#!/bin/bash 

copywriter_dir=$1

for i in `find *none -name "all_chrom.pdf"`;do
	j=`echo $i | sed 's/_[A-Z].*//'`
	rsync -a $i all_chrom_pdfs/$j
done
