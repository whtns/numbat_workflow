#!/bin/bash 

indir=$1
outdir=$2

for i in `find $indir -name "*removed_duplicates.bam" | sort -V`
do  
	echo $i
	j=$outdir/`basename $i | sed 's/.bam/.tdf/'`
	echo $j
	/usr/local/bin/TOOLS/IGVTools/igvtools count $i $j hg38
done

