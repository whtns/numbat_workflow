#!/bin/bash

# argv[1] = directory with the input fastq files
# argv[2] = organism
# argv[3] = output dir
# argv[4] = step to process

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
 	argv[1] = list of input fastq files
	argv[2] = organism
	argv[3] = reference genome (hg19 or hg38)
	argv[4] = output dir
	argv[5] = step to process'
    exit 0
fi

list_to_process=$1 #"list of input fastq files"
organism=$2
reference_genome=$3
outdir=$4
step_to_process=$5


for b in `cat $list_to_process | sort -V | grep -v "Undetermined"`;do 
	j=`echo $b | sed 's/R1/R2/'`
	cell_name=`echo $b | sed 's#.*/\(.*\)R1.*#\1#'`
	echo $b $j $cell_name
	./src/single_cell_pipeline.py -1 $b -2 $j -a AAGCAGTGGTATCAA -o \
	$organism -rg $reference_genome -d $outdir  -s $step_to_process
done

