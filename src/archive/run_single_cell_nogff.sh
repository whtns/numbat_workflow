#!/bin/bash

# argv[1] = directory with the input fastq files
# argv[2] = organism

directory_with_fastq_files=$1 #"input_fastq"
organism=$2
for i in `find $directory_with_fastq_files -name "*R1_001.SMARTer_trimmed.trimmed.fastq"`;do
	j=`echo $i | sed 's/R1/R2/'`
	file=`echo $i | sed 's/.*\///'`
	cell_name=`echo $file | sed 's;_.*;;'`
	echo $i $j $cell_name $file
	./single_cell_pipeline_nogff.py $i $j $organism $cell_name
done
