#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
 	argv[1] = directory with fastq files
	argv[2] = organism
	argv[3] = reference genome (hg19 or hg38)
	argv[4] = output dir
	argv[5] = step to process'
    exit 0
fi

directory_with_fastq_files=$1 #"list of input fastq files"
#~ sample_sheet=$2
organism=$2
reference_genome=$3
outdir=$4
step_to_process=$5

for i in `find $directory_with_fastq_files -name "*R1*fastq.gz" | sort -V | grep -v "Undetermined"`;do 
	j=`echo $i | sed 's/R1/R2/'`
	cell_name=`echo $i | sed 's#.*/\(.*\)R1.*#\1#'`
	echo $i $j $cell_name
	./src/single_cell_pipeline.py -1 $i -2 $j -a AAGCAGTGGTATCAA -o \
	$organism -rg $reference_genome -d $outdir  -s $step_to_process
done

#~ # compile stringtie output for all samples into 1) tpm gene x cell matrix 2) fpkm gene x cell matrix with name prefix as argv 2
#~ find $outdir -name "*.gtf" | sort -V | xargs ./src/gtfs_to_TPM-FPKM.py transcripts

# convert transcript counts to census counts, save as <input_prefix>_census_matrix.csv
#~ ./src/convert_to_census.R $outdir/transcripts.tpm.csv $sample_sheet
