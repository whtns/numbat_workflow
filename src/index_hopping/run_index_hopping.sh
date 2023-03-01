#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
 	argv[1] = cell_to_index_and_lane.csv; stored in ~/single_cell_pipeline/data/project/data/metadata
 	argv[2] = number of reads in each cell; from multiqc report; stored in ~/single_cell_pipeline/data/project/metadata
 	argv[3] = path to undetermined reads; if paired end supply R1 and R2 in sequence
 	'
    exit 0
fi

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

# argv[3:] = Undetermined.fastq[.gz] file(s)

directory_with_fastq_files=$3 #"list of input fastq files"

./index_hopping_rate.py cell_to_index_and_lane.csv per_cell_read_counts.csv Undetermined_reads.fastq.gz [Undetermined_reads.fastq.gz]*


