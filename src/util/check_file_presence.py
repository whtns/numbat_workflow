#!/usr/bin/python

import sys
import os.path

# argv[] = fastq R1 files to be checked
# FASTQ/_EGAR00001405855_PoolC-I-b16_GCCAAT_L003_R1_001.fastq.gz

header = [""]
header += ["R1.fastq","R2.fastq"]
header += ["cutadapt_R1.fastq", "cutadapt_R2.fastq"]
header += ["R1_trimmed.fastq","R2_trimmed.fastq", "R1_trimmed.unpaired.fastq","R2_trimmed.unpaired.fastq"]
header += ["fastqc_R1", "fastqc_R2"]
header += [".bam"]
header += ["_marked_duplicates.bam","duplicate_matrics"]
header += ["cufflinks_genes", "cufflinks_isoforms"]


print "\t".join(header)
output_directory_root = "output/"

for f1 in sys.argv[1:]:
	fastq_r1_location = f1
	fastq_r2_location = f1.replace("_R1_","_R2_")
	expected_files = [fastq_r1_location,fastq_r2_location]
	
	r1_base_filename = fastq_r1_location.replace(".fastq.gz","").split("/")[-1]
	r2_base_filename = fastq_r2_location.replace(".fastq.gz","").split("/")[-1]
	fastq_r1 = fastq_r1_location.split("/")[-1]
	fastq_r2 = fastq_r2_location.split("/")[-1]
	cell_name = r1_base_filename.split("_")[0]
	cell_output_dir = output_directory_root+"/"+cell_name+"/"
	cutadapt_output_dir = cell_output_dir +"cutadapt/"	
	fastq_r1_cutadapt = cutadapt_output_dir+ r1_base_filename + ".cutadapt.fastq"
	fastq_r2_cutadapt = cutadapt_output_dir+ r2_base_filename + ".cutadapt.fastq"
	expected_files += [fastq_r1_cutadapt, fastq_r2_cutadapt]
	
	trimmomatic_output_dir = cell_output_dir + "trimmomatic/"
	fastq_r1_trimmed = trimmomatic_output_dir+ r1_base_filename + ".trimmed.fastq"
	fastq_r2_trimmed = trimmomatic_output_dir+ r2_base_filename + ".trimmed.fastq"
	fastq_r1_trimmed_unpaired = trimmomatic_output_dir+ r1_base_filename + ".trimmed.unpaired.fastq"
	fastq_r2_trimmed_unpaired = trimmomatic_output_dir+ r2_base_filename + ".trimmed.unpaired.fastq"
	expected_files += [fastq_r1_trimmed, fastq_r2_trimmed, fastq_r1_trimmed_unpaired, fastq_r2_trimmed_unpaired]
	
	fastqc_r1 = fastq_r1_trimmed.replace(".fastq","_fastqc.html")
	fastqc_r2 = fastq_r2_trimmed.replace(".fastq","_fastqc.html")
	expected_files += [fastqc_r1, fastqc_r2]
	
	hisat2_output_dir = cell_output_dir + "hisat2/"
	hisat2_sorted_bam = hisat2_output_dir + cell_name + ".bam"
	expected_files += [hisat2_sorted_bam]
	
	picard_output_dir = cell_output_dir+"picard/"
	remove_duplicates_bam= picard_output_dir+ r1_base_filename + "_removed_duplicates.bam"
	remove_dup_metrics = picard_output_dir+r1_base_filename+"_removed_dup_metrics.txt"
	expected_files += [remove_duplicates_bam, remove_dup_metrics]
	
	cufflinks_output_dir = cell_output_dir + "cufflinks/"
	cufflinks_genes_file = cufflinks_output_dir + "genes.fpkm_tracking"
	cufflinks_isoforms = cufflinks_output_dir + "isoforms.fpkm_tracking"
	expected_files += [cufflinks_genes_file, cufflinks_isoforms]

	#print f1,f2
	
	out = [r1_base_filename]
	for i in expected_files:
		if(os.path.isfile(i)):
			out.append(str(os.path.getsize(i)))
		else:
			out.append("N")
	print "\t".join(out)
