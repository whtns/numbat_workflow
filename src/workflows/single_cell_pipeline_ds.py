#!/usr/bin/python

# argv[1,2] = R1/R2.fastq files
# argv[3] = organism [human/mouse]

import subprocess
import sys
import re
import os
import datetime

steps_to_process = ["all"] # 
num_threads = "4"
# *********************************************************************
# DEFINITION OF PATHS

gtf_juncs        = "gtf_juncs"
tophat_binary    = "tophat"
prep_reads_binary= "prep_reads"
cutadapt_binary  = "cutadapt"
cufflinks_binary = "cufflinks"
cuffmerge_binary = "cuffmerge"
cuffquant_binary = "cuffquant"
cuffnorm_binary = "cuffnorm"
cuffdiff_binary = "cuffdiff"
trimmomatic_cmd  = ["java", "-jar", "/usr/local/bin/trimmomatic-0.36.jar"]
min_base_quality = "15"
min_read_length  = "25"
sliding_window   = "6:20" # window_size:min_average_quality
adapter_sequence = "AAGCAGTGGTATCAA"

fastq_r1_location = sys.argv[1]
fastq_r2_location = sys.argv[2]
print sys.argv[1]
print sys.argv[2]  
print sys.argv[3]

if(len(sys.argv)>3):
	if(sys.argv[3]=="human"):
		gtf_location = "/dataVolume/storage/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
		bowtie2index_location = "/dataVolume/storage/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
		output_directory_root = "./human/" 
		reference_genome = "/dataVolume/storage/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
	elif(sys.argv[3]=="mouse"):
		gtf_location = "/dataVolume/storage/Mus_musculus/NCBI/GRCm38/Annotation/Genes/genes.gtf" # /dataVolume/storage/Mus_musculus/UCSC/GRCm38/Annotation/Genes/genes.gtf
		bowtie2index_location = "/dataVolume/storage/Mus_musculus/NCBI/GRCm38/Sequence/Bowtie2Index/genome" # /dataVolume/storage/Mus_musculus/UCSC/GRCm38/Sequence/Bowtie2Index/genome
		output_directory_root = "./mouse/" 
	else:
		print "organism not recognised! Please use human or mouse"
		exit()
	transcriptome_index_location = gtf_location.rsplit(".",1)[0]
else:
	print "Please specify species"
	exit()
 
print bowtie2index_location
print gtf_location
r1_base_filename = fastq_r1_location.replace(".fastq.gz","").split("/")[-1]
r2_base_filename = fastq_r2_location.replace(".fastq.gz","").split("/")[-1]
fastq_r1 = fastq_r1_location.split("/")[-1]
fastq_r2 = fastq_r2_location.split("/")[-1]
cell_name = r1_base_filename.split(".R")[0]
print cell_name
output_directory = "/dataVolume/storage/single_cell_pipeline/dshayler/"
cell_output_dir = output_directory+cell_name+"/"
subprocess.call(["mkdir",cell_output_dir])

# /DEFINITION OF PATHS
# *********************************************************************

# CUTADAPT
# =====================================================================
#cutadapt -b AAGCAGTGGTATCAA -B AAGCAGTGGTATCAA -m 20 --overlap 5 -o DC-0001_S1_R1_001.trimmed.fastq --paired-output DC-0001_S1_R2_001.trimmed.fastq DC-0001_S1_R1_001.fastq.gz DC-0001_S1_R2_001.fastq.gz > DC-0001.trim.log
cutadapt_output_dir = cell_output_dir +"cutadapt/"
fastq_r1_cutadapt = cutadapt_output_dir+ r1_base_filename + ".cutadapt.fastq"
fastq_r2_cutadapt = cutadapt_output_dir+ r2_base_filename + ".cutadapt.fastq"
if("cutadapt" in steps_to_process) or ("all" in steps_to_process):
	subprocess.call(["mkdir",cutadapt_output_dir])
	out_file = cell_output_dir +"cutadapt/"+cell_name+"_cutadapt.log"
	print "fastq_r1_cutadapt:",fastq_r1_cutadapt
	print "fastq_r2_cutadapt:",fastq_r2_cutadapt
	print "out_file: ",out_file
	cutadapt_arguments = ["-b", adapter_sequence, "-B", adapter_sequence, "-m", "20", "--overlap", "5"]
	cmd = [cutadapt_binary]
	cmd += cutadapt_arguments
	cmd += ["-o",fastq_r1_cutadapt,"--paired-output",fastq_r2_cutadapt]
	cmd += [fastq_r1_location, fastq_r2_location]
	
	print " ".join(cmd)
	with open(out_file,"w")as outf:
	 	subprocess.call(cmd,stdout=outf)
	del cmd

# TRIMMOMATIC
# =====================================================================
trimmomatic_output_dir = cell_output_dir + "trimmomatic/"
if("trimmomatic" in steps_to_process) or ("all" in steps_to_process):
	subprocess.call(["mkdir",trimmomatic_output_dir])
	print "running trimmomatic"
	fastq_r1_trimmed = trimmomatic_output_dir+ r1_base_filename + ".trimmed.fastq"
	fastq_r2_trimmed = trimmomatic_output_dir+ r2_base_filename + ".trimmed.fastq"
	fastq_r1_trimmed_unpaired = trimmomatic_output_dir+ r1_base_filename + ".trimmed.unpaired.fastq"
	fastq_r2_trimmed_unpaired = trimmomatic_output_dir+ r2_base_filename + ".trimmed.unpaired.fastq"
	cmd = trimmomatic_cmd[:]
	#print cmd
	cmd.append("PE")
	cmd.append("-threads")
	cmd.append("4")	
	#cmd.append("-trimlog")
	#cmd.append(log)
	cmd.append(fastq_r1_location)
	cmd.append(fastq_r2_location)
	cmd.append(fastq_r1_trimmed)
	cmd.append(fastq_r1_trimmed_unpaired)
	cmd.append(fastq_r2_trimmed)
	cmd.append(fastq_r2_trimmed_unpaired)
	cmd.append("LEADING:"+min_base_quality)
	cmd.append("TRAILING:"+min_base_quality)
	cmd.append("SLIDINGWINDOW:"+sliding_window)
	cmd.append("MINLEN:"+min_read_length)
	print " ".join(cmd)
	subprocess.call(cmd)
	del cmd

# FASTQC
# =====================================================================
if("fastqc" in steps_to_process) or ("all" in steps_to_process):
	print "running fastqc on the trimmed file"
	cmd = ["fastqc",fastq_r1_trimmed,fastq_r2_trimmed,fastq_r1_location,fastq_r2_location]
	print " ".join(cmd)
	subprocess.call(cmd)
	del cmd

# TOPHAT
# =====================================================================
#/usr/local/bin/tophat --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 --max-multihits 10 --library-type fr-unstranded --GTF /Lab_Share/iGenomes/hg19/Annotation/Genes/genes.gtf --transcriptome-index /Lab_Share/iGenomes/hg19/Annotation/Genes/genes --num-threads 4 --output-dir /Lab_Share/David_Cobrinik/011816/human/tophat/DC-0084 /Lab_Share/iGenomes/hg19/Sequence/Bowtie2Index/genome /Lab_Share/David_Cobrinik/011816/fastq/DC-0084_S84_R1_001.SMARTer_trimmed.trimmed.fastq /Lab_Share/David_Cobrinik/011816/fastq/DC-0084_S84_R2_001.SMARTer_trimmed.trimmed.fastq
fastq_r1_trimmed = trimmomatic_output_dir+ r1_base_filename + ".trimmed.fastq"
fastq_r2_trimmed = trimmomatic_output_dir+ r2_base_filename + ".trimmed.fastq"
fastq_r1_trimmed_unpaired = trimmomatic_output_dir+ r1_base_filename + ".trimmed.unpaired.fastq"
fastq_r2_trimmed_unpaired = trimmomatic_output_dir+ r2_base_filename + ".trimmed.unpaired.fastq"
tophat_output_dir = cell_output_dir + "tophat/"
tophat_output_bam = tophat_output_dir + "accepted_hits.bam"
if("tophat" in steps_to_process) or ("all" in steps_to_process):
	print "processing tophat"
	out_file = tophat_output_dir +cell_name+"_tophat.log"
	subprocess.call(["mkdir",tophat_output_dir])
	tophat_arguments = ["--read-mismatches", "2", "--read-gap-length", "2", "--read-edit-dist", "2", "--max-multihits", "10", "--library-type", "fr-unstranded"]
	print bowtie2index_location
	cmd = [tophat_binary]
	cmd += tophat_arguments
	cmd += ["--GTF",gtf_location]
	cmd += ["--transcriptome-index",transcriptome_index_location]
	cmd += ["--num-threads", num_threads]
	cmd += ["--output-dir",tophat_output_dir]
	cmd += [bowtie2index_location]
	cmd += [fastq_r1_trimmed,fastq_r2_trimmed]
	print " ".join(cmd)
	with open(out_file,"w")as outf:
	 	subprocess.call(cmd,stdout=outf)
	del cmd

# CUFFLINKS
# ====================================================================
cufflinks_output_dir = cell_output_dir + "cufflinks/"
cufflinks_gtf = cufflinks_output_dir+"transcripts.gtf"
if("cufflinks" in steps_to_process) or ("all" in steps_to_process):
	if not os.path.exists(cufflinks_output_dir):
		os.makedirs(cufflinks_output_dir)
	print "processing cufflinks"
	cmd = [cufflinks_binary]
	cmd += ["-p",num_threads]
	cmd += ["--GTF",gtf_location]
	cmd += ["-o",cufflinks_output_dir]
	cmd += [tophat_output_dir+"accepted_hits.bam"]
	out_file = cufflinks_output_dir + "cufflinks.log"
	err_file = cufflinks_output_dir + "cufflinks.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr, open(out_file,"w")as outf:
		 subprocess.call(cmd,stdout=outf,stderr=outerr)
	del cmd
	


# CUFFMERGE
# =====================================================================
cuffmerge_output_dir = output_directory + "cuffmerge/"
now = datetime.datetime.today()
tuxedo_time_dir = now.strftime('%H_%M_%d_%m_%Y/')
cuffmerge_assembly_list = output_directory+"cuffmerge_assembly_list.txt"
cuffmerge_gtf = cuffmerge_output_dir+"merged.gtf"
if("cuffmerge" in steps_to_process) or ("all" in steps_to_process) and (not os.path.isfile(cuffmerge_gtf)):
	if not os.path.exists(cuffmerge_output_dir):
		os.makedirs(cuffmerge_output_dir)
	with open(cuffmerge_assembly_list, 'w') as f:
		for root, dirs, files in os.walk("/dataVolume/storage/single_cell_pipeline/dominic"):
			for file in files:
				if file.endswith("transcripts.gtf"):
					f.write(file+'\n')
		f.truncate()
"""
	print "processing cuffmerge"
	cmd = [cuffmerge_binary]
	cmd += ["-p",num_threads]
	cmd += ["-g",gtf_location]
	cmd += ["-o",cuffmerge_output_dir]
	cmd += [cuffmerge_assembly_list]
	out_file = cufflinks_output_dir + "cufflinks.log"
	err_file = cufflinks_output_dir + "cufflinks.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	subprocess.call(cmd,stdout=outf,stderr=outerr)
	del cmd
"""
# CUFFQUANT
# =====================================================================
cuffquant_output_dir = cell_output_dir + "cuffquant/"
cuffquant_cxb = cuffquant_output_dir+"abundances.cxb"
cuffquant_processed_log = output_directory+"processed_cxb.txt" 
if("cuffquant" in steps_to_process) or ("all" in steps_to_process):
	if not os.path.exists(cuffquant_output_dir):
		os.makedirs(cuffquant_output_dir)
	print "processing cuffquant"
	cmd = [cuffquant_binary]
	cmd += ["-p",num_threads]
	cmd += ["-o",cuffquant_output_dir]
	cmd += ["-b",reference_genome]
	cmd += ["-u",cuffmerge_gtf]
	cmd += [tophat_output_bam]
	out_file = cuffquant_output_dir + "cuffquant.log"
	err_file = cufflinks_output_dir + "cuffquant.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	subprocess.call(cmd,stdout=outf,stderr=outerr)
		 	"""
	with open(cuffquant_processed_log, 'w') as f:
		f.write(cuffquant_cxb+"\n")
		"""
	del cmd

# CUFFNORM
# =====================================================================
cuffnorm_output_dir = output_directory + "cuffnorm/"
sample_sheet = output_directory + "SAMPLE_SHEET.ds_test.human.txt"
if("cuffnorm" in steps_to_process) or ("all" in steps_to_process):
	if not os.path.exists(cuffnorm_output_dir):
		os.makedirs(cuffnorm_output_dir)
	with open(cuffquant_processed_log, 'r') as f:
		cxb_list = [line.rstrip('\n') for line in f]
	print "processing cuffnorm"
	cmd = [cuffnorm_binary]
	cmd += ["-p",num_threads]
	cmd += ["-o",cuffnorm_output_dir]
	cmd += ["-use-sample-sheet"]
	cmd += [cuffmerge_gtf]
	cmd += [sample_sheet]
	for i in cxb_list:
		cmd += [i]
	out_file = cuffnorm_output_dir + "cuffnorm.log"
	err_file = cuffnorm_output_dir + "cuffnorm.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	subprocess.call(cmd,stdout=outf,stderr=outerr)
	del cmd

# CUFFDIFF
# =====================================================================
cuffdiff_output_dir = cell_output_dir + "cuffdiff/"
if("cuffdiff" in steps_to_process) or ("all" in steps_to_process):
	if not os.path.exists(cuffdiff_output_dir):
		os.makedirs(cuffdiff_output_dir)
	with open(cuffquant_processed_log, 'r') as f:
		cxb_list = [line.rstrip('\n') for line in f]
	print "processing cuffdiff"
	cmd = [cuffdiff_binary]
	cmd += ["-p",num_threads]
	cmd += ["-o",cuffdiff_output_dir]
	cmd += [gtf_location]
	for i in cxb_list:
		cmd += [i]
	out_file = cuffdiff_output_dir + "cuffdiff.log"
	err_file = cuffdiff_output_dir + "cuffdiff.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	subprocess.call(cmd,stdout=outf,stderr=outerr)
	del cmd

