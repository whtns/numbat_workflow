#!/usr/bin/python

import subprocess # runs shell commands within python
import sys
import re
import os
import datetime
import argparse


# list here all steps of the pipeline. Pipeline will run all these steps, if not requested otherwise in -s argument

steps_to_process_all = ["fastqc","cutadapt","trimmomatic", "hisat2","removeduplicates", "uniqify", "threshold_bam","stringtie"]

parser = argparse.ArgumentParser(description="runs single cell pipeline")
parser.add_argument("-se", "--fastq-se", dest="fse", help="fastq SE file REQUIRED", metavar="FILE[.gz]", required=True)
parser.add_argument("-a", "--adapter", dest="adapter", help="comma separated list of adapter sequences to trim REQUIRED", metavar="ADAPTER", required=True)
parser.add_argument("-o", "--organism", dest="organism", default= "human", help="organism (human/mouse)", metavar="organism")
parser.add_argument("-d", "--out", dest="out_dir", help="output directory REQUIRED", metavar="DIRECTORY", required=True)
parser.add_argument("-s", "--steps-to-run", dest="steps", help="steps of pipeline to run")
parser.add_argument("-rg", "--reference_genome", dest="reference_genome", help="steps of pipeline to run")
parser.add_argument("-m", "--merge_files", dest="partner_dir", default = "none", help="dataset with which to merge", metavar="MERGE")
parser.add_argument("--overwrite", dest="overwrite")

options = parser.parse_args()
fastq_SE_location = options.fse
print fastq_SE_location
adapter_sequence = options.adapter#.split(",")
output_directory_root = options.out_dir
if(options.partner_dir == "none"):
	merge_partner_dir = "no twin directory"
else:
	merge_partner_dir = options.partner_dir
if(options.steps == "All" or None):
	steps_to_process = steps_to_process_all
else:
	steps_to_process = options.steps.split(",")

print "Will run steps:", steps_to_process

num_threads = "7"

# *********************************************************************
# DEFINITION OF PATHS
gtf_juncs        = "gtf_juncs"
stringtie_binary = "stringtie"
tophat_binary    = "tophat"
hisat2_binary    = "hisat2"
MergeSamFiles = ["java", "-jar", "/usr/share/java/picard.jar", "MergeSamFiles"]
MarkDuplicates = ["java", "-jar", "/usr/share/java/picard.jar", "MarkDuplicates"]
CollectRnaSeqMetrics = ["java", "-jar", "/usr/share/java/picard.jar", "CollectRnaSeqMetrics"] 
prep_reads_binary= "prep_reads"
samtools_binary  = "samtools"
cutadapt_binary  = "cutadapt"
kallisto_binary = "kallisto"
Mosdepth = ["mosdepth"]
trimmomatic_cmd  = ["java", "-jar", "/home/skevin/TOOLS/Trimmomatic-0.36/trimmomatic-0.36.jar"]
min_base_quality = "15"
min_read_length  = "25"
sliding_window   = "6:20" # window_size:min_average_quality
#adapter_sequence = "AAGCAGTGGTATCAA"
cuffmerge_assembly_list = output_directory_root+"/cuffmerge_assembly_list.txt"

if(options.organism=="human"):
	if(options.reference_genome in ("hg19", "grch37")):
		hg19_gtf           = "/home/skevin/storage/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
		bowtie2index_location   = "/home/skevin/storage/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
		grch37_hisat2_reference   = "/home/skevin/storage/Homo_sapiens/grch37_tran/genome_tran"
		hg19_hisat2_reference   = "/home/skevin/storage/Homo_sapiens/UCSC/hg19/Sequence/hisat2-index/genome"
		exonic_bed_hg19 = "~/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes_exonic.bed"
		exonic_bed = exonic_bed_hg19
		hisat2_reference_genome = hg19_hisat2_reference
		gtf_location = hg19_gtf
	elif(options.reference_genome in ("hg38", "grch38")):
		bowtie2index_grch38     = "/home/skevin/storage/Homo_sapiens/grch38_NCBI_bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
		grch38_gtf                = "/home/skevin/storage/Homo_sapiens/grch38_tran/Homo_sapiens.GRCh38.87.gtf"
		grch38_hisat2_reference   = "/home/skevin/storage/Homo_sapiens/grch38_tran/genome_tran"
		hg38_reference_genome   = "/home/skevin/Homo_sapiens/grch38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
		exonic_bed_hg38 = "~/Homo_sapiens/grch38_tran/Homo_sapiens.GRCh38.87.exonic.bed"
		exonic_bed = exonic_bed_hg38
		hisat2_reference_genome = grch38_hisat2_reference
		gtf_location = grch38_gtf

elif(options.organism=="mouse"):
	grcm38_gtf = "/home/skevin/storage/Mus_musculus/NCBI/GRCm38/Annotation/Genes/genes.gtf" 
	grcm38_reference_genome = "/home/skevin/storage/Mus_musculus/NCBI/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
	hisat2_reference_genome = "/home/skevin/Mus_musculus/grcm38_tran/genome_tran"
	gtf_location = grcm38_gtf
	bowtie2index_location = "/home/skevin/storage/Mus_musculus/NCBI/GRCm38/Sequence/Bowtie2Index/genome" # /media/skevin/storage/Mus_musculus/UCSC/GRCm38/Sequence/Bowtie2Index/genome
else:
	print 'organism not recognised! Please use "human" or "mouse"'

transcriptome_index_location = gtf_location.rsplit(".",1)[0]

SE_base_filename = fastq_SE_location.replace(".fastq.gz","").split("/")[-1]
fastq_SE = fastq_SE_location.split("/")[-1]
cell_name = SE_base_filename.split("_")[1]

print cell_name
cell_output_dir = output_directory_root+"/"+cell_name+"/"
if not os.path.exists(cell_output_dir):
		os.makedirs(cell_output_dir)

# /DEFINITION OF PATHS
# *********************************************************************

# FASTQC
# =====================================================================
fastqc_output_dir = output_directory_root +"fastqc/"
if("fastqc" in steps_to_process):
	if not os.path.exists(fastqc_output_dir):
		os.makedirs(fastqc_output_dir)
	cmd = ["fastqc",fastq_SE_location] #,fastq_SE_location
	cmd += ["-t", "3"]
	cmd += ["-o", fastqc_output_dir]
	print " ".join(cmd)
	subprocess.call(cmd)
	del cmd

# CUTADAPT
# =====================================================================
#cutadapt -b AAGCAGTGGTATCAA -B AAGCAGTGGTATCAA -m 20 --overlap 5 -o DC-0001_S1_R1_001.trimmed.fastq --paired-output DC-0001_S1_R2_001.trimmed.fastq DC-0001_S1_R1_001.fastq.gz DC-0001_S1_R2_001.fastq.gz > DC-0001.trim.log
cutadapt_output_dir = cell_output_dir +"cutadapt/"
fastq_SE_cutadapt = cutadapt_output_dir+ SE_base_filename + ".cutadapt.fastq"
if("cutadapt" in steps_to_process):
	if not os.path.exists(cutadapt_output_dir):
		os.makedirs(cutadapt_output_dir)
	out_file = cell_output_dir +"cutadapt/"+cell_name+"_cutadapt.log"
	print "fastq_SE_cutadapt:",fastq_SE_cutadapt
	print "out_file: ",out_file
	cutadapt_arguments = ["-a", adapter_sequence, "-m", "20", "--overlap", "5"]
	cmd = [cutadapt_binary]
	cmd += cutadapt_arguments
	cmd += ["-o",fastq_SE_cutadapt]
	cmd += [fastq_SE_location]
	print " ".join(cmd)
	with open(out_file,"w")as outf:
		subprocess.call(cmd,stdout=outf)
	del cmd

# TRIMMOMATIC
# =====================================================================
trimmomatic_output_dir = cell_output_dir + "trimmomatic/"

fastq_SE_trimmed = trimmomatic_output_dir+ SE_base_filename + "_trimmed.fastq"

if("trimmomatic" in steps_to_process):
	if not os.path.exists(trimmomatic_output_dir):
		os.makedirs(trimmomatic_output_dir)
	print "running trimmomatic"
	cmd = trimmomatic_cmd[:]
	#print cmd
	cmd.append("SE")
	cmd.append("-threads")
	cmd.append(str(num_threads))
	#cmd.append("-trimlog")
	#cmd.append(log)
	cmd.append(fastq_SE_cutadapt)
	cmd.append(fastq_SE_trimmed)
	cmd.append("LEADING:"+min_base_quality)
	cmd.append("TRAILING:"+min_base_quality)
	cmd.append("SLIDINGWINDOW:"+sliding_window)
	cmd.append("MINLEN:"+min_read_length)
	print " ".join(cmd)
	subprocess.call(cmd)
	del cmd

# RESYNCHRONIZE FASTQ	
# =====================================================================
if("resync_fastq" in steps_to_process):

	cmd = ["./src/fastqCombinePairedEnd.py"]
	cmd += [fastq_r1_trimmed, fastq_r2_trimmed] #not fixed yet. Unlikely to be used again
	cmd += ["None"]
	print " ".join(cmd)
	subprocess.call(cmd)

# HISAT2[*tweaked]
# =====================================================================
# hisat2 -x /home/skevin/storage/Homo_sapiens/grch38/genome -1 76/trimmomatic/Hu_76_S127_R1_001.trimmed.fastq -2 76/trimmomatic/Hu_76_S127_R2_001.trimmed.fastq -U 76/trimmomatic/Hu_76_S127_R1_001.trimmed.unpaired.fastq -U 76/trimmomatic/Hu_76_S127_R2_001.trimmed.unpaired.fastq -S Hu_76_hisat2.sam
hisat2_output_dir = cell_output_dir + "hisat2/"
hisat2_sorted_bam = hisat2_output_dir + cell_name + ".bam"
if not os.path.exists(hisat2_output_dir):
	os.makedirs(hisat2_output_dir)
if("hisat2" in steps_to_process):
	print "processing hisat2-tweaked"
	log_file = hisat2_output_dir +cell_name+".hisat2.log"
	cmd = [hisat2_binary]
	cmd += ["--threads", num_threads]
	cmd += ["--new-summary"]
	cmd += ["--pen-noncansplice", "20"] #<PENALITY_NONCANONICAL> 0;3;12;20 Default = 3
	cmd += ["--mp", "1,0"] #<MAX_MISMATCH_PENALITY> 1; 2; 3; 6 Default = 6, <MIN_MISMATCH_PENALITY> 0;1;2 Default = 2
	cmd += [ "--sp", "3,1"] #<MAX_SOFTCLIPPING_PENALITY> 1;2;3 Default = 2, <MIN_SOFTCLIPPING_PENALITY> 0;1;2 Default = 1 
	cmd += ["-x",hisat2_reference_genome]
	cmd += ["-U",fastq_SE_trimmed]
		#cmd += ["-S",hisat2_output_sam]
	print datetime.datetime.today()
	print " ".join(cmd)
	with open(hisat2_sorted_bam,"wb") as sorted_bam, open(log_file, "wb") as outerr:
		hisat2 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=outerr)
		sort = subprocess.Popen([samtools_binary,"sort", "-O", "bam", "-T", hisat2_output_dir+"samtools_sort.tmp"], stdout=sorted_bam, stdin=hisat2.stdout, stderr=outerr)
		hisat2.stdout.close()
		hisat2.wait()
		sort.wait()
		sort.communicate()
	print datetime.datetime.today()
	del cmd	


#~ # MERGE FILES
#~ # =====================================================================
#~ #Merge I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
picard_output_dir = cell_output_dir+"picard/"
merge_partner_bam = merge_partner_dir+ cell_name + "hisat2" + SE_base_filename + ".bam"
merged_bam= hisat2_output_dir+ cell_name + "_merged.bam"
merge_log = picard_output_dir+SE_base_filename+"_merged.txt"
picard_tmp_dir = picard_output_dir+"tmp/"

if((("merge_files" in steps_to_process) or ("all" in steps_to_process)) & (merge_partner_dir != "none")):
	merged_bam= hisat2_output_dir+ cell_name + "_merged.bam"
	if not os.path.isdir(picard_output_dir):
		os.makedirs(picard_output_dir)
	cmd = MergeSamFiles[:]
	cmd.append("I=")
	cmd.append(merge_partner_bam)
	cmd.append("I=")
	cmd.append(hisat2_sorted_bam)
	cmd.append("O=")
	cmd.append(merged_bam)
	out_file = picard_output_dir + SE_base_filename+"_merged.log"
	err_file = picard_output_dir + SE_base_filename+"_merged.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
			errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "mergesamfiles finished successfully"
			else:
				print "mergesamfiles failed !!!!"
	del cmd

# MARKDUPLICATES (and remove)
# =====================================================================
#MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
if not os.path.isdir(picard_output_dir):
	os.makedirs(picard_output_dir)
picard_output_dir = cell_output_dir+"picard/"
removed_duplicates_bam= picard_output_dir+ SE_base_filename + "_removed_duplicates.bam"
removed_dup_metrics = picard_output_dir+SE_base_filename+"_removed_dup_metrics.txt"
picard_tmp_dir = picard_output_dir+"tmp/"

if(("removeduplicates" in steps_to_process) or ("all" in steps_to_process)):
	print merged_bam
	if not os.path.isdir(picard_output_dir):
		os.makedirs(picard_output_dir)
	cmd = MarkDuplicates[:]
	cmd.append("I=")
	cmd.append(hisat2_sorted_bam)
	cmd.append("O=")
	cmd.append(removed_duplicates_bam)
	cmd.append("M=")
	cmd.append(removed_dup_metrics)
	cmd.append("REMOVE_DUPLICATES=true")
	out_file = picard_output_dir + SE_base_filename+"_markduplicates.log"
	err_file = picard_output_dir + SE_base_filename+"_markduplicates.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
			errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "markduplicates finished successfully"
			else:
				print "markduplicates failed !!!!"
	del cmd

#~ # COLLECT RNASEQ METRICS 
#~ # =====================================================================
#~ #MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
#~ picard_output_dir = cell_output_dir+"picard/"
#~ removed_duplicates_bam= picard_output_dir+ SE_base_filename + "_removed_duplicates.bam"
#~ rnaseq_metrics = picard_output_dir+SE_base_filename+"_rna_metrics.txt"
#~ picard_tmp_dir = picard_output_dir+"tmp/"
#~ ref_flat = os.path.expanduser("~/Homo_sapiens/grch38_tran/grch38_refFlat.txt")
#~ ribosomal_interval_list = os.path.expanduser("~/Homo_sapiens/grch38_tran/intervalList.txt")

#~ if(("rnaseq_metrics" in steps_to_process) or ("all" in steps_to_process)):
	#~ if not os.path.isdir(picard_output_dir):
		#~ os.makedirs(picard_output_dir)
	#~ cmd = CollectRnaSeqMetrics[:]
	#~ cmd += ["I=", removed_duplicates_bam]
	#~ cmd += ["O=", rnaseq_metrics]
	#~ cmd += ["REF_FLAT=", ref_flat]
	#~ cmd += ["STRAND=NONE"]
	#~ cmd += ["RIBOSOMAL_INTERVALS=", ribosomal_interval_list]
	#~ print " ".join(cmd)
	#~ with open(err_file,"w")as outerr:
		#~ with open(out_file,"w")as outf:
			#~ errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			#~ if(errcode == 0):
				#~ print "collectrnaseqmetrics finished successfully"
			#~ else:
				#~ print "collectrnaseqmetrics failed !!!!"
	#~ del cmd

# RSEQC 
# =====================================================================
#MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
picard_output_dir = cell_output_dir+"picard/"
removed_duplicates_bam= picard_output_dir+ SE_base_filename + "_removed_duplicates.bam"
picard_tmp_dir = picard_output_dir+"tmp/"
ref_gene_model = os.path.expanduser("~/Homo_sapiens/grch38_tran/Homo_sapiens.GRCh38.87.bed")
rseqc_out = picard_output_dir+SE_base_filename
if(("rseqc" in steps_to_process) or ("all" in steps_to_process)):
	if not os.path.isdir(picard_output_dir):
		os.makedirs(picard_output_dir)
	genebody_coverage = ["geneBody_coverage.py", "-i", removed_duplicates_bam, "-r", ref_gene_model, "-o", rseqc_out] 
	inner_distance = ["inner_distance.py", "-i", removed_duplicates_bam, "-o", rseqc_out, "-r", ref_gene_model]
	junction_saturation =  ["junction_saturation.py", "-i", removed_duplicates_bam, "-o", rseqc_out, "-r", ref_gene_model]
	read_distribution = ["read_distribution.py", "-i", removed_duplicates_bam, "-r", ref_gene_model]
	out_file = picard_output_dir + SE_base_filename+"_rnaseq_metrics.log"
	err_file = picard_output_dir + SE_base_filename+"_rnaseq_metrics.err"
	read_distout = picard_output_dir + SE_base_filename+"_read_dist.log"
	read_disterr = picard_output_dir + SE_base_filename+"_read_dist.err" 
	print datetime.datetime.today()
	print " ".join(read_distribution)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
			with open(read_distout, "w")as rdout:
				with open(read_disterr, "w")as rderr:
					#subprocess.call(genebody_coverage,stdout=outf,stderr=outerr)
					subprocess.call(inner_distance, stdout=outf, stderr=outerr)
					subprocess.call(junction_saturation, stdout=outf, stderr=outerr)
					subprocess.call(read_distribution, stdout=rdout, stderr=rderr)
					#~ if(errcode == 0):
						#~ print "rseqc finished successfully"
					#~ else:
						#~ print "rseqc failed !!!!"
	del genebody_coverage
	del inner_distance
	del junction_saturation
	del read_distribution



# GET UNIQUE INTRONIC READS
# =====================================================================
exonic_bed  = '/dataVolume/storage/Homo_sapiens/grch38_tran/Homo_sapiens.GRCh38.87.exonic.bed'
# hisat2 -x /home/thor/storage/Homo_sapiens/grch38/genome -1 76/trimmomatic/Hu_76_S127_R1_001.trimmed.fastq -2 76/trimmomatic/Hu_76_S127_R2_001.trimmed.fastq -U 76/trimmomatic/Hu_76_S127_R1_001.trimmed.unpaired.fastq -U 76/trimmomatic/Hu_76_S127_R2_001.trimmed.unpaired.fastq -S Hu_76_hisat2.sam
if not os.path.exists(hisat2_output_dir):
 	os.makedirs(hisat2_output_dir)
intronic_bam = picard_output_dir + SE_base_filename + "_unique_intronic.bam"
if("uniqify" in steps_to_process):
	print "getting unique reads"
	if not os.path.exists(hisat2_output_dir):
		os.makedirs(hisat2_output_dir)
	cmd = ["/home/skevin/single_cell_pipeline/src/kstachelek_src/bioinf_oneliners/grep_unique_reads.sh"]
	cmd += [hisat2_sorted_bam]
	print datetime.datetime.today()
	print " ".join(cmd)
	err_file = hisat2_output_dir + SE_base_filename+"_uniqify.log"
	with open(intronic_bam,"wb")as outf, open(err_file, "wb") as outerr:
		uniq = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=outerr)
		get_intron = subprocess.Popen(["bedtools", "intersect", "-wa", "-abam", "-", "-b", exonic_bed, "-v"], stdin=uniq.stdout, stdout=outf, stderr=outerr)
		uniq.wait()
		get_intron.wait()
		get_intron.communicate()
	del cmd

# THRESHOLD BAM BELOW DEPTH N PILEUPS
# =====================================================================

# hisat2 -x /home/thor/storage/Homo_sapiens/grch38/genome -1 76/trimmomatic/Hu_76_S127_R1_001.trimmed.fastq -2 76/trimmomatic/Hu_76_S127_R2_001.trimmed.fastq -U 76/trimmomatic/Hu_76_S127_R1_001.trimmed.unpaired.fastq -U 76/trimmomatic/Hu_76_S127_R2_001.trimmed.unpaired.fastq -S Hu_76_hisat2.sam

intronic_bam = picard_output_dir + SE_base_filename + "_unique_intronic.bam"
threshold_bam = hisat2_output_dir + SE_base_filename + "_depth_under_3.bam"
threshold_bed  = hisat2_output_dir + SE_base_filename + "_depth_over_3.bed"

if("threshold_bam" in steps_to_process):
	print "filtering bam at max pileup (default 3)"
	if not os.path.exists(hisat2_output_dir):
		os.makedirs(hisat2_output_dir)
	cmd = ["samtools", "depth", intronic_bam]
	print datetime.datetime.today()
	print " ".join(cmd)
	err_file = hisat2_output_dir + SE_base_filename+"_threshold_bam.log"
	with open(threshold_bam,"wb")as outbam, open(threshold_bed, "wb") as outbed, open(err_file, "wb") as outerr:
		print " ".join(["samtools", "depth", intronic_bam])
		samtools_depth = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=outerr)
		bed_above_N = subprocess.Popen(["/home/skevin/single_cell_pipeline/src/bed_from_areas_covered_above_N_v2.py", "3"], stdin=samtools_depth.stdout, stdout=subprocess.PIPE, stderr=outerr)
		bam_below_N = subprocess.Popen(["bedtools", "intersect", "-wa", "-abam", intronic_bam, "-b", "-", "-v"], stdin=bed_above_N.stdout, stdout=outbam, stderr=outerr)
		samtools_depth.wait()
		bed_above_N.wait()
		bam_below_N.wait()
		bam_below_N.communicate()
	print datetime.datetime.today()
	
	del cmd

# STRINGTIE
# =====================================================================
stringtie_output_dir = output_directory_root + "/stringtie/"
stringtie_cell_dir = stringtie_output_dir + cell_name+"/"
stringtie_gtf = stringtie_cell_dir + cell_name+".gtf"
if("stringtie" in steps_to_process):
	if not os.path.exists(stringtie_output_dir):
		os.makedirs(stringtie_output_dir)
	if not os.path.exists(stringtie_cell_dir):
		os.makedirs(stringtie_cell_dir)
	print "processing stringtie"
	cmd = [stringtie_binary]
	cmd += ["-G", gtf_location] 
	cmd += ["-x", "MT"]
	cmd += ["-eB"]
	cmd += ["-p",num_threads]
	cmd += ["-o",stringtie_gtf]
	cmd += [removed_duplicates_bam]
	out_file = stringtie_cell_dir + "stringtie.log"
	err_file = stringtie_cell_dir + "stringtie.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	subprocess.call(cmd,stdout=outf,stderr=outerr)
	del cmd

# MOSDEPTH COVERAGE
# =====================================================================
#MOSDEPTH
coverage_table= picard_output_dir + SE_base_filename + "_coverage.txt"
coverage_dist= picard_output_dir+ SE_base_filename + "_cumul.dist"

if(not os.path.isfile(coverage_table) and "mosdepth" in steps_to_process):
	cmd = Mosdepth[:]
	cmd += ["--distribution", coverage_dist]
	cmd += [removed_duplicates_bam]
	print " ".join(cmd)
	log_file = picard_output_dir +SE_base_filename+"_mosdepth.log"
	with open(log_file,"w")as logf:
		with open(coverage_table, "w") as outf:
			errcode = subprocess.call(cmd, stdout = outf)
	del cmd
	
#~ PLOT ALIGNMENT STATS
#~ =========================================================================


#~ # KALLISTO
#~ # =====================================================================
#~ kallisto_output_dir = output_directory_root + "/kallisto/"
#~ kallisto_cell_dir = kallisto_output_dir + cell_name+"/"
#~ kallisto_spike_index = "/home/skevin/single_cell_pipeline/stats/kallisto_spike_index"
#~ if("kallisto" in steps_to_process):
	#~ if not os.path.exists(kallisto_output_dir):
		#~ os.makedirs(kallisto_output_dir)
	#~ if not os.path.exists(kallisto_cell_dir):
		#~ os.makedirs(kallisto_cell_dir)
	#~ print "processing kallisto"
	#~ cmd = [kallisto_binary]
	#~ cmd += ["quant"] 
	#~ cmd += ["-i",kallisto_spike_index]
	#~ cmd += ["-o",kallisto_cell_dir]
	#~ cmd += [fastq_r1_trimmed]
	#~ cmd += [fastq_r2_trimmed]
	#~ out_file = kallisto_cell_dir + "kallisto.log"
	#~ err_file = kallisto_cell_dir + "kallisto.err"
	#~ print " ".join(cmd)
	#~ with open(err_file,"w")as outerr:
		#~ with open(out_file,"w")as outf:
		 	#~ subprocess.call(cmd,stdout=outf,stderr=outerr)
	#~ del cmd
	


