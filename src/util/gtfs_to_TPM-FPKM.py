#!/usr/bin/python

import pandas as pd
import sys
import re

# argv[1] = output file (will create <name>.fpkm.csv and <name>.tpm.csv)
# argv[2:] = gtf files

files = sys.argv[2:]
fpkm_csv_filename = sys.argv[1]+".fpkm.csv"
tpm_csv_filename = sys.argv[1]+".tpm.csv"

RE_FPKM=re.compile('FPKM "([^"]+)"')
RE_TPM =re.compile('TPM "([^"]+)"')
RE_GENE_ID=re.compile('gene_id "([^"]+)"')
#RE_GENE_NAME=re.compile('gene_name "([^"]+)"')
RE_TRANSCRIPT_ID=re.compile('transcript_id "([^"]+)"')

#transcript2gene = {}

cells = [f.split("/")[-1].split(".")[0] for f in files]
# ~ print cells
fpkm_table = pd.DataFrame()
tpm_table  = pd.DataFrame()

for f in files:
	# ~ print f
	cell = f.split("/")[-1].split(".")[0]
	cell_fpkm = []
	cell_tpm = []
	cell_transcripts = []
	with open(f) as gtf:
		for line in gtf:
			if(line[0] == '#'):
				continue
			x=line.split("\t")
			if x[2]=="transcript":
				fpkm = RE_FPKM.search(x[8]).group(1)
				tpm = RE_TPM.search(x[8]).group(1)
				#print tpm
				#exit()
				#gene_id = RE_GENE_ID.search(x[8])
				transcript_id = RE_TRANSCRIPT_ID.search(x[8]).group(1)
				#transcript2gene[transcript_id]=gene_id
				cell_fpkm.append(fpkm)
				cell_tpm.append(tpm)
				cell_transcripts.append(transcript_id)
	cell_fpkm_frame = pd.DataFrame(cell_fpkm, index=cell_transcripts, columns=[cell])
	fpkm_table = fpkm_table.join(cell_fpkm_frame, how="outer")
	
	cell_tpm_frame = pd.DataFrame(cell_tpm, index=cell_transcripts, columns=[cell])
	tpm_table = tpm_table.join(cell_tpm_frame, how="outer")

#print fpkm_table
fpkm_table.index.name = "TRANSCRIPT_ID"
tpm_table.index.name = "TRANSCRIPT_ID"


#print fpkm_table

fpkm_table.to_csv(fpkm_csv_filename, sep="\t")
tpm_table.to_csv(tpm_csv_filename, sep="\t")
