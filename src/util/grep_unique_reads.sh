#!/bin/bash

input_bam=$1 #"input bam"
output_bam=`echo $input_bam | sed 's/.bam/_intronic.bam/'`
(samtools view -H $input_bam; samtools view -F 2308 $input_bam | grep -w 'NH:i:1') | samtools view -bS - 
