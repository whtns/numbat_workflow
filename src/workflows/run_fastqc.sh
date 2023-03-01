#!/bin/bash

find ../sc_RB_devel/20170407_sunlee_H_sapiens/ -name "*.fastq.gz" | parallel -j 6 fastqc

