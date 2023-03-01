#!/bin/bash
./copywriter_single_cell.R -b 250000 -i ~/single_cell_pipeline/output/GT_20180109_SHL_H_sapiens_RB_31_output -r *under_3.bam$
./copywriter_single_cell.R -b 500000 -i ~/single_cell_pipeline/output/GT_20180109_SHL_H_sapiens_RB_31_output -r *under_3.bam$
./copywriter_single_cell.R -b 2000000 -i ~/single_cell_pipeline/output/GT_20180109_SHL_H_sapiens_RB_31_output -r *under_3.bam$
