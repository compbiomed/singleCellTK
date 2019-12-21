#!/bin/bash -l

#$ -S /bin/bash

#$ -P camplab

#$ -cwd

#$ -j y

#$ -o import_qsub.log

#$ -N importfunction

#$ -l h_rt=150:00:00


module load R/3.6.0
module load gcc

Rscript ImportScript.R \
 -u raw_feature_bc_matrix \
 -f filtered_feature_bc_matrix \
 -p CellRanger \
 -g TRUE \
 -s pbmc \
 -d Directory

