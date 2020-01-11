#!/bin/bash -l

#$ -S /bin/bash

#$ -P ProjectName

#$ -cwd

#$ -j y

#$ -o import_qsub.log

#$ -N importfunction

#$ -l h_rt=150:00:00


module load R/3.6.0
module load gcc

Rscript ImportScript.R \
 -d raw_feature_bc_matrix \
 -c filtered_feature_bc_matrix \
 -p Preprocessing \
 -g TRUE \
 -s SampleName \
 -d Directory

