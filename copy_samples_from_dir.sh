#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -A ams-754-aa
#PBS -q lm
#PBS -N copy-data-and-prep
#PBS -o msprime-pipe-prep.out
#PBS -e msprime-pipe-prep.err

source_folder="/sb/project/ams-754-aa/zaf_projects/deep/data/"
target_folder="/sb/project/ams-754-aa/zaf_projects/compare_impute/data/"

touch ${source_folder}msprime-gz.inprogress
gzip ${source_folder}msprime-data166668.vcf
rm ${source_folder}msprime-gz.inprogress
touch ${source_folder}msprime-copy.inprogress
mkdir ${target_folder}
cp ${source_folder}msprime-data166668.vcf.gz ${target_folder}
rm ${source_folder}msprime-copy.inprogress
touch ${source_folder}msprime-preparing.inprogress
bash /sb/project/ams-754-aa/zaf_projects/compare_impute/prep_samples_for_pipeline.sh ${target_folder}msprime-data166668.vcf.gz
rm ${source_folder}msprime-preparing.inprogress
touch ${source_folder}msprime-preparing.finished
