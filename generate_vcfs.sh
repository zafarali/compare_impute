#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -A ams-754-aa
#PBS -q sw
#PBS -N generate-vcfs
#PBS -o generate-vcfs.out
#PBS -e generate-vcfs.err

datafolder="/sb/project/ams-754-aa/zaf_projects/compare_impute/data/"
vcffile="/sb/project/ams-754-aa/zaf_projects/compare_impute/data/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
vcftools="/cvmfs/soft.mugqic/CentOS6/software/vcftools/vcftools-0.1.14/bin/vcftools"
cd $PBS_O_WORKDIR


#### CREATE REFERENCE FILE
${vcftools} --gzvcf ${vcffile} --keep ${datafolder}ref.inds --remove-indels --recode --stdout | gzip -c - > ${datafolder}ref.vcf.gz.inprogress &&
cd .. &&
mv ${datafolder}ref.vcf.gz.inprogress ${datafolder}ref.vcf.gz &&


### CHUNK THIS VCF
MIN_AND_MAX=$(zcat ${vcffile} | grep -v -E "^#" | awk '{print $2}' | sort |  awk 'NR==1; END{print}')
MAX=$(echo "$MIN_AND_MAX" | tail -n 1)
MIN=$(echo "$MIN_AND_MAX" | head -n 1)

echo "Using min:${MIN} and max:${MAX}"

# Calculates the actual size of a chunk in each VCF
CORESIZE=5000000 # 5Mb region
CHUNKSIZE=$(python -c "from math import ceil; range=${MAX}-${MIN}; chunks=ceil(range/float(${CORESIZE})); print(int(ceil(range / chunks)))") &&


SPLITTER="/sb/project/ams-754-aa/zaf_projects/imputation/data/splitVCFref.jar"

java -Xmx4g -jar $SPLITTER --vcfref $VCFFILE --coresize ${CHUNKSIZE}.b


#### NOW CREATE STUDY GROUND TRUTH.
${vcftools} --gzvcf ${vcffile} --keep ${datafolder}test.inds --remove-indels --recode --stdout | gzip -c - > ${datafolder}study_truth.vcf.gz.inprogress &&
mv ${datafolder}study_truth.vcf.gz.inprogress ${datafolder}study_truth.vcf.gz &&


#### NOW CREATE THE AFFY STUDY FILE
${vcftools} --gzvcf ${vcffile} --snps ${datafolder}affy6_markers.txt --remove-indels --recode --stdout | gzip -c - > ${datafolder}study_affy.vcf.gz.inprogress &&
mv ${datafolder}study_affy.vcf.gz.inprogress ${datafolder}study_affy.vcf.gz
