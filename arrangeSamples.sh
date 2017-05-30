popFolder="PopInds"

#Make sure the vcf file has snp id's of the form chrom:pos if it is simulated and the separator is tab
#zcat myvcf.vcf.gz | awk '{if ($3==".") $3=$1":"$2; print$0}'
#separate head and body to make space into tabs


zcat $1 |head -20 | grep '#CHROM' | cut -f10- | sed 's/\t/\n/g' > temp.names
lines=`cat temp.names | wc -l`
echo "number of individuals is:"
echo $lines
s=$(( lines / 3 ))
split -l$s temp.names

rm temp.names
mv xaa $popFolder/YRI
mv xab $popFolder/CEU
mv xac $popFolder/CHB

python src/generate_vcfs.py --simulated
python src/generate_impute_results.py --runonly= --arrays=affy,omni,illumina --types=c --steps=3 --dry_run 0 --phased

#python src/phase_impute_pipeline.py --step 1


#python /sb/home/apoursha/imputation/compare_impute/src/phase_impute_pipeline.py --ind_ref=/home/apoursha/imputation/compare_impute/samplesize/test=5/ref=5/ref.inds --ind_count=10 --ref_gzvcf=/home/apoursha/imputation/compare_impute/samplesize/test=5/ref=5/ch22-ref.vcf.gz --test_gzvcf=/home/apoursha/imputation/compare_impute/samplesize/test=5/ref=5/ch22-test_affy.vcf.gz --step=1 --logs=/home/apoursha/imputation/compare_impute/samplesize/test=5/ref=5/imputation/step1-ch22-affy/logs --script_log=/home/apoursha/imputation/compare_impute/samplesize/test=5/ref=5/imputation/step1-ch22-affy/logs/script.log --script_dir=/sb/home/apoursha/imputation/compare_impute/src --job_name_prefix=samplesize-0005-0005-ch22-run-01 --intermediate_dir=/home/apoursha/imputation/compare_impute/samplesize/test=5/ref=5/imputation/step1-ch22-affy/tmp --outdir=/home/apoursha/imputation/compare_impute/samplesize/test=5/ref=5/imputation/step1-ch22-affy/results --completed_file=/home/apoursha/imputation/compare_impute/samplesize/test=5/ref=5/imputation/step1-ch22-affy/results/completed.flag --chr=22 --qsub_id_file=/home/apoursha/imputation/compare_impute/samplesize/test=5/ref=5/imputation/step1-ch22-affy/logs/job_id.txt --qsub_use_testq=0 --qsub_run_locally=1
