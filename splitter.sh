vcf=${1}

### CHUNK THIS VCF
MIN_AND_MAX=$(zcat ${vcf} | grep -v -E "^#" | awk '{print $2}' | sort |  awk 'NR==1; END{print}')
MAX=$(echo "$MIN_AND_MAX" | tail -n 1)
MIN=$(echo "$MIN_AND_MAX" | head -n 1)

echo "Using min:${MIN} and max:${MAX}"

# Calculates the actual size of a chunk in each VCF
CORESIZE=5000000 # 5Mb region
CHUNKSIZE=$(python -c "from math import ceil; range=${MAX}-${MIN}; chunks=ceil(range/float(${CORESIZE})); print(int(ceil(range / chunks)))") &&


SPLITTER="/sb/project/ams-754-aa/zaf_projects/imputation/data/splitVCFref.jar"

java -Xmx4g -jar $SPLITTER --vcfref ${vcf} --coresize ${CHUNKSIZE}.b
