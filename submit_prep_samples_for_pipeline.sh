#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=24:00:00
#PBS -A ams-754-aa
#PBS -q sw
#PBS -N preparing-msprime-for-piepline
#PBS -o msprime-prep.out
#PBS -e msprime-prep.err

cd $PBS_O_WORKDIR
bash prep_samples_for_pipeline.sh 