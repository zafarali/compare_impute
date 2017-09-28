from optparse import OptionParser
import os
import glob
import numpy as np
import random
from lib.dephase_vcf import dephase_vcf

create_folder = lambda f: [ os.makedirs(f) if not os.path.exists(f) else False ]
HOME = test_dir = os.path.realpath(os.path.join(__file__, "..", '..'))

parser = OptionParser()

parser.add_option('--vcf', default='/sb/project/ams-754-aa/zaf_projects/compare_impute/data/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
parser.add_option('--outfolder', default=os.path.join(HOME, 'results'))
parser.add_option('--vcftools', default='/cvmfs/soft.mugqic/CentOS6/software/vcftools/vcftools-0.1.14/bin/vcftools')
parser.add_option('--java', default='/usr/bin/java')
parser.add_option('--split_vcf_jar', default=HOME + '/src/lib/splitVCFref.jar')
parser.add_option('--markers', default=HOME + '/data/affy6_markers.txt')
parser.add_option('--simulated', action='store_true')
parser.add_option('--dephase', action='store_true')
parser.add_option('--queue', default='debug')
parser.add_option('--nodes', default=16)

class Sampler(object):
	def __init__(self, vcf):
		individual_ids = np.genfromtxt(vcf, dtype=str, max_rows=1, comments='##').tolist()
		self.all_inds = individual_ids[9:]
	def train_test_split(self, train_count):
		train_count = min(train_count,len(self.all_inds))
		random.shuffle(self.all_inds)
		return self.all_inds[:train_count], self.all_inds[train_count:]

def write_individuals_file(filename, individuals):
	with open(filename, 'w') as f:
		for i in individuals:
			f.write(str(i)+'\n')

def create_individuals_files(options):
	sampler = Sampler(options.vcf)
	train_inds, test_inds = sampler.train_test_split(train_count=1204)
	reference_file = os.path.join(options.outfolder, 'ref.inds')
	study_file = os.path.join(options.outfolder, 'study.inds')
	write_individuals_file(reference_file, train_inds)
	write_individuals_file(study_file, test_inds)


def main():
	pass

def create_vcf(vcf, arguments, outname, options):
	template="""
	{vcftools} --gzvcf {vcf} {vcftools_args} --recode --stdout | gzip -c - > {outname}.inprogress
	mv {outname}.inprogress {outname}"""
	return template.format(vcftools_args=arguments,
						   outname=outname,
						   vcf=options.vcf,
						   vcftools=options.vcftools)

def split_vcfs(vcf, options, coresize=5000000):
	splitter="""
	bash {home}/splitter.sh {vcf}"""
	return splitter

def simulate_affy_panels(options):
	cmd = """
	python3.5 {home}/src/simulate_affy_panels.py --vcf {vcf} --out {outfolder} --vcftools {vcftools}"""
	return cmd

TEMPLATE = """#!/bin/bash
#PBS -l nodes=1:ppn={nodes}
#PBS -l walltime=10:00:00
#PBS -A ams-754-aa
#PBS -N pipeline-{vcf}
#PBS -o {outfolder}/out.txt
#PBS -e {outfolder}/error.txt
#PBS -q {queue}
"""
if __name__ == '__main__':
	(options,args) = parser.parse_args()


	qsub_script = TEMPLATE

	create_folder(options.outfolder)
	create_individuals_files(options)

	if options.dephase:
		qsub_script += '{home}/src/dephase.py --vcf {vcf} --out {outfolder}'
		options.vcf = os.path.join(options.outfolder, '*.dephased.*')
	qsub_script += create_vcf(options.vcf,
							  '--keep {outfolder}/ref.inds --remove-indels',
							  '{outfolder}/ref.vcf.gz',
							  options)
	qsub_script += create_vcf(options.vcf,
							  '--keep {outfolder}/test.inds --remove-indels',
							  '{outfolder}/study_ground_truth.vcf.gz',
							  options)

	if options.simulated:
		qsub_script += simulate_affy_panels(options)
		options.markers = os.path.join(options.outfolder, 'panelaffy', 'simulated_affy.txt')

	qsub_script += create_vcf('{outfolder}/study_ground_truth.vcf.gz',
							  '--snps {markers} --remove-indels',
							  '{outfolder}/study.vcf.gz',
							  options)
	qsub_script += split_vcfs('{outfolder}/ref.vcf.gz', options)


	print(qsub_script)
	print(qsub_script.format(
				home=HOME,
				outfolder=options.outfolder,
				markers=options.markers,
				vcf=options.vcf,
				nodes=options.nodes,
				queue=options.queue))
	# print(HOME)