import subprocess
import os
import argparse
import pandas as pd
import numpy as np

LD_SUBSAMPLER=3

# no idea where these numbers came from
# @TODO: ask @apoursh to fill me in on where they name from
PANELS = {}
PANELS['omni2.5', 1] = int(2388927/60/25)
PANELS['affy6', 1]   = int(885501/60/25)
PANELS['illumina3', 1] = int(308330/60/25)


create_folder = lambda f: [ os.makedirs(f) if not os.path.exists(f) else False ]

def make_simulated_panel(vcf, outfolder, vcftools, size=PANELS['affy6', 1], panel_name='affy'):
	"""
	:param vcf: location of the vcf
	:param outfolder: location to save outputs to
	:param vcftools: location of vcftools
	"""
	vcf_file_name = os.path.basename(vcf)

	af_file_name =  os.path.join(outfolder, 'panel' + panel_name, vcf_file_name + ".frq")

	outfolder = os.path.join(outfolder, "panel" + panel_name)

	create_folder(outfolder)
	panelfile = os.path.join(outfolder, "simulated_" + panel_name +'.txt')

	if os.path.exists(panelfile):
		print("Output already exists, skipping")
		return None
	# check if allele frequencies exist, if not regenerate them

	command = "{vcftools} --gzvcf {vcf} --freq2 --out {outfile}"


	if not os.path.exists(af_file_name):
		subprocess.Popen(command.format(vcftools=vcftools,
			vcf=vcf,
			outfile=af_file_name[:-4]),
		shell=True).wait()

	# wait for jobs to finish.


	freq = pd.read_table(filepath_or_buffer=af_file_name, sep='\t', header=0)
	chrom = np.array(freq.index)
	pos = np.array(freq.CHROM)
	freq = np.array(freq['{FREQ}'])
	freq[freq > 0.5] = 1.0 - freq[freq > 0.5] # get minor allele frequency
	freq = freq / np.sum(freq)
	locations = np.random.choice(range(freq.shape[0]),
								 size=LD_SUBSAMPLER * size,
								 replace=False,
								 p = freq)
	locations = locations[0::LD_SUBSAMPLER]
	locations = np.sort(locations)
	del freq
	with open(panelfile, 'w') as fp:
		for location in locations:
			fp.write("""%(ch)s:%(pos)s\n""" % {
				'ch' : chrom[location],
				'pos' : pos[location]})


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Takses phased vcf and dephases it')
	parser.add_argument('--vcf', metavar='<VCF>', help='VCF to be dephased')
	parser.add_argument('--vcftools', metavar='<VCFTOOLS>', help='location of vcftools')
	parser.add_argument('--out', metavar='<FOLDER>', help='folder to store data')
	args = parser.parse_args()
	create_folder(args.out)
	make_simulated_panel(args.vcf,
					     args.out,
					     args.vcftools)

