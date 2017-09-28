import argparse
import os
import sys
from lib.dephase_vcf import dephase_vcf
create_folder = lambda f: [ os.makedirs(f) if not os.path.exists(f) else False ]

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Takses phased vcf and dephases it')
	parser.add_argument('--vcf', metavar='<VCF>', help='VCF to be dephased')
	parser.add_argument('--out', metavar='<FOLDER>', help='folder to store data')
	args = parser.parse_args()
	vcf_name = args.vcf

	if not os.path.isfile(vcf_name):
		print("No such file: '{}'".format(vcf_name), file=sys.stderr)

	create_folder(args.out)
	dephase_vcf(vcf_name, args.out)