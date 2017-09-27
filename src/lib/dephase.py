import pyvcf
import argparse
import os.path


def main():
    parser = argparse.ArgumentParser(description='Takses phased vcf and dephases it')
    parser.add_argument('--vcf', metavar=<VCF>, help='VCF to be dephased')
    args = parser.parse_args()
    vcf_name = args.vcf

    if not os.path.isfile(vcf_name):
        print("No such file: '{}'".format(vcf_name), file=sys.stderr)


if __name__=='__main__':

