"""Comparison of imputation results"""
import argparse
import numpy as np
import pandas as pd
from cyvcf2 import VCF
import subprocess
import pdb
from sklearn.metrics import confusion_matrix



TABIX = "/cvmfs/soft.mugqic/CentOS6/software/tabix/tabix-0.2.6/tabix"
BGZIP = "/cvmfs/soft.mugqic/CentOS6/software/tabix/tabix-0.2.6/bgzip"

class imputeOutputs(object):
    @profile
    def __init__(self, snps, testFile):
        self.snps = snps
        self.i2   = None
        self.NN   = None
        cmd = """gunzip {filename} &&
        {bgzip} {unzipped} &&
        {tabix} {filename}""".format(
                tabix = TABIX,
                bgzip = BGZIP,
                unzipped = testFile[:-3],
                filename = testFile)
        vcf = VCF(testFile, gts012=True)
        snps, names = [], []
        for v in vcf:
            snps.append(v.gt_types.tolist())
            names.append(v.ID)
        try:
            subprocess.call(cmd,shell=True)
        except Exception as e:
            print("bad stuff happened!\n{}".format(e))

        self.test = pd.DataFrame(snps, index=names)

    def impute2(self, filename):
        i2 = pd.read_csv(filename, delimiter=' ', header=None)
        n = i2.shape[0]
        p = int((i2.shape[1]-5)/3)
        gt = np.array([0,1,2])
        i2.index = i2[1]
        for i in range(0,p):
            i2[i] = gt[np.argmax(np.array([i2[5+3*i],i2[6+3*i],i2[7+3*i]]), axis=0)]
        colNames = list(range(p))
        self.i2 = i2[colNames].astype(int)

    def evaluate(self, df):
        disc = np.mean(np.mean(np.abs(df-self.test)))
        cm = confusion_matrix(self.test.values.ravel(), df.values.ravel())
        return(disc, cm)

    def NN(self, filename):
        pass




def main():
    parser = argparse.ArgumentParser(description='compare two sets of imputation results')
    parser.add_argument('--snps', help="File with relevant snp names", required=True)
    parser.add_argument('--test', help="Test vcf (actual answers")
    parser.add_argument('--i2_file', help="path the impute2 results file", required=True)
    parser.add_argument('--NN_file', help="path to NN imputation results")
    args = parser.parse_args()

    outputs = imputeOutputs(args.snps, args.test)
    outputs.impute2(args.i2_file)
    disc, cm = outputs.evaluate(outputs.i2)
    print(disc)
    print(cm)
    df = outputs.i2.copy()
if __name__=='__main__':
    main()
