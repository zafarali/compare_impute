"""
generate_vcfs.py
author: Gerard Tse (gerardtse@gmail.com)

Generates test and reference sets in the format of VCF from the 1000 Genome data. This script
generates sets sample size and population based analysis. This data this script produces is
consumed by generate_impute_results.py.
"""
import errno
from optparse import OptionParser
import os
import random
import re
import stat
import sys
import tempfile
import pdb
from numpy.random import choice
import numpy as np
import pandas as pd
# Our own libraries
from lib import analysis
from lib import common
_EPS=2e-4
LD_SUBSAMPLER=3
LD_SUBSAMPLER=3



USAGE = """
generate_vcf.py
No flags required by default.
Please refer to source code for command-line options.
"""

DEFAULT_SCRIPT_LOG_LOCATION = analysis.SAMPLED_DATA_BASE + "/logs/generate_vcfs_script.log"
parser = OptionParser(USAGE)
parser.add_option('--script_log', default=DEFAULT_SCRIPT_LOG_LOCATION, help='script log file')
parser.add_option('--vcftools', default='/cvmfs/soft.mugqic/CentOS6/software/vcftools/vcftools-0.1.14/bin/vcftools')
parser.add_option('--inds_data', default='/home/apoursha/imputation/compare_impute/PopInds')
parser.add_option('--java', default='/usr/bin/java')
parser.add_option('--split_vcf_jar', default='/home/apoursha/imputation/compare_impute/src/lib/splitVCFref.jar')
parser.add_option('--markers_affy', default='/home/apoursha/imputation/compare_impute/data/affy6_markers.txt')
parser.add_option('--markers_illumina', default='/home/apoursha/imputation/compare_impute/data/illumina3_markers.txt')
parser.add_option('--markers_omni', default='/home/apoursha/imputation/compare_impute/data/omni2_5_markers.txt')
parser.add_option('--markers_exome', default='/home/apoursha/imputation/compare_impute/data/annotatedList.txt')
#parser.add_option('--markers_affy_exome', default='/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/data/exome_chip/affy/Axiom_bbv09_content.core_50k_YRI.rfmt.txt')
parser.add_option('--simulated', action='store_true')

PANELS = {}
PANELS['omni2.5', 22] = int(2388927/60)
PANELS['affy6', 22]   = int(885501/60)
PANELS['illumina3', 22] = int(308330/60)


(options,args) = parser.parse_args()
if options.script_log is None:
    parser.error('no script log file given')
common.initialize_config(options)

# Helper class for loading individuals and sampling them
class Reservior(object):
    def __init__(self, items):
        self.items = random.sample(list(items), len(items))

    def consume(self, count):
        count = min(count, len(self.items))
        consumed = self.items[:count]
        self.items = self.items[count:]
        return consumed

    def half(self):
        return self.consume(len(self.items) / 2)

    def all(self):
        return self.consume(len(self.items))

class Individuals(object):
    def __init__(self):
        self.load()

    def load(self): #TODO fix this. This puts individuals in a samples dictionary based on their ances. and creates an array of all individuals
        filenames = os.listdir(options.inds_data)
        samples = {}
        all_inds = []
        for fn in filenames:
            pop = re.sub(r'\.inds$', '', fn)
            f = open(options.inds_data + "/" + fn, 'r')
            inds = [line.strip() for line in f]
            samples[pop] = inds
            all_inds += inds

        self.samples = samples
        self.pops = samples.keys()
        self.all_inds = all_inds

    def create_sampling(self, include_pops=None, exclude_pops=None):
        assert include_pops is None or exclude_pops is None
        if include_pops is not None:
            for pop in include_pops:
                assert pop in self.pops
        if exclude_pops is not None:
            for pop in exclude_pops:
                assert pop in self.pops

        if include_pops is None and exclude_pops is None:
            return Reservior(self.all_inds)
        pops = None
        if include_pops:
            pops = include_pops
        elif exclude_pops:
            pops = [pop for pop in self.pops if pop not in exclude_pops]
        inds = []
        for pop in pops:
            inds += self.samples[pop]
        return Reservior(inds)

def gzinput(chrom):
    return '/home/apoursha/imputation/with_22.vcf.gz'#'/srv/gs1/projects/bustamante/reference_panels/1kG_DATA/integrated/20120317_new_phase1_integrated_genotypes_version_3/orig_files/ALL.chr%d.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz' % chrom

def process_sample_size_analysis(individuals):
    for item in analysis.SampleSizeAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "samplesize-%04d-%04d-ch%02d-run-%02d" % (item.testsize(), item.refsize(), chrom, iteration)
            sampling = individuals.create_sampling()
            ref_inds = sampling.consume(item.refsize())
            test_inds = sampling.consume(item.testsize())
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)
def simulated_sample_size_analysis(individuals):
    for item in analysis.SampleSizeAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "samplesize-%04d-%04d-ch%02d-run-%02d" % (item.testsize(), item.refsize(), chrom, iteration)
            sampling = individuals.create_sampling()
            ref_inds = sampling.consume(item.refsize())
            test_inds = sampling.consume(item.testsize())
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)



def process_population_analysis(individuals):  #TODO this is unchecked and should not be used
    print("process_population_analysi has not been checked. therefore I'm stopping this")
    assert False
    for item in analysis.HoldOnePopAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "pop-excluded-%s-ch%02d-run-%02d" % (item.pop(), chrom, iteration)
            test_inds = individuals.create_sampling(include_pops = [item.pop()]).all()
            ref_inds = individuals.create_sampling(exclude_pops = [item.pop()]).consume(900)
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)

    for item in analysis.HoldNonePopAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "pop-included-%s-ch%02d-run-%02d" % (item.pop(), chrom, iteration)
            pop_inds = individuals.create_sampling(include_pops = [item.pop()])
            other_inds = individuals.create_sampling(exclude_pops = [item.pop()])
            ref_inds = (pop_inds.half() + other_inds.all())[:900]
            test_inds = pop_inds.all()  # All the remaining, which is the other half
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)

    # Americas is the smallest group with 3 pops, 60, 66, 55 individuals each.
    # We leave up to 33 in the test set and 60 + 66/2 + 55 = 148 in the ref set and apply it all over.
    # For Iberian people the set has size 14. In this case all of it becomes test case, and we take
    # 148 from the other pops.
    # so we are comparing apples to apples
    for item in analysis.RegionBiasedControlAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "pop-unbiased-%s-ch%02d-run-%02d" % (item.pop(), chrom, iteration)
            pop_inds = individuals.create_sampling(include_pops = [item.pop()])
            other_inds = individuals.create_sampling(exclude_pops = [item.pop()])
            test_inds = pop_inds.consume(30)
            ref_inds = random.sample(pop_inds.all() + other_inds.all(), 148)
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)

    for item in analysis.RegionBiasedAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "pop-biased-%s-ch%02d-run-%02d" % (item.pop(), chrom, iteration)
            pop_inds = individuals.create_sampling(include_pops = [item.pop()])
            other_inds = individuals.create_sampling(include_pops = item.regional_pops())
            test_inds = pop_inds.consume(33)
            ref_inds = random.sample(pop_inds.all() + other_inds.all(), 148)
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)

def create_inds_file(filename, samples):
    f = open(filename, 'w')
    f.writelines(["%s\n" % sample for sample in samples])
    f.close()

def create_vcf(jobname, outfile, vcftools_args, wait_jobs = []):
    if os.path.exists(outfile):
        # No need to recreate
        print("Output already exists, skipping: {}".format(jobname))
        return None
    # Make a log dir if it doesn't exists
    logdir = os.path.join(os.path.dirname(outfile), 'logs')
    common.mkdirp(logdir)
    return common.qsub_shell_commands(
        jobname, "Create VCF: " + jobname,
        """
        mkdir %(jobname)s &&
        cd %(jobname)s &&
        %(vcftools)s %(vcftools_args)s --recode --stdout | gzip -c - > %(outfile)s.inprogress &&
        cd .. &&
        mv %(outfile)s.inprogress %(outfile)s &&
        rm -Rf %(jobname)s
        """ % {
            'jobname' : jobname,
            'vcftools' : options.vcftools,
            'vcftools_args' : vcftools_args,
            'outfile' : outfile
        },
        wait_jobs = wait_jobs,
        logdir = logdir
        )

def chunk_vcf(jobname, vcffile, wait_jobs=[], logdir=None, core_size=5000*1000, flanking_size=250*1000):
    if os.path.exists(vcffile + ".splitlog.csv"):
        # No need to recreate
        print("Output already exists, skipping: %s" % jobname)
        return None

    logdir = os.path.join(os.path.dirname(vcffile), 'logs')
    cmd = """
      #
      # To determine the chunk size
      # First find first and last locations in the VCF
      # Then calculate the range size and divide by desired chunk size for number of chunks
      # Use the ceil of that to recalculate actual chunk size with that many chunks
      #
      MIN_AND_MAX=$(zcat %(vcf)s | grep -v -E "^#" | awk '{print $2}' | sort -g |  awk 'NR==1; END{print}')  &&
      MAX=$(echo "$MIN_AND_MAX" | tail -n 1)  &&
      MIN=$(echo "$MIN_AND_MAX" | head -n 1)  &&
      echo $MAX &&
      echo $MIN &&
      CHUNKSIZE=$(python -c "from math import ceil; range=$MAX-$MIN; chunks=ceil(range/float(%(core_size)s)); print(int(ceil(range / chunks)))") &&
      echo "Using min:$MIN max:$MAX desired_chunk_size:%(core_size)s actual_chunk_size:$CHUNKSIZE" &&
      # To work around a bug in base-pair count parsing in the jar, add period (It normally expects Kb Mb)
      %(java)s -Xmx4g -jar %(split_vcf_jar)s --vcfref %(vcf)s --coresize ${CHUNKSIZE}.b --flankingsize %(flanking_size)s.b
    """
    return common.qsub_shell_commands(
        jobname, "Chunk VCF: " + jobname,
        cmd % {
            'jobname' : jobname,
            'java': options.java,
            'split_vcf_jar' : options.split_vcf_jar,
            'vcf' : vcffile,
            'core_size' : core_size,
            'flanking_size' : flanking_size
        },
        wait_jobs = wait_jobs,
        logdir = logdir,
        vmem = 7
        )

def create_vcfs(job_name_suffix, item, chrom, iteration, src_vcf, ref_inds, test_inds):
    job_dir = item.dirname(iteration)
    common.mkdirp(job_dir)
    # Drop in a file called job_name.cfg
    f = open(os.path.join(job_dir, "job_name.cfg"), 'w')
    f.write(job_name_suffix)
    f.close()

    ref_vcf = item.vcf_filename(analysis.T_REF, chrom, iteration)
    test_vcf = item.vcf_filename(analysis.T_TEST_TRUTH, chrom, iteration)

    state = os.path.exists(ref_vcf) + os.path.exists(test_vcf)
    if state == 1:
        if os.path.exists(ref_vcf):
            os.remove(ref_vcf)
        if os.path.exists(test_vcf):
            os.remove(test_vcf)

    ref_inds_file = None
    test_inds_file = None
    if state != 2:
        ref_inds_file = item.inds_filename(analysis.T_REF)
        test_inds_file = item.inds_filename(analysis.T_TEST_TRUTH)
        create_inds_file(ref_inds_file, ref_inds)
        create_inds_file(test_inds_file, test_inds)
    make_ref_jid = create_vcf(
        "make-ref-" + job_name_suffix,
        item.vcf_filename(analysis.T_REF, chrom, iteration),
        "--gzvcf %s --keep %s --remove-indels" % (src_vcf, ref_inds_file))

    chunk_vcf(
        "chunk-ref-" + job_name_suffix,
        item.vcf_filename(analysis.T_REF, chrom, iteration),
        wait_jobs = [make_ref_jid])
    make_test_jid = create_vcf(
        "make-test-" + job_name_suffix,
        test_vcf,
        "--gzvcf %s --keep %s --remove-indels" % (src_vcf, test_inds_file))

    if options.simulated:
        sample_simulated_test_data_vcfs(
                job_name_suffix, chrom, ref_vcf, item, iteration, test_vcf,
                make_ref_jid, make_test_jid)
    else:
        sample_real_test_data_vcfs(
                job_name_suffix, item, chrom, iteration, test_vcf, make_test_jid)

def make_simulated_panel(job_name_suffix, panel, ref_vcf, wait_jobs, size, panel_name):
    panel_directory = os.path.dirname(panel)
    panel_file_name = os.path.basename(panel)
    outfile = panel_directory + "/simulated_" + panel_file_name
    setattr(options, "markers_" + panel_name, outfile)
    if os.path.exists(outfile):
        print("Output already exists, skipping: make-fake-{}-{}".format(panel_name, job_name_suffix))
        return
    # check if allele frequencies exist, if not regenerate them
    vcf_file_name = os.path.basename(ref_vcf)
    af_file_name = panel_directory + "/" + vcf_file_name + ".frq"
    logdir = os.path.join(panel_directory, 'logs')
    common.mkdirp(logdir)
    allele_freq_jid = None
    if not os.path.exists(af_file_name):
        jobname = 'Frequency' + job_name_suffix
        allele_freq_jid = common.qsub_shell_commands(
        jobname, "Create allele frequency file: " + jobname,
        """
        %(vcftools)s --gzvcf %(infile)s --freq2 --out %(outfile)s
        """ % {
            'jobname' : job_name_suffix,
            'vcftools' : options.vcftools,
            'infile' : ref_vcf,
            'outfile' : panel_directory + "/" + vcf_file_name # vcftools adds a .frq
        },
        wait_jobs = wait_jobs,
        logdir = logdir,
        sync = True
        )
        # wait for jobs to finish.
    freq = pd.read_table(filepath_or_buffer=af_file_name, sep='\t', header=0)
    chrom = np.array(freq.index)
    pos = np.array(freq.CHROM)
    freq = np.array(freq['{FREQ}'])
    freq[freq > 0.5] = 1.0 - freq[freq > 0.5] # get minor allele frequency
    #freq[freq == 0]  = _EPS
    freq = freq / np.sum(freq)
    locations = np.random.choice(range(freq.shape[0]), size=LD_SUBSAMPLER * size, replace=False, p = freq)
    locations = locations[0::LD_SUBSAMPLER]
    locations = np.sort(locations)
    del freq
    with open(outfile, 'w') as fp:
        for location in locations:
            fp.write("""%(ch)s:%(pos)s\n""" % {
                'ch' : chrom[location],
                'pos' : pos[location]})
        fp.close()

def sample_simulated_test_data_vcfs(job_name_suffix, chrom, ref_vcf, item,
        iteration, test_vcf, make_ref_jid, make_test_jid):
    #TODO actually make this work with multiple chroms
    make_simulated_panel(job_name_suffix,options.markers_affy, ref_vcf,
        make_ref_jid, PANELS['affy6', chrom], 'affy')
    make_simulated_panel(job_name_suffix,options.markers_illumina, ref_vcf,
        make_ref_jid, PANELS['illumina3', chrom], "illumina")
    make_simulated_panel(job_name_suffix,options.markers_omni, ref_vcf,
        make_ref_jid, PANELS['omni2.5', chrom], "omni")


    create_vcf(
        "make-test-affy-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_AFFY, chrom, iteration),
        "--gzvcf %s --snps %s" % (test_vcf, options.markers_affy),
        wait_jobs = [make_test_jid])

    create_vcf(
        "make-test-illumina-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_ILLUMINA, chrom, iteration),
        "--gzvcf %s --snps %s" % (test_vcf, options.markers_illumina),
        wait_jobs = [make_test_jid])

    create_vcf(
        "make-test-omni-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_OMNI, chrom, iteration),
        "--gzvcf %s --snps %s" % (test_vcf, options.markers_omni),
        wait_jobs = [make_test_jid])




#def sample_simulated_test_data_vcfs(
#        job_name_suffix, item, chrom, iteration, test_vcf, make_test_jid):
#    down_sample_vcf(
#            "make-fake-affy-" + job_name_suffix, item.vcf_filename(
#                analysis.T_TEST_AFFY, chrom, iteration),
#            test_vcf, PANELS['affy6', chrom], wait_jobs = [make_test_jid])
#
#    down_sample_vcf(
#            "make-fake-illumina-" + job_name_suffix, item.vcf_filename(
#                analysis.T_TEST_ILLUMINA, chrom, iteration),
#            test_vcf, PANELS['illumina3', chrom], wait_jobs = [make_test_jid])
#
#    down_sample_vcf(
#            "make-fake-omni" + job_name_suffix, item.vcf_filename(
#                analysis.T_TEST_OMNI, chrom, iteration),
#            test_vcf, PANELS['omni2.5', chrom], wait_jobs = [make_test_jid])
#

def sample_real_test_data_vcfs(
        job_name_suffix, item, chrom, iteration, test_vcf, make_test_jid):
    create_vcf(
        "make-test-affy-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_AFFY, chrom, iteration),
        "--gzvcf %s --snps %s" % (test_vcf, options.markers_affy),
        wait_jobs = [make_test_jid])

    create_vcf(
        "make-test-illumina-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_ILLUMINA, chrom, iteration),
        "--gzvcf %s --snps %s" % (test_vcf, options.markers_illumina),
        wait_jobs = [make_test_jid])

    create_vcf(
        "make-test-omni-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_OMNI, chrom, iteration),
        "--gzvcf %s --snps %s" % (test_vcf, options.markers_omni),
        wait_jobs = [make_test_jid])

    create_vcf(
        "make-test-exome-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_EXOME, chrom, iteration),
        "--gzvcf %s --positions %s" % (test_vcf, options.markers_exome),
        wait_jobs = [make_test_jid])

    create_vcf(
        "make-test-affy-exome-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_AFFY_EXOME, chrom, iteration),
        "--gzvcf %s --positions %s" % (test_vcf, options.markers_affy_exome),
        wait_jobs = [make_test_jid])

def down_sample_vcf(jobname, outfile, test_vcf, size, wait_jobs = []):
    if os.path.exists(outfile):
        # No need to recreate
        print("Output already exists, skipping: {}".format(jobname))
        return None
    # Make a log dir if it doesn't exists
    logdir = os.path.join(os.path.dirname(outfile), 'logs')
    common.mkdirp(logdir)
    print(size)
    return common.qsub_shell_commands(
            jobname, "Create VCF: " + jobname,
            """
            mkdir %(jobname)s &&
            cd %(jobname)s &&
            zcat %(vcf_file)s | head -20 | grep '#' > %(outfile)s.inprogress &&
            zcat %(vcf_file)s | grep -v '#' | shuf -n %(panel_size)d | sort -g -k2 >> %(outfile)s.inprogress &&
            gzip %(outfile)s.inprogress
            cd .. &&
            mv %(outfile)s.inprogress.gz %(outfile)s &&
            rm -Rf %(jobname)s
            """ % {
            'jobname' : jobname,
            'vcf_file' : test_vcf,
            'panel_size' : int(size),
            'outfile' : outfile
        },
        wait_jobs = wait_jobs,
        logdir = logdir
        )




def main():
    individuals = Individuals()
    process_sample_size_analysis(individuals)
    #process_population_analysis(individuals)

if __name__ == "__main__":
    main()
