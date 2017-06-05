"""
phase_impute_pipeline.py
author: Alicia Martin (armartin@stanford.edu)
author: Gerard Tse (gerardtse@gmail.com)

Schedule one phase and impute pipeline. Four pipelines are implemented in this script,
and can be selected by changing the --step parameter. These pipelines are
--step=1 : Beagle
--step=2 : MaCH-Admix
--step=3 : ShapeIt / Impute2
--step=4 : SLRP / Impute2
"""

from optparse import OptionParser
import os
import stat
import subprocess
import sys
import tempfile
import numpy as np
import re
import pdb

from lib import common
from lib import split_vcf_info
HOME='/sb/project/ams-754-aa/apoursh_projects/imputation/compare_impute'

USAGE = """
phase_impute_pipeline.py    --ind_ref <file containing reference individuals>
                            --ind_count <number of individuals in reference plus test panel>
                            --ref_gzvcf <input reference panel vcf file for phasing>
                            --test_gzvcf <input test set vcf file for phasing>
                            --step <step to perform>
                            --logs <directory for cluster log files>
                            --script_dir <script dir file>
                            --script_log <script log file>
                            --job_name_prefix <prefix for all job names>
                            --intermediate_dir <directory for intermediate files>
                            --completed_file <a file to be written when processing is complete>
                            --outdir <output directory>
                            --chr <chromosome>
                            --qsub_id_file <file to store created job ids>
            (OPTIONAL):
                            --recomb_map <genetic recombination map>
                            --java <java path>
                            --gatk <gatk path>
                            --beagle <beagle path>
                            --mach <mach path>
                            --shapeit <shapeit path>
                            --impute2 <impute2 path>
                            --slrp <slrp path>
                            --vcftools <vcftools path>
                            --plink <plink path>
                            --pseq <pseq path>
                            --perl <perl path>
                            --phased

NOTE: steps are as follows:
1   Beagle phasing and imputation
2   MaCH phasing and imputation
3   SHAPEIT2 phasing followed by IMPUTE2 imputation
4   SLRP phasing followed by IMPUTE2 imputation  Currently not working
"""
parser = OptionParser(USAGE)
parser.add_option('--ind_ref', help='reference individuals')
parser.add_option('--ind_count', help='number of individuals in ref and test panel total')
parser.add_option('--ref_gzvcf', help='Zipped VCF file of reference panel')
parser.add_option('--test_gzvcf', help='Zipped VCF file of test set')
parser.add_option('--step', help='which steps to perform') # 1, 2, 3, or 4
parser.add_option('--logs', help='directory for cluster log files')
parser.add_option('--script_log', help='script log file')
parser.add_option('--script_dir', help='script directory file')
parser.add_option('--job_name_prefix', help='prefix for all job names')
parser.add_option('--intermediate_dir', help='directory for intermediate files')
parser.add_option('--outdir', help='output directory')
parser.add_option('--completed_file', help='A file to be written when the processing is complete')
parser.add_option('--chr', help='chromosome')
parser.add_option('--qsub_id_file', help='file for job ids')
parser.add_option('--qsub_use_testq', default='0', help='use the test queue for scheduling')
parser.add_option('--qsub_run_locally', default='0', help='run jobs locally instead of using qsub')
parser.add_option('--recomb_map', help='genetic recombination map', default=HOME + '/../recombinationMaps/genetic_map_GRCh37_chr22.txt')
parser.add_option('--ref', default='/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa')
parser.add_option('--java', default='/cvmfs/soft.mugqic/CentOS6/software/java/openjdk-jdk1.8.0_72/bin/java')
parser.add_option('--gatk', default='/cvmfs/soft.mugqic/CentOS6/software/GenomeAnalysisTK/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar')
#parser.add_option('--beagle', default='/cvmfs/soft.mugqic/CentOS6/software/beagle/beagle-03May16.862/beagle.03May16.862.jar')
#parser.add_option('--beagle', default='/cvmfs/soft.mugqic/CentOS6/software/beagle/beagle-r1399/beagle.r1399.jar')
parser.add_option('--beagle', default=HOME + '/src/bin/beagle.jar')
parser.add_option('--vcf2beagle', default=HOME + '/src/vcf2beagle.jar')
parser.add_option('--mach', default=HOME + '/src/bin/mach-admix')
parser.add_option('--shapeit', default=HOME + '/src/bin/shapeit')
parser.add_option('--impute2', default=HOME + '/src/bin/impute2')
#parser.add_option('--slrp', default='/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/tools/bin/SLRP')  #this one is hard to install
parser.add_option('--vcftools', default='/cvmfs/soft.mugqic/CentOS6/software/vcftools/vcftools-0.1.14/bin/vcftools')

parser.add_option('--plink', default=HOME + '/src/bin/plink')
parser.add_option('--pseq', default=HOME + '/src/bin/pseq')
parser.add_option('--perl', default='/usr/bin/perl')
parser.add_option('--phased', action='store_true')
parser.add_option('--first_last', default=HOME+'/data/hg19_chromInfo_view.csv')

(options,args) = parser.parse_args()
#perform cursory checks that we have all required arguments
if options.ind_ref is None:
    parser.error('ind_ref is not given')
if options.ind_count is None:
    parser.error('ind_count is not given')
if options.ref_gzvcf is None:
    parser.error('ref_gzvcf path not given')
if options.test_gzvcf is None:
    parser.error('test_gzvcf path not given')
if options.step is None:
    parser.error('no steps to run')
if options.logs is None:
    parser.error('no logs directory given')
if options.script_log is None:
    parser.error('no script log file given')
if options.script_dir is None:
    parser.error('no script_dir file given')
if options.job_name_prefix is None:
    parser.error('no job name prefix given')
if options.intermediate_dir is None:
    parser.error('no intermediate directory given')
if options.completed_file is None:
    parser.error('no completed_file is given')
if options.outdir is None:
    parser.error('no output directory given')
if options.chr is None:
    parser.error('no chr given')
if options.qsub_id_file is None:
    parser.error('no qsub_id_file given')
common.initialize_config(options)

def jobname(name):
    return "%s-step%s-%s" % (options.job_name_prefix, options.step, name)

# Short hand for assigning vmem value to jobs.
def vmem(small, medium, large):
    if options.ind_count < 100:
        return small if small > 3 else None
    elif options.ind_count < 500:
        return medium if medium > 3 else None
    return large if large > 3 else None


########################Step 1) Perform population phasing and imputation with Beagle######################
def beagle():

    prefix = os.path.join(options.intermediate_dir, '')
    test_beagle = prefix + "_unphased_test"
    ref_beagle_pattern = prefix + "_phased_ref.region%02d"
    ref_marker_pattern = ref_beagle_pattern + ".markers"
    beagle_combined_out = prefix + "_impute_test"
    beagle_out_pattern = prefix + "_out%02d"
    beagle_out_combined = prefix + "_out"

    split_info = split_vcf_info.SplitVcfInfo(options.ref_gzvcf)
    assert split_info.valid()
    convert_cmd = """
      VCF=$(mktemp %(vcf)s.unzip.XXXXXXXX)
      zcat %(vcf)s > $VCF
      trap "rm $VCF" EXIT   # Make sure it is cleaned up

      #make ref.bgl.gz file
      cat $VCF | %(java)s -Xmx4g -jar %(vcf2beagle_jar)s %(beagle_prefix)s_missing_ref %(beagle_prefix)s

      #now make markers file for beagle run (required since we have two inputs)!
      #need to only print SNPs where there is one ref and one alt allele
      if [ %(marker_file)s != "/dev/null" ]; then
        grep -v '^#' $VCF | awk 'length($4)==1 && length($5)==1 {print $3 "\t" $2 "\t" $4 "\t" $5}' > %(marker_file)s
      fi

      trap - EXIT
      rm $VCF"""

    convert_test_jid = common.qsub_shell_commands(
        jobname("beagle-convert-test"),
        "Step 1b) Convert input files to beagle (test)",
        convert_cmd % {
            'vcf': options.test_gzvcf,
            'java': options.java,
            'vcf2beagle_jar': options.vcf2beagle,
            'beagle_prefix': test_beagle,
            'marker_file': '/dev/null',
            },
        vmem = 7)

    process_region_jobs = []
    for i in range(split_info.num_regions()):
        convert_region_jid = common.qsub_shell_commands(
            jobname("beagle-convert-region-%d" % i),
            "Step 1b) Convert input files to beagle (region %d)" % i,
            convert_cmd % {
                'vcf': split_info.filename(i),
                'java': options.java,
                'vcf2beagle_jar': options.vcf2beagle,
                'beagle_prefix': ref_beagle_pattern % i,
                'marker_file': ref_marker_pattern % i,
            },
            vmem = 7)
#changed beagle prefix from PHASED_OUT_PREFIX=%(beagle_out_prefix)s.$(basename %(phased)s)
        if options.phased:
            beagle_cmd = """{java} -Xmx8g -jar {beagle_jar} phased={to_impute} missing=? phased={phased} markers={markers} out={beagle_out_prefix} lowmem=true
            """.format(
                    java      = options.java,
                    beagle_jar= options.beagle,
                    to_impute = test_beagle + ".bgl.gz",
                    phased    = (ref_beagle_pattern % i) + ".bgl.gz",
                    markers   = ref_marker_pattern % i,
                    beagle_out_prefix= beagle_out_pattern % i)
        else:
            beagle_cmd = """%(java)s -Xmx8g -jar %(beagle_jar)s unphased=%(unphased)s missing=? phased=%(phased)s markers=%(markers)s out=%(beagle_out_prefix)s lowmem=true
            """.format(
                java= options.java,
                beagle_jar= options.beagle,
                unphased= test_beagle + ".bgl.gz",
                phased= (ref_beagle_pattern % i) + ".bgl.gz",
                markers= ref_marker_pattern % i,
                beagle_out_prefix= beagle_out_pattern % i)

        beagle_cmd += """
          PHASED_OUT_PREFIX=%(beagle_out_prefix)s.$(basename %(unphased)s)
          # Add position data from marker file to second column of results and filter out markers not in core range.
          paste <(cut -f 2 %(markers)s) ${PHASED_OUT_PREFIX}.r2 |
              awk '$1 >= %(core_starts)s && $1 <= %(core_ends)s' | awk '{t=$1;$1=$2;$2=t;print}' > %(out_prefix)s.r2
          zcat ${PHASED_OUT_PREFIX}.phased.gz | head -n 1 | sed 's/ / pos /' | gzip > %(out_prefix)s.phased.gz
          paste -d ' ' <(cut -f 2 %(markers)s) <(zcat ${PHASED_OUT_PREFIX}.phased.gz | tail -n +2) |
              awk '$1 >= %(core_starts)s && $1 <= %(core_ends)s' |
              awk '{new1=$2; new2=$3; new3=$1; $1=new1; $2=new2; $3=new3; print}' | gzip >> %(out_prefix)s.phased.gz
        """
        beagle_jid = common.qsub_shell_commands(
            jobname("beagle-impute-region-%d" % i),
            "Step 1c) Perform phasing and imputation with Beagle",
            beagle_cmd % {
                'unphased': test_beagle + ".bgl.gz",
                'markers': ref_marker_pattern % i,
                'beagle_out_prefix': beagle_out_pattern % i,
                'core_starts': split_info.core_boundaries(i)[0],
                'core_ends': split_info.core_boundaries(i)[1],
                'out_prefix': beagle_out_pattern % i,   # Reusing it is fine. Beagle put a whole lot of suffix.
            },
            vmem = 12,
            wait_jobs = [ convert_region_jid, convert_test_jid ],
            queue = "extended")
        process_region_jobs.append(beagle_jid)

    region_phased_gz = [ (beagle_out_pattern % i) + ".phased.gz" for i in range(split_info.num_regions()) ]
    region_r2 = [ (beagle_out_pattern % i) + ".r2" for i in range(split_info.num_regions()) ]
    combined_phased_gz = beagle_out_combined + ".phased.gz"
    combined_r2 = beagle_out_combined + ".r2"
    cmd = """
      zcat %(first_phased_gz)s | head -n 1 | gzip - > %(combined_phased_gz)s
      for phased_gz in %(all_phased_gz)s; do zcat ${phased_gz} | tail -n +2 | gzip - >> %(combined_phased_gz)s; done
      cat %(all_r2)s >> %(combined_r2)s
    """
    combine_results = common.qsub_shell_commands(
        jobname("beagle-combine-results"),
        "Step 2d) Combine beagle results",
        cmd % {
            'first_phased_gz': region_phased_gz[0],
            'all_phased_gz': " ".join(region_phased_gz),
            'combined_phased_gz': combined_phased_gz,
            'all_r2': " ".join(region_r2),
            'combined_r2': combined_r2,
        },
        wait_jobs = process_region_jobs)

    return {
        'jobs': [ combine_results ],
        'files': [ combined_phased_gz, combined_r2 ]
    }

########################Step 3) Perform population phasing and imputation with Shapeit+Impute2######################
def shapeit_impute():
    #TODO the string formatting here is pretty horrible.
    #TODO these should all share their ref phasing
    #same for run_impute2

    #turn test vcfs into plink files for shapeit
    prefix = os.path.join(options.intermediate_dir, '')
    plink_test = prefix + '_test'
    leghap_prefix = prefix + '_ref.haplotypes'
    shapeit= prefix + '_shapeit_test'
    recomb_map = re.sub("_chr(\d+)_", '_chr' + options.chr + '_', options.recomb_map)
    test_impute2 = '_impute2_test'

    if  not options.phased:
        #generate shapeit reference panel
        cmd = ('%s %s -vcf %s -leghap %s -chr %s'
               % (options.perl, options.script_dir + '/lib/vcf2impute_legend_haps.pl', options.ref_gzvcf, leghap_prefix, options.chr))
        job2 = common.qsub(
            jobname("make_shapeit_ref"),
            'Step 3b) Perform population phasing and imputation with Shapeit+Impute2', cmd)
        # does this need qsubbing?
        #make minimal .sample file to use reference panel in phasing
        cmd = ('sh %s/lib/make_sample.sh %s %s'
               % (options.script_dir, leghap_prefix + '.sample', options.ind_ref))
        job3 = common.qsub(
            jobname('make_sample'),
            'Step 3c) Perform population phasing and imputation with Shapeit+Impute2', cmd)

        #run shapeit
        cmd=('%s --input-vcf %s --input-map %s --input-ref %s.hap.gz %s.legend.gz %s.sample --output-max %s %s'
                % (options.shapeit, options.test_gzvcf, recomb_map, leghap_prefix, leghap_prefix, leghap_prefix, shapeit + '.hap', shapeit + '.sample'))
        job4 = common.qsub(
            jobname('shapeit2'),
            'Step 3d) Perform population phasing and imputation with Shapeit+Impute2', cmd,
            wait_jobs = [job2, job3])

    else:
        ref_vcf = options.ref_gzvcf
        test_vcf = options.test_gzvcf
        job4 = extract_phase_from_vcf(test_vcf, ref_vcf, shapeit+'.hap', leghap_prefix)

    #chunk and submit impute2 jobs
    job_ids = execute_impute2(test_hap=shapeit+'.hap', ref_hap=leghap_prefix+'.hap.gz',
            ref_legend=leghap_prefix+'.legend.gz', recomb_map=recomb_map,
            output=prefix+test_impute2, chrom=options.chr, first_last=options.first_last, jobs=job4,
            wd=prefix)

    f = open(common.CONFIG.QSUB_ID_FILE, 'w')
    for jid in job_ids:
        f.write("%s\n" % jid)
    f.close()
    return {
        'jobs': job_ids,
        'files': [ os.path.join(prefix, test_impute2 + ".gz"), os.path.join(prefix, test_impute2 + '.gz.info') ]
    }

def extract_phase_from_vcf(test_vcfName, ref_vcfName, test_outname,
        ref_outname,format='shapeit2'):
    """Extracts the phasing and writes files in a particular format with outname as the name of the file. The function returns the job_id for the submitted job. Takes gzipped vcfs, returns gzipped ref file.
    """
    if format == 'shapeit2':
        if not os.path.exists(test_outname):
            cmd = """zcat {gzvcf} | grep -v '^#' | cut -f 6-9 --complement | awk ' {{ t = $2; $2 = $3; $3 = t; print; }} ' | sed 's/|\|\t/ /g' > {outfile}"""
            jid1 = common.qsub_shell_commands(
            jobname("extract-phase-info"),
            "Step 3b) Extract phasing info from test vcf: ",
            cmd.format(
                    gzvcf = test_vcfName,
                    outfile    = test_outname)
            )
        #check if the reference phasing exists
        if not os.path.exists(ref_outname + ".hap.gz"):
            cmd = """zcat {gzvcf} | grep -v '^#' | cut -f 1-9 --complement | sed 's/|\|\t/ /g' | gzip -c > {outfile}.gz """

            legend_cmd = """echo -e id \t position\t a0\t a1 > {outname} &&
            zcat {gzvcf} | grep -v '^#' | awk '{{print $3, $2, $4, $5}}' >> {outname} &&
            gzip {outname}"""
            cmd += '&&\n' + legend_cmd
            jid2 = common.qsub_shell_commands(
            jobname("extract-phase-info"),
            "Step 3c) Extract phasing info from test vcf: ",
            cmd.format(
                    gzvcf      = ref_vcfName,
                    outfile    = ref_outname + '.hap',
                    outname    = ref_outname + '.legend')
            )

        return(jid1,jid2)






########################Step 2) Perform population phasing and imputation with MaCH######################
def mach():
    # Intermediate files
    prefix = os.path.join(options.intermediate_dir, '')
    test_ped_prefix = prefix + "_test_ped"
    test_ped = test_ped_prefix + ".ped"
    test_map = test_ped_prefix + ".map"
    test_merlin_ped = prefix + "_test_merlin.ped"
    test_merlin_dat = prefix + "_test_merlin.dat"

    mach_out_combined = prefix + "_mach_imputed"
    mach_out_pattern = mach_out_combined + ".region%02d"

    split_info = split_vcf_info.SplitVcfInfo(options.ref_gzvcf)
    assert split_info.valid()

    cmd = """
       %(vcftools)s --gzvcf %(gzvcf)s --plink --out %(ped_prefix)s &&
       cut -f1-5,7- %(ped)s > %(ped_out)s && awk '{print \"M\", $2}' %(map)s > %(dat_out)s
    """
    make_test_merlin = common.qsub_shell_commands(
        jobname("make_test_merlin"),
        "Step 2a) Create Merlin input files for reference set",
        cmd % {
            'vcftools' : options.vcftools, 'gzvcf' : options.test_gzvcf, 'ped_prefix' : test_ped_prefix,
            'ped': test_ped, 'ped_out': test_merlin_ped, 'map': test_map, 'dat_out': test_merlin_dat})

    cmd = "%(mach)s --geno --compact --forceImputation --vcfRef --outvcf --datfile %(test_dat)s --pedfile %(test_ped)s --haps %(ref_vcf)s --outputstart %(core_starts)s --outputend %(core_ends)s --prefix %(out_prefix)s"
    if not options.phased:
        cmd += " --phase"
    process_region_jobs = []
    for i in range(split_info.num_regions()):
        jobid = common.qsub(
            jobname("mach-region-%d" % i),
            "Step 2b) MaCH: Phase and Impute region %d" % i,
            cmd % {
                'mach': options.mach,
                'test_dat': test_merlin_dat,
                'test_ped': test_merlin_ped,
                'ref_vcf': split_info.filename(i),
                'core_starts': split_info.core_boundaries(i)[0],
                'core_ends': split_info.core_boundaries(i)[1],
                'out_prefix': mach_out_pattern % i
            },
            wait_jobs = [ make_test_merlin ])
        process_region_jobs.append(jobid)

    region_output_vcfs = [ (mach_out_pattern % i) + ".vcf.gz" for i in range(split_info.num_regions()) ]
    region_output_infos = [ (mach_out_pattern % i) + ".info" for i in range(split_info.num_regions()) ]
    combined_vcf = mach_out_combined + ".vcf.gz"
    combined_info = mach_out_combined + ".info"
    cmd = """
      zcat %(first_vcf_gz)s | grep -E "^#" | gzip - > %(combined_vcf_gz)s  &&
      zcat %(all_vcf_gz)s | grep -v -E "^#" | gzip - >> %(combined_vcf_gz)s  &&
      head -n 1 %(first_info)s  > %(combined_info)s   &&
      tail -n +2 %(all_info)s | sed '/^$/d' | grep -v == >> %(combined_info)s
    """
    combine_results = common.qsub_shell_commands(
        jobname("mach-combine-results"),
        "Step 2c) Combine MaCH results",
        cmd % {
            'first_vcf_gz': region_output_vcfs[0],
            'all_vcf_gz': " ".join(region_output_vcfs),
            'combined_vcf_gz': combined_vcf,
            'first_info': region_output_infos[0],
            'all_info': " ".join(region_output_infos),
            'combined_info': combined_info,
        },
        wait_jobs = process_region_jobs)

    return { 'jobs' : [combine_results],
             'files' : [ combined_vcf, combined_info ] }

########################Step 4) Perform population phasing and imputation with SLRP+Impute2######################
def slrp_impute():
    pass   #fuck slrp
    # Intermediate files
    prefix = options.outdir + options.outbase + "_slrp"
    test_vcf = prefix + "_test.recode.vcf"
    ref_vcf = prefix + "_ref.recode.vcf"
    slrp_out = prefix + "_slrp"

    # Doesn't quite work. Suffers from MemoryError
    cmd = "%(slrp)s --vcf %(vcf)s --fastPreProc --outVCF %(out)s --verbose"
    make_ref_ped = common.qsub(
        jobname("slrp"),
        "Step 4a) SLRP Phasing",
        cmd % {'slrp' : options.slrp, 'vcf' : ref_vcf, 'out' : slrp_out},
        is_script = True)

def execute_impute2(test_hap, ref_hap, ref_legend, recomb_map,
        output, jobs, chrom, first_last, wd, int_size=5e6):
    chrom_dict = {}
    with open(first_last, 'r') as fl:
        first_last = fl.readlines()
        for line in first_last:
            line = line.strip().split(',')
            chrom_dict[line[0].lstrip('chr')] = line[1]
    pos1 = 1
    last_pos     = int(chrom_dict[chrom])
    grid_spacing = (last_pos - pos1) / float(int_size)
    intervals    = np.linspace(pos1 - 1, last_pos, grid_spacing)
    all_jobs, all_iters, all_info = [], [], []
    #pdb.set_trace()
    for i in range(len(intervals)-1):
        outfile = output + '_iter' + str(i)
        cmd = '{impute2} -known_haps_g {test_haps} -m {map} -int {start} {end} -h {ref_hap} -l {legend} -o {outfile}'.format(
            impute2  = options.impute2,
            test_haps= test_hap,
            map      = recomb_map,
            start    = int(intervals[i] + 1),
            end      = int(intervals[i + 1]),
            ref_hap  = ref_hap,
            legend   = ref_legend,
            outfile  = outfile
        )
        job_id = common.qsub("impute2_{}".format(i),
            "submitting region {}".format(i),
            cmd,
            wait_jobs = jobs,
            logdir = wd
        )
        all_jobs.append(job_id)
        all_iters.append(outfile)
        all_info.append(outfile + '_info')

    cmd = 'cat {} | gzip > {}.gz'.format(
            ' '.join(all_iters),
            output
            )
    job_id_impute = common.qsub_shell_commands(
            "cat-job",
            "Combining various impute2 files:\n",
            cmd,
            wait_jobs = all_jobs,
            logdir = wd
            )
    # combine info files
    cmd = """tail -n +2 {info} | sed '/^$/d' | grep -v == > {out}.gz.info""".format(
            out = output,
            info = ' '.join(all_info)
            )
    job_id_info = common.qsub_shell_commands(
            "info-concat",
            "Combining various impute2 info files:\n ",
            cmd,
            wait_jobs = all_jobs,
            logdir = wd
            )
    return job_id_impute, job_id_info


########################Run pipeline steps########################
def main():
    if options.step == '1':
        output = beagle()
    elif options.step == '2':
        output = mach()
    elif options.step == '3':
        output = shapeit_impute()
    elif options.step == '4':
        job = slrp_impute()
    common.qsub_shell_commands(
        jobname("finalize"),
        "Copy results and cleanup",
        """
        OUTPUT_FILES="$(ls %(outputs)s)"

        if [ $(echo $OUTPUT_FILES | wc -w) -lt %(expected_output_count)d ]; then
           echo Missing output files. Pipeline failed >&2
           exit 125
        fi

        for f in $OUTPUT_FILES; do
            size=$(du -b $f | cut -f 1)
            if [ "$size" == "0" ]; then
              echo Required output file is empty: $f. Pipeline failed >&2
              exit 125
            fi
        done

        # Copy results out
        for f in $OUTPUT_FILES; do cp $f %(outdir)s/$(basename $f); done
        # Delete inetermediate
        # rm -Rf %(intermediate_dir)s
        # Set marker file
        date > %(completed_file)s
        """ % {
            'outputs': " ".join(output['files']),
            'expected_output_count': len(output['files']),
            'outdir': options.outdir,
            'intermediate_dir': options.intermediate_dir,
            'completed_file': options.completed_file
        },
        wait_jobs = output['jobs'])

if __name__ == "__main__":
    main()
