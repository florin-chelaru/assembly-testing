'''
Created on Nov 10, 2013

@author: kzampog
'''

import argparse, sys, os.path, subprocess
from samfile import SamFile
import coverage as cov
import statistics as stat

def parse_args():
    parser = argparse.ArgumentParser(description='Detects assembly regions of extreme coverage')
    parser.add_argument('-a', '--assembly', required=True, help='Assembly file in FASTA format', metavar='ASSEMBLY_FASTA', dest='assembly_file')
    parser.add_argument('-r', '--reads', required=True, help='Read sequences file in FASTA format', metavar='READS_FASTA', dest='reads_file')
    parser.add_argument('-w', '--window_length', default=1, type=int, help='Sliding window length (integer)', metavar='WIN_LEN')
    parser.add_argument('-s', '--window_slide_step', default=1, type=int, help='Window sliding step (integer)', metavar='WIN_STEP')
    parser.add_argument('-t', '--test_type', required=True, help='Statistical test type used to determine extreme coverage; can be either Gaussian or Percentile', metavar='TEST_TYPE')
    parser.add_argument('-p', '--test_param', type=float, required=True, help='Statistical test parameter (float)', metavar='TEST_PARAM')
    parser.add_argument('-o', '--output_dir', default='./out/', required=False, help='Output file directory', metavar='OUT_DIR')
    args = parser.parse_args()
    # File and path existence checks
    err = ''
    if not os.path.exists(args.assembly_file):
        err = 'assembly file'
    if not os.path.exists(args.reads_file):
        err = 'read sequences file'
    if not os.path.isdir(args.output_dir):
        err = 'output directory'
    if err != '':
        print 'Non-existent ' + err + '.'
        sys.exit(1)
    if len(args.output_dir) > 0 and args.output_dir[-1] != '/':
        args.output_dir += '/'
    base_name = os.path.basename(args.assembly_file)[::-1]
    base_name = os.path.splitext(base_name)[-1][::-1][0:-1]
    args.base_name = base_name
    return args

def generate_sam_file(args):
    bt2name = args.output_dir + args.base_name + '.bt2'
    samname = args.output_dir + args.base_name + '.sam'
    command1 = 'bowtie2-build ' + args.assembly_file + ' ' + bt2name
    # bowtie2-align same as bowtie2 (script)?
    command2 = 'bowtie2-align -x ' + bt2name + ' -f -U ' + args.reads_file + ' -S ' + samname
    retval = subprocess.call(command1)
    if retval != 0:
        print 'BowTie2 indexing failed'
        sys.exit(1)
    retval = subprocess.call(command2)
    if retval != 0:
        print 'BowTie2 alignment failed'
        sys.exit(1)
    return samname


if __name__ == '__main__':
    args = parse_args()
    print 'Calculating alignments (SAM file) using BowTie2...'
    sam_filename = generate_sam_file(args)
    print 'Parsing SAM file...'
    samfile = SamFile.read(sam_filename)
    print 'Reading auxiliary data from assembly FASTA file...'
    cdata = cov.parse_fasta_file(args.assembly_file)
    print 'Calculating base-pair level coverage...'
    bp_cov = cov.ContigBPCoverage(samfile, cdata)
    print 'Calculating window-averaged coverage...'
    w_cov = cov.ContigWindowCoverage(bp_cov, args.window_length, args.window_slide_step)
    print 'Performing statistical test...'
    st = stat.CoverageStatistics(w_cov, args.test_type, args.test_param)
    print 'Results:'
    print st.to_string()
    print 'Writing test results to files...'
    st.write_all_files(args.output_dir + args.base_name)
    print 'Done'

