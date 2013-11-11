'''
Created on Nov 7, 2013

@author: jason & kostas
'''

import numpy as np
import matplotlib.pyplot as plt
# from samfile import SamFile

'''
Read a FASTA assembly file and store the {contig_id, contig_length}
pairs in a dictionary.
'''

def parse_fasta_file(filename):
    contig_length = dict()
    try:
        with open(filename) as f:
            current_key = None
            for line in f:
                line = line.strip();
                if line[0] == '>':
                    current_key = line[1:]  # update current key
                    contig_length[current_key] = 0
                else:
                    contig_length[current_key] += (len(line))
    except EnvironmentError as err:
        print "Unable to open FASTA file: {}".format(err);
    return contig_length

def write_coverage_plot(cov, contig_id, filename):
    plt.plot(cov.contig_coverage[contig_id]);
    plt.title('Contig {0} coverage'.format(contig_id))
    plt.xlabel('Window index');
    plt.ylabel('Coverage');
    plt.savefig(filename);
    # plt.show();

class ContigBPCoverage:

    def __init__(self, samfile, contig_length):
        self.contig_length = contig_length
        self.contig_coverage = dict()
        for contig_id in contig_length:
            cov = np.array([len(samfile.coverage(contig_id, pos, 1)) for pos in range(contig_length[contig_id])], dtype=np.float32)
            self.contig_coverage[contig_id] = cov

class ContigWindowCoverage(object):

    '''
    Computes the coverage of a contig given a
    specified window size and sliding step for the window.
    Note that binning is a special case of sliding windows when the step size
    just happens to be equal to window length.
    '''

    '''
    Constructor
    @param samfile a Samfile object
    @param contig_data a ContigData object
    @param window_length a positive integer representing the window length
    @param step_size a positive integer representing the step size for the window.
    '''
    def __init__(self, bp_coverage, window_length, step_size):
        # Sanity checking
        if window_length < 0 or step_size < 0:
            raise RuntimeError, "ContigCoverage constructor: Provided negative argument."
        if bp_coverage is None:
            raise RuntimeError, "ContigCoverage constructor: Provided None argument."
        # Main algorithm
        self.window_length = window_length
        self.step_size = step_size
        self.contig_length = bp_coverage.contig_length
        self.contig_coverage = {};  # associates contigs with their windows' coverages
        self.contig_window_start_index = {}
        window = np.ones(window_length) / window_length
        for contig_id in bp_coverage.contig_length:
            if window_length == 1:
                cov = bp_coverage.contig_coverage[contig_id]
            else:
                cov = np.convolve(bp_coverage.contig_coverage[contig_id], window, 'valid')
            starting_points = range(len(cov))
            if step_size > 1:
                cov = cov[0::step_size]
                starting_points = starting_points[0::step_size]
            self.contig_coverage[contig_id] = cov
            self.contig_window_start_index[contig_id] = starting_points

'''
if __name__ == '__main__':
    samfile = SamFile.read('../../../../../tutorial/read_coverage/influenza-A.sam')
    cdata = parse_fasta_file('../../../../../data/influenza-A/influenza-A.assembly.fasta')
    bp_cov = ContigBPCoverage(samfile, cdata)
    w_cov = ContigWindowCoverage(bp_cov, 100, 1)
    # write_coverage_plot(w_cov, '8', 'dok.png')
    print len(w_cov.contig_coverage)
    cov = w_cov
    # write_coverage_plot(bp_cov, '1', '')
    for c in cov.contig_coverage:
        print "Contig {0} has {1} bp's, {2} windows and coverage: {3}".format(c, cov.contig_length[c], len(cov.contig_coverage[c]), cov.contig_coverage[c]);
'''
