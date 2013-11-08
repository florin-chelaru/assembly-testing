'''
Created on Nov 7, 2013

@author: jason & kostas        
'''
import numpy as np
from samfile import SamFile
import matplotlib.pyplot as plt

class ContigData(object):

    '''
    Read a Fasta assembly file and store the {contig_id, contig_length}
    pairs in a dictionary.
    '''

    def __init__(self, filename):
        '''
        Constructor takes a FASTA filename and stores
        the data in the dictionary.
        '''
        self.contig_length = dict();
        self._fastafile = filename;

        # parse the FASTA file
        try :
            with open(filename) as f:
                current_key = None
                for line in f:
                    line = line.strip();
                    if line[0] == '>':
                        current_key = line[1:]  # update current key
                        self.contig_length[current_key] = 0
                    else:
                        self.contig_length[current_key] += (len(line))  # don't count '\n'
        except EnvironmentError as err:
                print "Unable to open FASTA file: {}".format(err); 


class ContigCoverage(object):

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
    def __init__(self, samfile, contig_data, window_length, step_size):

        # Sanity checking -- NEEDS FIX!
        
        if window_length < 0 or step_size < 0:
            raise RuntimeError, "ContigCoverage constructor: Provided negative argument."
        if samfile is None or contig_data is None:
            raise RuntimeError, "ContigCoverage constructor: Provided None argument."
        
        # Main algorithm
        self.window_length = window_length
        self.step_size = step_size
        self.contig_length = contig_data.contig_length
        self.contig_coverage = {};  # associates contigs with their windows' coverages.
        self.contig_window_starting_points = {}
        self.__calcBPCoverage__(samfile, contig_data)
        cl = contig_data.contig_length
        window = np.ones(window_length) / window_length
        for contig_id in cl:
            if window_length == 1:
                cov = self._bplevelcoverage_[contig_id]
            else:
                cov=np.convolve(self._bplevelcoverage_[contig_id], window, 'valid')
            starting_points = range(len(cov))
            if step_size > 1: # convolve for speed
                cov = cov[0::step_size]
                starting_points = starting_points[0::step_size]
            self.contig_coverage[contig_id] = cov
            self.contig_window_starting_points[contig_id] = starting_points
                

    def __calcBPCoverage__(self, samfile, contig_data):
        self._bplevelcoverage_ = dict()
        cl = contig_data.contig_length
        for contig_id in cl:
            cov = np.array([len(samfile.coverage(contig_id, pos, 1)) for pos in range(cl[contig_id])], dtype=np.float32)
            self._bplevelcoverage_[contig_id] = cov
        
    
    def write_coverage_plot(self, contig_id, filename):
        plt.plot(coverage.contig_coverage[contig_id]);
        plt.xlabel('Window index');
        plt.ylabel('Coverage');
        plt.savefig(filename);
        plt.show();

if __name__ == '__main__':
    cdata = ContigData('../../../../../data/influenza-A/influenza-A.assembly.fasta');
    samfile = SamFile.read('../../../../../tutorial/read_coverage/influenza-A.sam');
    coverage = ContigCoverage(samfile, cdata, 100, 1);
    print "asdfasdfasdf"
    for contig in coverage.contig_coverage:
        print "Contig with id {0} has {1} base-pairs and {2} windows.".format(contig, cdata.contig_length[contig], len(coverage.contig_coverage[contig]));
        print "Contig with id {0} has a window coverage of {1}.".format(contig, str(coverage.contig_coverage[contig]));

    print 'asdfasdf {0} asdf'.format(len(samfile.coverage('1', 3, 3)));
    #coverage.write_coverage_plot('1', 'dok.png');
