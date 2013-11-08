'''
Created on Nov 7, 2013

@author: jason & kostas
'''
import inspect
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
        '''
        if window_length < 0 or step_size < 0:
            raise RuntimeError, "{0}: Provided negative argument.".format(inspect.stack()[0][3])
        if samfile is None or contig_data is None:
            raise RuntimeError, "{0}: Provided None argument.".format(inspect.stack()[0][3])
        '''

        # Main algorithm
        cl = contig_data.contig_length;
        self.window_length = window_length;
        self.step_size = step_size;
        self.contig_coverage = {};  # associates contigs with their windows' coverages.
        self.contig_window_starting_points = {}
        for contig_id in cl:
            window_coverage = [];
            window_starting_points = [];
            current_pos = 1;
            while current_pos + window_length - 1 <= cl[contig_id]:
                window_starting_points.append(current_pos);
                curr_window_coverage = self.__calcwindowcoverage__(samfile, contig_id, current_pos, window_length)
                window_coverage.append(curr_window_coverage);
                current_pos += step_size;
            # Do I have anything remaining in the contig?
            if current_pos < cl[contig_id]:
                window_starting_points.append(current_pos);
                last_window_coverage = self.__calcwindowcoverage__(samfile, contig_id, current_pos, cl[contig_id] - current_pos + 1)
                window_coverage.append(last_window_coverage);

            self.contig_coverage[contig_id] = window_coverage;
            self.contig_window_starting_points[contig_id] = window_starting_points;

    # Calculate the coverage over an entire window
    def __calcwindowcoverage__(self, samfile, contig_id, current_pos, window_length):
        alignments = samfile.coverage(contig_id, current_pos, window_length);
        total_window_coverage = 0;
        for aln in alignments:
            starting_overlap_index = max(aln.start(), current_pos);
            ending_overlap_index = min(aln.end(), current_pos + window_length - 1);
            total_window_coverage += (1.0 * (ending_overlap_index - starting_overlap_index + 1)) / window_length;
        return total_window_coverage;

    def write_coverage_plot(self, contig_id, filename):
        plt.plot(coverage.contig_coverage[contig_id]);
        plt.xlabel('Window index');
        plt.ylabel('Coverage');
        plt.savefig(filename);
        # plt.show();

if __name__ == '__main__':
    cdata = ContigData('../../../../../data/influenza-A/influenza-A.assembly.fasta');
    samfile = SamFile.read('../../../../../tutorial/read_coverage/influenza-A.sam');
    coverage = ContigCoverage(samfile, cdata, 1, 1);
    print "asdfasdfasdf"
    for contig in coverage.contig_coverage:
        print "Contig with id {0} has {1} base-pairs and {2} windows.".format(contig, cdata.contig_length[contig], len(coverage.contig_coverage[contig]));
        print "Contig with id {0} has a window coverage of {1}.".format(contig, str(coverage.contig_coverage[contig]));

    # print 'asdfasdf {0} asdf'.format(len(samfile.coverage('1', 3, 3)));
    # coverage.write_coverage_plot('1', 'dok.png');
