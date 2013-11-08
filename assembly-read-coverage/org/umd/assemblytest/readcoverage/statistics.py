'''
Created on Nov 7, 2013

@author: jason & kostas
'''

import numpy
from samfile import SamFile
from coverage import ContigData, ContigCoverage
# import matplotlib.pyplot as plt

class StatisticalTest(object):

    '''
    Calculates overall statistics (of certain type) from assembly and
    applies tests on coverage functions.
    '''

    '''
    Constructor docs
    '''
    def __init__(self, samfile, contig_data, type_of_test, threshold):

        # Sanity checking first....

        if type_of_test == "Gaussian":
            self.type = type_of_test;
            cl = contig_data.contig_length;
            cov = [];
            for contig_id in cl:
                pos = 0;
                while pos < cl[contig_id]:
                    cov.append(len(samfile.coverage(contig_id, pos, 1)));
                    pos += 1;
            print '%d' % len(cov)
            mu = numpy.mean(cov);
            sigma = numpy.std(cov);
            # print 'mean = %f, sigma = %f' % (mu, sigma)
            self.thresh_low = mu - threshold * sigma;
            self.thresh_high = mu + threshold * sigma;
            # plt.plot(cov);
            # plt.show();

    def apply(self, contig_cov):
        print "brb"


if __name__ == '__main__':
    cdata = ContigData('../../../../../data/influenza-A/influenza-A.assembly.fasta');
    samfile = SamFile.read('../../../../../tutorial/read_coverage/influenza-A.sam');
    st = StatisticalTest(samfile, cdata, "Gaussian", 3);
    print 'T_low = {0}, T_high = {1}'.format(st.thresh_low, st.thresh_high);
