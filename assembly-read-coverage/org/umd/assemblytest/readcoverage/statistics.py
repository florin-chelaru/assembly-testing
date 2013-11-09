'''
Created on Nov 7, 2013

@author: jason & kostas
'''
import numpy as np
from coverage import parse_fasta_file, ContigBPCoverage, ContigWindowCoverage
from samfile import SamFile
import matplotlib.pyplot as plt

class CoverageStatistics(object):

    def __get_all_values__(self, cov):
        # Beautiful code
        allvals = [cov.contig_coverage[cid].tolist() for cid in cov.contig_coverage]
        allvals = np.array([item for sublist in allvals for item in sublist])
        return allvals

    def __map_win_to_bp__(self, bool_mask, cov, cid):
        bp = np.zeros(cov.contig_length[cid], dtype=np.bool)
        mask_ind = np.flatnonzero(bool_mask)
        for i in range(len(mask_ind)):
            bp_pos = cov.contig_window_start_index[cid][mask_ind[i]]
            bp[bp_pos:(bp_pos + cov.window_length)] = True
        return bp

    def __init__(self, cov, test_type, param):
        self.test_type = test_type
        if test_type == 'Gaussian':
            if param < 0:
                raise RuntimeError, "CoverageStatistics(): invalid standard deviation multiplier (expected non-negative value)."
            allvals = self.__get_all_values__(cov)
            self.mu = np.mean(allvals)
            self.sigma = np.std(allvals)
            self.t_high = self.mu + param * self.sigma
            self.t_low = self.mu - param * self.sigma
        elif test_type == 'Percentile':
            if param < 0:
                raise RuntimeError, "CoverageStatistics(): invalid percentile (expected float value in [0,1])."
            allvals = self.__get_all_values__(cov)
            allvals = np.sort(allvals)
            f = param / 2.0
            ind_high = np.floor((1 - f) * len(allvals))
            ind_low = np.ceil(f * len(allvals))
            self.t_high = allvals[ind_high]
            self.t_low = allvals[ind_low]
        else:
            raise RuntimeError, "CoverageStatistics(): unknown test type."
        # Test coverage values against t_high and t_low
        self.contig_coverage = cov.contig_coverage
        self.contig_length = cov.contig_length
        self.contig_overcovered_windows = {}
        self.contig_undercovered_windows = {}
        self.contig_overcovered_bps = {}
        self.contig_undercovered_bps = {}
        for cid in cov.contig_coverage:
            win_over = cov.contig_coverage[cid] > self.t_high
            win_under = cov.contig_coverage[cid] < self.t_low
            bp_over = self.__map_win_to_bp__(win_over, cov, cid)
            bp_under = self.__map_win_to_bp__(win_under, cov, cid)
            # binary vectors
            self.contig_overcovered_windows[cid] = win_over
            self.contig_undercovered_windows[cid] = win_under
            self.contig_overcovered_bps[cid] = bp_over
            self.contig_undercovered_bps[cid] = bp_under

if __name__ == '__main__':
    samfile = SamFile.read('../../../../../tutorial/read_coverage/influenza-A.sam')
    cdata = parse_fasta_file('../../../../../data/influenza-A/influenza-A.assembly.fasta')
    bp_cov = ContigBPCoverage(samfile, cdata)
    w_cov = ContigWindowCoverage(bp_cov, 200, 1)
    st = CoverageStatistics(w_cov, 'Percentile', 0.05)
    c = '8'
    # plt.plot(w_cov.contig_coverage[c])
    # plt.savefig('cov.png')
    # plt.plot(st.contig_overcovered_windows[c])
    # plt.savefig('win_over.png')
    # plt.plot(st.contig_undercovered_bps[c])
    # plt.savefig('bp_under.png')
    plt.plot(st.contig_overcovered_bps[c])
    plt.savefig('bp_over.png')
