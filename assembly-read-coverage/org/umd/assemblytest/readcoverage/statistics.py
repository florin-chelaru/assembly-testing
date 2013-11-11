'''
Created on Nov 7, 2013

@author: jason & kostas
'''

import numpy as np
# from coverage import parse_fasta_file, ContigBPCoverage, ContigWindowCoverage
# from samfile import SamFile
# import matplotlib.pyplot as plt

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

    def __map_bp_to_int__(self, bp):
        starts = []
        ends = []
        i = 1
        if bp[0] == True:
            starts = [0]
        while i < len(bp):
            if bp[i] == True and bp[i - 1] == False:
                starts += [i]
            if bp[i] == False and bp[i - 1] == True:
                ends += [i - 1]
            i += 1
        if bp[-1] == True:
            ends = ends + [len(bp) - 1]
        return zip(starts, ends)

    def __tuple_list_to_string__(self, template, tl):
        retval = ''
        for t in tl:
            retval += template.format(*t)
        return retval

    def write_bv_to_file(self, bvm, file_name):
        try :
            with open(file_name, 'w') as f:
                for cid in bvm:
                    f.write('>' + cid + '\n')
                    f.write(''.join(np.char.mod('%d', bvm[cid])) + '\n')
                f.close()
        except EnvironmentError as err:
            print "Unable to open write binary vector to file".format(err);
        return

    def write_all_files(self, base_name):
        infix = '_{0}_P{1}'.format(self.test_type, self.test_param)
        self.write_bv_to_file(self.contig_overcovered_bps, base_name + infix + '_OVER_BP.cov')
        self.write_bv_to_file(self.contig_undercovered_bps, base_name + infix + '_UNDER_BP.cov')
        self.write_bv_to_file(self.contig_overcovered_windows, base_name + infix + '_OVER_WIN.cov') # fixed suffix
        self.write_bv_to_file(self.contig_undercovered_windows, base_name + infix + '_UNDER_WIN.cov')
        return

    def to_string(self): # Typically one would overwrite _str()_, but still ok
        param_template = '{0:>30}{1:>15}\n'
        if self.test_type == 'Gaussian':
            title = 'Gaussian test with a standard deviation multiplier threshold of {0}\n'.format(self.test_param)
            params = 'Estimated parameters:\n'
            params += param_template.format('Mean (mu):', self.mu)
            params += param_template.format('Standard deviation (sigma):', self.sigma)
            params += param_template.format('High threshold (t_high):', self.t_high)
            params += param_template.format('Low threshold (t_low):', self.t_low)
        elif self.test_type == 'Percentile':
            title = 'Percentile test for {0}% most extreme values (top/bottom {1}%)\n'.format(100 * self.test_param, 50 * self.test_param)
            params = 'Estimated parameters:\n'
            params += param_template.format('High threshold (t_high):', self.t_high)
            params += param_template.format('Low threshold (t_low):', self.t_low)
        else:
            raise RuntimeError, "CoverageStatistics(): Internal error."
        template = '{0:>13} {1:>13} {2:>13}\n'
        table_head = self.__tuple_list_to_string__(template, [('CONTIG_ID', 'START_INDEX', 'END_INDEX')])
        over_int_list = [(cid, t[0], t[1]) for cid in self.contig_overcovered_intervals for t in self.contig_overcovered_intervals[cid]]
        under_int_list = [(cid, t[0], t[1]) for cid in self.contig_undercovered_intervals for t in self.contig_undercovered_intervals[cid]]
        table_over = 'The following regions have coverage above t_high:\n'
        table_over += table_head
        table_over += self.__tuple_list_to_string__(template, over_int_list)
        table_under = 'The following regions have coverage below t_low:\n'
        table_under += table_head
        table_under += self.__tuple_list_to_string__(template, under_int_list)
        return title + params + '\n' + table_over + '\n' + table_under

    def __init__(self, cov, test_type, test_param):
        self.test_type = test_type
        self.test_param = test_param
        if test_type == 'Gaussian':
            if test_param < 0:
                raise RuntimeError, "CoverageStatistics(): invalid standard deviation multiplier (expected non-negative value)."
            allvals = self.__get_all_values__(cov)
            self.mu = np.mean(allvals)
            self.sigma = np.std(allvals)
            self.t_high = self.mu + test_param * self.sigma
            self.t_low = self.mu - test_param * self.sigma
        elif test_type == 'Percentile':
            if test_param < 0:
                raise RuntimeError, "CoverageStatistics(): invalid percentile (expected float value in [0,1])."
            allvals = self.__get_all_values__(cov)
            allvals = np.sort(allvals)
            f = test_param / 2.0
            ind_high = np.floor((1 - f) * len(allvals))
            ind_low = np.ceil(f * len(allvals))
            self.t_high = np.float32(allvals[ind_high])
            self.t_low = np.float32(allvals[ind_low])
        else:
            raise RuntimeError, "CoverageStatistics(): unknown test type."
        # Test coverage values against t_high and t_low and extract problematic
        # regions for each contig
        self.contig_coverage = cov.contig_coverage
        self.contig_length = cov.contig_length
        self.contig_overcovered_windows = {}
        self.contig_undercovered_windows = {}
        self.contig_overcovered_bps = {}
        self.contig_undercovered_bps = {}
        self.contig_overcovered_intervals = {}
        self.contig_undercovered_intervals = {}
        for cid in cov.contig_coverage:
            win_over = cov.contig_coverage[cid] > self.t_high
            win_under = cov.contig_coverage[cid] < self.t_low
            bp_over = self.__map_win_to_bp__(win_over, cov, cid)
            bp_under = self.__map_win_to_bp__(win_under, cov, cid)
            int_over = self.__map_bp_to_int__(bp_over)
            int_under = self.__map_bp_to_int__(bp_under)
            # binary vectors
            self.contig_overcovered_windows[cid] = win_over
            self.contig_undercovered_windows[cid] = win_under
            self.contig_overcovered_bps[cid] = bp_over
            self.contig_undercovered_bps[cid] = bp_under
            # lists of intervals
            self.contig_overcovered_intervals[cid] = int_over
            self.contig_undercovered_intervals[cid] = int_under

'''
if __name__ == '__main__':
    samfile = SamFile.read('../../../../../tutorial/read_coverage/influenza-A.sam')
    cdata = parse_fasta_file('../../../../../data/influenza-A/influenza-A.assembly.fasta')
    bp_cov = ContigBPCoverage(samfile, cdata)
    w_cov = ContigWindowCoverage(bp_cov, 100, 1)
    st = CoverageStatistics(w_cov, 'Gaussian', 2.0)
    # st = CoverageStatistics(w_cov, 'Percentile', 0.05)
    print st.to_string()
    # for c in st.contig_overcovered_intervals:
    #    print st.contig_overcovered_intervals[c]
    # bp = np.array([True, False, True, True, False, True, True, True])
    # print st.__map_bp_to_int__(bp)
    # c = '8'
    # plt.plot(w_cov.contig_coverage[c])
    # plt.savefig('cov.png')
    # plt.plot(st.contig_overcovered_windows[c])
    # plt.savefig('win_over.png')
    # plt.plot(st.contig_undercovered_bps[c])
    # plt.savefig('bp_under.png')
    # plt.plot(st.contig_overcovered_bps[c])
    # plt.savefig('bp_over.png')
'''
