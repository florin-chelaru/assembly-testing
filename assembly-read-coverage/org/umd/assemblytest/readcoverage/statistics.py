'''
Created on Nov 7, 2013

@author: jason & kostas
'''
import numpy as np
#from samfile import SamFile
from coverage import ContigCoverage
# import matplotlib.pyplot as plt

class CoverageStats(object):
    
    def gaussiantest(self, contig_coverage, threshold):
    
        if threshold < 0:
            raise RuntimeError, "gaussiantest(): Threshold should be non-negative."
        idstocov = contig_coverage.contig_coverage
        idstolength = contig_coverage.contig_length
        win_length = contig_coverage.window_length
        idstowindowstartpts = contig_coverage.contig_window_starting_points
        allcovs = np.ndarray.flatten([idstocov[ids] for ids in idstocov])
        mu = np.mean(allcovs)
        stdev = np.std(allcovs)
        t_high = mu + threshold * stdev
        t_low = mu - threshold * stdev
        bad_windows = {} # maps ids to numpy arrays
        for c_id in idstocov:
            indices_high = idstocov[c_id] > t_high
            indices_low = idstocov[c_id] < t_low
            a = np.array(range(indices_high))
            b = np.array(range(indices_low))
            #bad_windows_above needs to be a list of triples
            bad_windows_above = [(idstowindowstartpts[c_id][index], 
                                  idstowindowstartpts[c_id][index] + win_length - 1,
                                  idstocov[c_id][index]) 
                                 for index in a[indices_high]]
            bad_windows[c_id] = 