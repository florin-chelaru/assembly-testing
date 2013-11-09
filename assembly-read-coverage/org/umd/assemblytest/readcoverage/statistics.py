'''
Created on Nov 7, 2013

@author: jason & kostas
'''

import numpy as np

class StatisticalTest(object):
    
    '''
    Computes statistical data (mean, stdev), coverage percentiles
    for a given sequence<->contig mapping.
    '''
    
    '''
    Constructor docs
    '''
    def __init__(self, threshold, type_of_test="Gaussian"):
        
        # Sanity checking first....
        
        if type_of_test == "Gaussian":
            # return "bad" intervals with their μ and σ?
        else: # percentiles?
            # return "bad intervals iff they are above or below
            # certain percentiles?
            
            
    
     
        
        
        