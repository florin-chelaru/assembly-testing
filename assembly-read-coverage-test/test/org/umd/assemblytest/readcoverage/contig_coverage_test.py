'''
Created on Nov 7, 2013

@author: jason & kostas
'''
import unittest
from org.umd.assemblytest.readcoverage.coverage import ContigData as CData;
from org.umd.assemblytest.readcoverage.samfile import SamFile;
from org.umd.assemblytest.readcoverage.coverage import ContigCoverage;

# It's not very straightforward to
# implement unit cases for full coverage tests.
class ContigCovTest(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def testName(self):
        pass

    def lengthsOfMaps(self):
        cdata = CData('../../../../../data/influenza-A/influenza-A.assembly.fasta');
        samfile = SamFile.read('../../../../../tutorial/read_coverage/influenza-A.sam');
        coverage = ContigCoverage(samfile, cdata, 100, 100);
        self.assertEquals(len(coverage.contig_coverage), len(coverage.contig_window_starting_points));

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
