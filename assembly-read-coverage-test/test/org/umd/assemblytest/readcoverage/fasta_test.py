'''
Created on Nov 7, 2013

@author: jason & kostas
'''
import unittest
from org.umd.assemblytest.readcoverage.coverage import ContigData as CData;

class FastaTest(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def testName(self):
        pass

    def testParsing(self):
        cd = CData('resources/dummy_fasta.txt');
        self.assertTrue('1' in cd.contig_length)
        self.assertTrue('2' in cd.contig_length)
        self.assertTrue('non_numerical_id' in cd.contig_length)
        self.assertTrue('10' in cd.contig_length)

        # for key in cd.contig_length:
        #   print "Key {0} maps to: {1}".format(key, cd.contig_length[key]);
        self.assertEquals(cd.contig_length['1'], 5)
        self.assertEquals(cd.contig_length['2'], 20)
        self.assertEquals(cd.contig_length['non_numerical_id'], 0)
        self.assertEquals(cd.contig_length['10'], 10)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
