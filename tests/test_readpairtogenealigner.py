"""Tests mapmuts.align.ReadPairToGeneAligner

This test generates simulated data and then analyzes it with
ReadPairToGeneAligner and makes sure that the returned
resultsa are reeasonable.

Written by Jesse Bloom, 2012.
"""


import time
import sys
import unittest
import mapmuts.simulate
import mapmuts.align
import mapmuts.sequtils


class TestReadPairToGeneAligner(unittest.TestCase):
    """Tests mapmuts.align.ReadPairToGeneAligner."""

    def setUp(self):
        """Performs set up."""
        sys.stderr.write("\nTesting align.ReadPairToGeneAligner.")
        self.libsize = int(2e4)
        mapmuts.simulate.Seed(1)
        self.r1_trim3 = mapmuts.sequtils.ReadFASTA('R1_trim3.fasta')[0][1].upper()
        self.r2_trim3 = mapmuts.sequtils.ReadFASTA('R2_trim3.fasta')[0][1].upper()
        self.r2_trim3_rc =\
                mapmuts.sequtils.ReverseComplement(self.r2_trim3)
        self.fullgene = mapmuts.sequtils.ReadFASTA(\
                'Aichi68-NP_amplicon.fasta')[0][1].upper()
        self.generange = (61, 1555)
        self.readlength = 50
        self.insertlength = ('uniform', 30, 60)
        self.maxn = 2

    def test_AlignmentStatistics(self):
        """Tests if reads with errors can be aligned by ReadPairToGeneAligner."""
        sys.stderr.write('\nTesting alignment of reads with some errors.')
        sys.stderr.write('\nFirst simulating reads...')
        mutrate = 0.01
        (r1_reads, r2_reads) = mapmuts.simulate.SimulateIlluminaPE(\
                self.fullgene, self.r1_trim3, self.r2_trim3, self.readlength,
                self.insertlength, self.libsize, mutrate=mutrate)
        r1_reads = [s.upper() for s in r1_reads]
        r2_reads = [s.upper() for s in r2_reads]
        r2_reads_rc = mapmuts.sequtils.ReverseComplement(r2_reads)
        aligner = mapmuts.align.ReadPairToGeneAligner(self.fullgene,\
                self.generange, self.r1_trim3, self.r2_trim3_rc, self.maxn, 30,\
                maxrm=1, maxa1m=1, maxa2m=1, maxgenem=1, upcase=False)
        sys.stderr.write("\nNow using ReadPairToGeneAligner to align reads... ")
        start = time.clock()
        for (r1, r2) in zip(r1_reads, r2_reads_rc):
            a = aligner.Align(r1, r2)
        t = time.clock() - start
        sys.stderr.write("processed all %d reads in %.3f seconds.\n"\
                % (self.libsize, t))
        self.assertTrue(aligner.Statistics('nattempted') > aligner.Statistics(\
                'npaired'), msg='npaired >= nattempted')
        self.assertTrue(aligner.Statistics('npaired') > aligner.Statistics(\
                'nalignedfullgene'), msg='nalignedfullgene >= npaired')
        self.assertTrue(aligner.Statistics('nalignedfullgene') > aligner.Statistics(\
                'nalignedgene'), msg='nalignedgene >= nalignedfullgene')
        sys.stderr.write('Alignment statistics seem reasonable.')

    def test_PerfectReads(self):
        """Tests that perfect reads can all be aligned by ReadPairToGeneAligner."""
        sys.stderr.write("\nTesting alignment of perfect reads.")
        sys.stderr.write("\nFirst simulating reads...")
        gene = self.fullgene[self.generange[0] : self.generange[1]]
        (r1_reads, r2_reads) = mapmuts.simulate.SimulateIlluminaPE(\
                gene, self.r1_trim3, self.r2_trim3, self.readlength,
                self.insertlength, self.libsize, mutrate=0.0)
        r1_reads = [s.upper() for s in r1_reads]
        r2_reads = [s.upper() for s in r2_reads]
        r2_reads_rc = mapmuts.sequtils.ReverseComplement(r2_reads)
        aligner = mapmuts.align.ReadPairToGeneAligner(gene,\
                (0, len(gene)), self.r1_trim3, self.r2_trim3_rc, self.maxn, 30,
                maxrm=0, maxa1m=0, maxa2m=0, maxgenem=0, upcase=False)
        sys.stderr.write("\nNow using ReadPairToGeneAligner to align reads... ")
        start = time.clock()
        for (r1, r2) in zip(r1_reads, r2_reads_rc):
            a = aligner.Align(r1, r2)
            if not a:
                msg = 'FAILED. Could not align %s and %s.\n' % (r1, r2)
                sys.stderr.write(msg)
                self.assertTrue(a, msg=msg)
                break
        t = time.clock() - start
        sys.stderr.write("test passed, all %d reads aligned in %.3f seconds."\
                % (self.libsize, t))
        for key in ['nattempted', 'npaired', 'nalignedfullgene', 'nalignedgene']:
            self.assertTrue(aligner.Statistics(key) == self.libsize, msg=\
                    'Statistic for %s does not match library size.' % key)
        for insertlength in range(self.insertlength[1], self.insertlength[2] + 1):
            nexpected = self.libsize / float(self.insertlength[2] -\
                    self.insertlength[1] + 1)
            self.assertTrue(nexpected * 0.8 <= aligner.Statistics(\
                    'insertlength')[insertlength] <= nexpected * 1.2,
                    msg='Unexpected number of reads of length %d' % insertlength)

    def test_AlignmentIndices(self):
        """Tests that ReadPairToGeneAligner alignment indices are correct."""
        sys.stderr.write("\nTesting alignment indices.")
        sys.stderr.write("\nFirst simulating reads...")
        (r1_reads, r2_reads) = mapmuts.simulate.SimulateIlluminaPE(\
                self.fullgene, self.r1_trim3, self.r2_trim3, self.readlength,
                self.insertlength, self.libsize, mutrate=0.0)
        r1_reads = [s.upper() for s in r1_reads]
        r2_reads = [s.upper() for s in r2_reads]
        r2_reads_rc = mapmuts.sequtils.ReverseComplement(r2_reads)
        aligner = mapmuts.align.ReadPairToGeneAligner(self.fullgene,\
                self.generange, self.r1_trim3, self.r2_trim3_rc, self.maxn, 30,
                maxrm=0, maxa1m=0, maxa2m=0, maxgenem=0, upcase=False)
        sys.stderr.write("\nNow using ReadPairToGeneAligner to align reads... ")
        start = time.clock()
        naligned = 0
        gene = self.fullgene[self.generange[0] : self.generange[1]]
        for (r1, r2, r2_rc) in zip(r1_reads, r2_reads, r2_reads_rc):
            a = aligner.Align(r1, r2_rc)
            if isinstance(a, str):
                naligned += 1
                (ga, r1a, r2a) = a.split('\n')
                slice = gene[int(ga.split()[0]) - 1 : int(ga.split()[2])]
                self.assertTrue(ga.split()[1] == slice, msg='FAILED. Incorrect gene '\
                        + 'indices for:\n%s' % a)
                (r1i1, r1slice, r1i2) = (int(r1a.split()[0]), r1a.split()[1], int(r1a.split()[2]))
                r1slice = mapmuts.align.RemoveDots(slice, r1slice)
                (r2i1, r2slice, r2i2) = (int(r2a.split()[0]), r2a.split()[1], int(r2a.split()[2]))
                r2slice = mapmuts.align.RemoveDots(slice, r2slice)
                if r1i1 < r1i2:
                    # not reverse complement
                    self.assertTrue(r1[r1i1 - 1 : r1i2] == r1slice, msg='FAILED. '\
                            + "Incorrect r1 indices for:\n%s" % a)
                    r2slice_rc = mapmuts.sequtils.ReverseComplement(r2slice)
                    self.assertTrue(r2[r2i2 - 1 : r2i1] == r2slice_rc, msg='FAILED. '\
                            + 'Incorrect r2 indices for:\n%s' % a)
                elif r1i1 > r1i2:
                    # is reverse complement
                    r1slicerc = mapmuts.sequtils.ReverseComplement(r1slice)
                    self.assertTrue(r1[r1i2 - 1 : r1i1] == r1slicerc, msg='FAILED. '\
                            + 'Incorrect r1 indices for:\n%s' % a)
                    self.assertTrue(r2[r2i1 - 1 : r2i2] == r2slice, msg='FAILED. '\
                            + 'Incorrect r2 indices for:\n%s' % a)
        t = time.clock() - start
        sys.stderr.write("all indices correct for the %d of %d reads that aligned."\
                % (naligned, self.libsize) + ' Alignments completed in %.3f seconds'\
                % t)



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
