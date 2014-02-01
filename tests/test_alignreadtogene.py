"""Tests the alignment of reads to a gene.

This script uses unittest to test the capability of the mapmuts
package to lign reads to a gene using the 
mapmuts.align.AlignReadToGene function. It first creates simulated
reads, and then aligns them to the gene. It also prints information
about the time needed to perform these operations.

Written by Jesse Bloom, 2012.
"""

import time
import sys
import math
import unittest
import random
import mapmuts.sequtils
import mapmuts.align
import mapmuts.simulate


class TestAlignReadToGene(unittest.TestCase):
    """Simulates Illumina paired-end reads and tests alignment.

    Writes information on timing to stderr.
    """

    def setUp(self):
        """Performs set up.

        Generates simulated reads for analysis by testing functions.
        """
        sys.stderr.write("\nTesting AlignReadToGene: alignment of reads"\
                + " to a gene.")
        self.libsize = 1e4
        mapmuts.simulate.Seed(1)
        self.gene = mapmuts.sequtils.ReadFASTA('Aichi68-NP.fasta')[0][1]
        self.gene = self.gene.upper()
        self.gene_rc = mapmuts.sequtils.ReverseComplement(self.gene)
        self.readlength = 50
        self.mutrate = 0.01
        self.insertlength = ('uniform', 50, 50)

    def test_AllShouldAlign(self):
        """Simulated reads align with sufficiently generous mismatch cutoffs.
        """
        sys.stderr.write("\nFirst simulating reads...")
        (self.r1_reads, self.r2_reads) = mapmuts.simulate.SimulateIlluminaPE(
                self.gene, '', '', self.readlength, 
                self.insertlength, self.libsize, self.mutrate)
        sys.stderr.write("\nTesting simulated reads alignment with " +\
                'sufficiently generous mismatch cutoffs... ')
        start_clock = time.clock()
        for r in self.r1_reads:
            a = mapmuts.align.AlignReadToGene(r, self.gene, self.gene_rc, 
                    7, upcase=False)
            self.assertTrue(a, msg='Alignment failed for %s.' % r)
            if not a:
                sys.stderr.write("FAILED to align read!")
                break
        else:
            t = time.clock() - start_clock
            sys.stderr.write("aligned all reads in %.3f seconds" % t)

    def test_NoneShouldAlign(self):
        """Simulated random reads should not align."""
        sys.stderr.write("\nFirst creating random reads...")
        random_reads = []
        for i in range(int(self.libsize)):
            s = ''.join([random.choice(['A', 'T', 'C', 'G']) for j\
                    in range(int(self.readlength))])
            random_reads.append(s)
        sys.stderr.write("\nTesting that random reads do not align... ")
        start_clock = time.clock()
        for r in random_reads:
            a = mapmuts.align.AlignReadToGene(r, self.gene, self.gene_rc,
                    5, upcase=False)
            if a:
                msg = 'FAILED: spurious alignment for %s and %s.' % (r1, r2)
                self.assertFalse(a, msg=msg)
                sys.stderr.write(msg)
        else:
            t = time.clock() - start_clock
            sys.stderr.write("examined reads in %.3f seconds"\
                    % t + ', test passed as none aligned.')

    def test_AlignStatistics(self):
        """Reasonable read alignment statistics for various cutoffs."""
        sys.stderr.write("\nFirst simulating reads...")
        (self.r1_reads, self.r2_reads) = mapmuts.simulate.SimulateIlluminaPE(
                self.gene, '', '', self.readlength, 
                self.insertlength, self.libsize, self.mutrate)
        sys.stderr.write("\nTesting simulated reads alignment " +\
                'statistics for various mismatch cutoffs...\n')
        sys.stderr.write("Reads of length %d were generated with %.3f"\
                % (self.readlength, self.mutrate) + ' error rate.')
        nexpected = 0
        fexpected = math.exp(-self.mutrate * self.readlength)
        for maxm in range(4):
            naligned = 0
            start_clock = time.clock()
            for r in self.r1_reads:
                a = mapmuts.align.AlignReadToGene(r, self.gene,
                        self.gene_rc, maxm, upcase=False)
                if a:
                    naligned += 1
            t = time.clock() - start_clock
            sys.stderr.write("\nIn %.3f seconds, aligned %d of %d"\
                    % (t, naligned, self.libsize) + ' with <= %s' %\
                    maxm + ' mismatches.')
            if maxm:
                fexpected *= (self.mutrate * self.readlength) / float(maxm)
            nexpected += fexpected * self.libsize
            if 0.95 * nexpected < naligned < 1.05 * nexpected:
                sys.stderr.write("\nThis is close to expected number of %d."\
                        % nexpected)
            else:
                msg = 'FAILED. Expected to align %d, '\
                        % nexpected + 'but instead aligned %d.' % naligned
                sys.stderr.write("\n%s" % msg)
                self.assertTrue(False, msg=msg)



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
