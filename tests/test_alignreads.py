"""Tests the alignment of paired-end reads.

This script uses unittest to test the capability of the mapmuts
package to align paired-end reads and remove adaptor sequences. It
first creates simulated reads, and then analyzes them to make sure that
the overlap and adaptors are correctly identified. It also prints
information about the time needed to perform these operations.

Written by Jesse Bloom, 2012.
"""

import time
import sys
import math
import random
import unittest
import mapmuts.sequtils
import mapmuts.align
import mapmuts.simulate


class TestAlignReads(unittest.TestCase):
    """Simulates Illumina paired-end reads and tests alignment.

    Writes information on timing to stderr.
    """

    def setUp(self):
        """Performs set up.

        Generates simulated reads for analysis by testing functions.
        """
        sys.stderr.write("\nTesting alignment of simulated Illumina " +\
                "paired-end reads.")
        self.libsize = 2e4
        sys.stderr.write('\nFirst simulating library of %d reads...' %\
                int(self.libsize))
        mapmuts.simulate.Seed(1)
        self.r1_trim3 = mapmuts.sequtils.ReadFASTA('R1_trim3.fasta')[0][1]
        self.r2_trim3 = mapmuts.sequtils.ReadFASTA('R2_trim3.fasta')[0][1]
        self.gene = mapmuts.sequtils.ReadFASTA('Aichi68-NP.fasta')[0][1]
        self.readlength = 50
        self.mutrate = 0.01
        self.insertlength = ('uniform', 30, 60)
        (self.r1_reads, self.r2_reads) = mapmuts.simulate.SimulateIlluminaPE(
                self.gene, self.r1_trim3, self.r2_trim3, self.readlength, 
                self.insertlength, self.libsize, self.mutrate)
        self.r2_reads_rc =\
                mapmuts.sequtils.ReverseComplement(self.r2_reads)
        self.r2_trim3_rc =\
                mapmuts.sequtils.ReverseComplement([self.r2_trim3])[0]

    def test_AllShouldAlign(self):
        """Simulated reads align with sufficiently generous mismatch cutoffs.
        """
        self.r1_reads_caps = [''.join([nt for nt in r1 if nt.istitle()])
                for r1 in self.r1_reads]
        self.r2_reads_rc_caps = [''.join([nt for nt in r2 if nt.istitle()])
                for r2 in self.r2_reads_rc]
        self.r1_reads = [s.upper() for s in self.r1_reads]
        self.r2_reads_rc = [s.upper() for s in self.r2_reads_rc]
        sys.stderr.write("\nTesting simulated reads alignment with " +\
                'sufficiently generous mismatch cutoffs... ')
        start_clock = time.clock()
        for (r1, r2, r1_caps, r2_caps) in zip(self.r1_reads,
                self.r2_reads_rc, self.r1_reads_caps,
                self.r2_reads_rc_caps):
            a = mapmuts.align.AlignReads(r1, r2, self.r1_trim3,
                    self.r2_trim3_rc, 30, 8, 4, 4, upcase=False)
            if not a:
                sys.stderr.write("FAILED to align reads:\nr1=%s\nr2=%s"\
                        % (r1, r2))
                self.assertTrue(a, msg='Alignment failed for %s and %s. ' %\
                        (r1, r2))
                break
            else:
                r1slice = r1[a[0][0] : a[0][1]]
                r2slice = r2[a[1][0] : a[1][1]]
                if len(r1_caps) > len(r1slice):
                    r1_caps = r1_caps[a[0][0] : ]
                if len(r2_caps) > len(r2slice):
                    r2_caps = r2_caps[ : a[1][1]]
                if not (r1slice == r1_caps and r2slice == r2_caps):
                    msg='Failed to match full read.\n' +\
                            '%s = r1_caps\n%s = r1slice\n' % (r1_caps, r1slice)+\
                            '%s = r2_caps\n%s = r2slice\n' % (r2_caps, r2slice) +\
                            '\nr1 = %s\nr2 = %s\nExtracted alignment was %s.'\
                            % (r1, r2, str(a))
                    sys.stderr.write(msg)
                    self.assertTrue(r1slice == r1_caps and r2slice == r2_caps,\
                            msg=msg)
        else:
            t = time.clock() - start_clock
            sys.stderr.write("aligned all reads in %.3f seconds" % t)

    def test_NoneShouldAlign(self):
        """Simulated random reads should not align."""
        random_reads = []
        for i in range(int(self.libsize)):
            s = ''.join([random.choice(['A', 'T', 'C', 'G']) for j\
                    in range(int(self.readlength))])
            random_reads.append(s)
        self.r1_reads = [s.upper() for s in self.r1_reads]
        sys.stderr.write("\nTesting that random reads do not align... ")
        start_clock = time.clock()
        for (r1, r2) in zip (self.r1_reads, random_reads):
            a = mapmuts.align.AlignReads(r1, r2, self.r1_trim3,
                    self.r2_trim3_rc, 30, 7, 4, 4, upcase=False)
            if a:
                msg = 'FAILED: spurious alignment for %s and %s.' % (r1, r2)
                self.assertFalse(a, msg=msg)
                sys.stderr.write(msg)
        else:
            t = time.clock() - start_clock
            sys.stderr.write("examined reads in %.3f seconds "\
                    % t + ', test passed as none aligned.')

    def test_AlignStatistics(self):
        """Reasonable read alignment statistics for various cutoffs."""
        sys.stderr.write("\nTesting simulated reads alignment " +\
                'statistics for various mismatch cutoffs...\n')
        sys.stderr.write("Reads of length %d were generated with %.3f"\
                % (self.readlength, self.mutrate) + ' error rate.\n')
        naligned = 0
        start_clock = time.clock()
        self.r1_reads = [s.upper() for s in self.r1_reads]
        self.r2_reads_rc = [s.upper() for s in self.r2_reads_rc]
        for (r1, r2) in zip(self.r1_reads, self.r2_reads_rc):
            a = mapmuts.align.AlignReads(r1, r2, self.r1_trim3,
                    self.r2_trim3_rc, 30, 0, 0, 0, upcase=False)
            if a:
                naligned += 1
                msg = 'Problem in specified overlap:\nr1=%s' % r1\
                        + '\nr2=%s\ngot %s\n' % (r2, a)
                self.assertTrue(r1[a[0][0] : a[0][1]] == r2[a[1][0] :\
                        a[1][1]], msg=msg)
                if r1[a[0][0] : a[0][1]] != r2[a[1][0] : a[1][1]]:
                    sys.stderr.write(msg)
        t = time.clock() - start_clock
        sys.stderr.write("In %.3f seconds, aligned %d of %d"\
                % (t, naligned, self.libsize) + ' with no mismatches.\n')
        nexpected = math.exp(-2 * self.mutrate * self.readlength)\
                * self.libsize
        if 0.9 * nexpected < naligned < 1.1 * nexpected:
            sys.stderr.write("This is close to expected number of %d."\
                    % nexpected)
        else:
            msg = 'FAILED. Expected to align %d with no mismatches, '\
                    % nexpected + 'but instead aligned %d.' % naligned
            sys.stderr.write(msg)
            self.assertTrue(False, msg=msg)
        for (nr, na1, na2) in [(1, 0, 0), (1, 1, 0), (1, 1, 1)]:
            nalignedprev = naligned
            naligned = 0
            start_clock = time.clock()
            for (r1, r2) in zip(self.r1_reads, self.r2_reads_rc):
                a = mapmuts.align.AlignReads(r1, r2, self.r1_trim3,
                        self.r2_trim3_rc, 30, nr, na1, na2, upcase=False)
                if a:
                    naligned += 1
            t = time.clock() - start_clock
            sys.stderr.write("\nIn %.3f seconds, aligned %d of %d" % (t,\
                    naligned, self.libsize) + ' with %d read ' % nr +\
                    'mismatches, %d adaptor 1 mismatches, and %d '\
                    % (na1, na2) + 'adaptor 2 mismatches.')
            msg = 'FAILED. Did not align more reads than with previous'\
                    + 'more lenient cutoff: %d versus %d.' % (naligned, \
                    nalignedprev)
            if naligned <= nalignedprev:
                sys.stderr.write("\n%s" % msg)
            self.assertTrue(naligned > nalignedprev, msg=msg)



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
