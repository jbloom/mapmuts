"""Tests main.MakeAlignments and main.MakeAlignmentsPlots.


Written by Jesse Bloom, 2012.
"""


import time
import sys
import os
import unittest
import mapmuts.align
import mapmuts.sequtils
import mapmuts.main
import mapmuts.plot
import mapmuts.latex


class TestMakeAlignments(unittest.TestCase):
    """Tests main.MakeAlignments and main.MakeAlignmentsPlots."""

    def test_MakeAlignments(self):
        """Tests main.MakeAlignments and main.MakeAlignmentsPlots."""
        sys.stderr.write("\nTesting main.MakeAlignments...")
        r1files = ['R1.fastq.gz']
        r2files = ['R2.fastq.gz']
        a1 = mapmuts.sequtils.ReadFASTA('R1_trim3.fasta')[0][1]
        a2 = mapmuts.sequtils.ReadFASTA('R2_trim3.fasta')[0][1]
        minoverlap = 20
        maxa1m = maxa2m = 2
        maxn = 2
        maxrm = 2
        maxgenem = 4
        minq = 25
        gzipped = True
        applyfilter = True
        fullgene = mapmuts.sequtils.ReadFASTA('Aichi68-NP_amplicon.fasta')[0][1]
        generange = (61, 1555)
        outfileprefix = 'test'
        mapmuts.main.MakeAlignments(r1files, r2files, gzipped, applyfilter,\
                minq, fullgene, generange, a1, a2, maxn, minoverlap, maxrm,\
                maxa1m, maxa2m, maxgenem, outfileprefix)
        expectprefix = 'expected_%s' % outfileprefix
        for (suffix, Readf) in \
                [('alignmentstatistics', mapmuts.io.ReadAlignmentStatistics),\
                 ('insertlengths', mapmuts.io.ReadInsertLengths),\
                 ('R1mismatches', mapmuts.io.ReadReadMismatches),\
                 ('R2mismatches', mapmuts.io.ReadReadMismatches)]:
            actual = Readf(open("%s_%s.txt" % (outfileprefix, suffix)))
            expected = Readf(open("%s_%s.txt" % (expectprefix, suffix)))
            msg = 'FAILED: Actual and expected output differ for %s_%s.txt.'\
                    % (outfileprefix, suffix)
            if actual != expected:
                sys.stderr.write(msg)
            self.assertEqual(actual, expected, msg=msg)
        sys.stderr.write(" testing success.\nNow testing main.MakeAlignmentsPlots... ")
        if mapmuts.plot.PylabAvailable():
            if mapmuts.latex.PdflatexAvailable():
                latexsummary = True
            else:
                latexsummary = False
                sys.stderr.write('\nWARNING, cannot test latexsummary '\
                        + 'as pdflatex is not available...\n')
            mapmuts.main.MakeAlignmentsPlots(outfileprefix, minq, maxn,\
                    minoverlap, maxrm, maxa1m, maxa2m, maxgenem, latexsummary,\
                    summarytitle=outfileprefix)
            sys.stderr.write('Test completed successfully.')
        else:
            sys.stderr.write("\nWARNING, cannot test MakeAlignmentsPlots as "\
                    + 'pylab is not available.')


    def tearDown(self):
        """Remove created files."""
        files = ['test_log.txt', 'test_insertlengths.txt', 'test_alignmentstatistics.txt', 'test_alignments.txt.gz', 'test_R2mismatches.txt', 'test_R1mismatches.txt', 'test_alignmentstatistics.pdf', 'test_insertlengths.pdf', 'test_R1mismatches.pdf', 'test_makealignments_summary.tex', 'test_R2mismatches.pdf', 'test_makealignments_summary.pdf', 'test_makealignments_summary.log', 'test_makealignments_summary.aux']
        for f in files:
            if os.path.isfile(f):
                os.remove(f)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
