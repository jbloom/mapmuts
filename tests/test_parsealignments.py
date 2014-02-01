"""Tests main.ParseNTCodonCounts.


Written by Jesse Bloom, 2012.
"""


import time
import sys
import os
import unittest
import mapmuts.sequtils
import mapmuts.main


class TestParseAlignments(unittest.TestCase):
    """Tests main.ParseNTCodonCounts."""

    def test_ParseNTCodonCounts(self):
        """Tests main.ParseNTCodonCounts."""
        sys.stderr.write("\nTesting main.ParseNTCodonCounts...")
        gzipped = True
        fullgene = mapmuts.sequtils.ReadFASTA('Aichi68-NP_amplicon.fasta')[0][1]
        generange = (61, 1555)
        gene = fullgene[generange[0] : generange[1]].upper()
        r1exclude = [1, 23, 27]
        r2exclude = [32]
        outfileprefix = 'test_parsealignments'
        alignmentfile = 'testparse_alignments.txt.gz'
        logfile = 'test_parsealignments.log'
        log = open(logfile, 'w')
        mapmuts.main.ParseNTCodonCounts(alignmentfile, outfileprefix, gene,\
                r1exclude, r2exclude, log=log)
        log.close()
        expectprefix = 'expected_%s' % outfileprefix
        for (suffix, Readf) in \
                [('ntcounts', mapmuts.io.ReadNTCounts),\
                 ('codoncounts', mapmuts.io.ReadCodonCounts)]:
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
            mapmuts.main.ParseNTCodonCountsPlots(outfileprefix, r1exclude,\
                    r2exclude, latexsummary, summarytitle=outfileprefix)
            sys.stderr.write('Test completed successfully.')
        else:
            sys.stderr.write("\nWARNING, cannot test MakeAlignmentsPlots as "\
                    + 'pylab is not available.')


    def tearDown(self):
        """Removes created files."""
        files = ['test_parsealignments_ntcounts.txt', 'test_parsealignments_codoncounts.txt', 'test_parsealignments.log', 'test_parsealignments_codondepth.pdf', 'test_parsealignments_syn-ns-dist.pdf', 'test_parsealignments_parsesummary.tex', 'test_parsealignments_nmutspercodon-dist.pdf', 'test_parsealignments_parsesummary.pdf', 'test_parsealignments_parsesummary.log', 'test_parsealignments_parsesummary.aux']
        for f in files:
            if os.path.isfile(f):
                os.remove(f)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
