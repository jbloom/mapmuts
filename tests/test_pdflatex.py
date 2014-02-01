"""Tests for availability of pdflatex.

Written by Jesse Bloom, 2012.
"""

import sys
import unittest
import mapmuts.latex


class TestPdflatexAvailable(unittest.TestCase):
    """Tests for availability of pdflatex."""

    def test_Import(self):
        """Is pdflatex available?
        """
        sys.stderr.write("\nTesting if pdflatex is available...")
        failmsg = 'FAILED to find pdflatex in current search path. '\
                + 'Creation of summary pages with pdflatex will not'\
                + ' be supported. Other aspects of '\
                + 'the mapmuts package can still be used.\n'
        succeedmsg = 'pdflatex is available, summary pages can be '\
                + 'created.'
        if mapmuts.latex.PdflatexAvailable():
            sys.stderr.write(succeedmsg)
        else:
            sys.stderr.write(failmsg)
        self.assertTrue(mapmuts.latex.PdflatexAvailable(), msg=failmsg)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
