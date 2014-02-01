"""Tests for availability of ``pyPdf``.

Written by Jesse Bloom, 2013.
"""

import sys
import unittest
import mapmuts.weblogo


class TestPyPdfAvailable(unittest.TestCase):
    """Tests for availability of ``pyPdf``."""

    def test_Import(self):
        """Is ``pyPdf`` available?
        """
        sys.stderr.write("\nTesting if pyPdf is available...")
        failmsg = 'FAILED to find pyPdf in current search path. '\
                + 'Creation of sequence logos with overlays will not'\
                + ' be supported. Other aspects of '\
                + 'the mapmuts package can still be used.\n'
        succeedmsg = 'weblogo is available, sequence logos with overlays can be '\
                + 'created.'
        if mapmuts.weblogo.PyPdfAvailable():
            sys.stderr.write(succeedmsg)
        else:
            sys.stderr.write(failmsg)
        self.assertTrue(mapmuts.weblogo.PyPdfAvailable(), msg=failmsg)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
