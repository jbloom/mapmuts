"""Tests for availability of ``weblogo``.

Written by Jesse Bloom, 2013.
"""

import sys
import unittest
import mapmuts.weblogo


class TestWebLogoAvailable(unittest.TestCase):
    """Tests for availability of ``weblogo``."""

    def test_Import(self):
        """Is ``weblogo`` available?
        """
        sys.stderr.write("\nTesting if weblogo is available...")
        failmsg = 'FAILED to find weblogo in current search path. '\
                + 'Creation of sequence logos with weblogo will not'\
                + ' be supported. Other aspects of '\
                + 'the mapmuts package can still be used.\n'
        succeedmsg = 'weblogo is available, sequence logos can be '\
                + 'created.'
        if mapmuts.weblogo.WebLogoAvailable():
            sys.stderr.write(succeedmsg)
        else:
            sys.stderr.write(failmsg)
        self.assertTrue(mapmuts.weblogo.WebLogoAvailable(), msg=failmsg)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
