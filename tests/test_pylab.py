"""Tests for availability of pylab / matplotlib.

Written by Jesse Bloom, 2012.
"""

import sys
import unittest
import mapmuts.plot


class TestPylabAvailable(unittest.TestCase):
    """Tests for availability of pylab / matplotlib."""

    def test_Import(self):
        """Can pylab / matplotlib be imported?
        """
        sys.stderr.write("\nTesting if pylab / matplotlib are present...")
        failmsg = 'FAILED to import pylab / matplotlib. Plotting '\
                + 'functions will not be supported. Other aspects of '\
                + 'the mapmuts package can still be used.\n'
        succeedmsg = 'successfully imported pylab / matplotlib. Plotting'\
                + ' is enabled.'
        if mapmuts.plot.PylabAvailable():
            sys.stderr.write(succeedmsg)
        else:
            sys.stderr.write(failmsg)
        self.assertTrue(mapmuts.plot.PylabAvailable(), msg=failmsg)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
