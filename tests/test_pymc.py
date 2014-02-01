"""Tests for availability of pymc and numpy.

Written by Jesse Bloom, 2013.
"""

import sys
import unittest
import mapmuts.bayesian


class TestPymcAvailable(unittest.TestCase):
    """Tests for availability of numpy and pymc."""

    def test_Import(self):
        """Can pymc and numpy be imported?
        """
        sys.stderr.write("\nTesting if pymc and numpy are present...")
        failmsg = 'FAILED to import pymc or numpy. MCMC inference will not be supported'
        succeedmsg = 'successfully imported pync and numpy. MCMC inference is supported'
        if mapmuts.bayesian.PymcAvailable():
            sys.stderr.write(succeedmsg)
        else:
            sys.stderr.write(failmsg)
        self.assertTrue(mapmuts.bayesian.PymcAvailable(), msg=failmsg)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
