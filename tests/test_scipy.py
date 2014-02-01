"""Tests for availability of scipy.

Written by Jesse Bloom, 2013.
"""

import sys
import unittest
import mapmuts.bayesian


class TestScipyAvailable(unittest.TestCase):
    """Tests for availability of scipy."""

    def test_Import(self):
        """Can scipy be imported?
        """
        sys.stderr.write("\nTesting if scipy is present...")
        failmsg = 'FAILED to import scipy. The MCMC inference will be much less efficient without this package.'
        succeedmsg = 'successfully imported scipy. MCMC inference should run efficiently.'
        if mapmuts.bayesian.ScipyAvailable():
            sys.stderr.write(succeedmsg)
        else:
            sys.stderr.write(failmsg)
        self.assertTrue(mapmuts.bayesian.ScipyAvailable(), msg=failmsg)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
