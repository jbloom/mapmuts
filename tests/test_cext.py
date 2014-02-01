"""Tests for availability of C extensions.

This script uses unittest to test of the C extensions to the
mapmuts package can be imported.

Written by Jesse Bloom, 2012.
"""

import sys
import unittest

class TestCExt(unittest.TestCase):
    """Tests for import of C extensions to mapmuts package."""

    def test_Import(self):
        """Can C extensions to mapmuts package be imported?
        """
        sys.stderr.write("\nTesting whether C extensions are present.")
        for ext in ['calign', 'csequtils']:
            try:
                sys.stderr.write("\nTrying to import %s..." % ext)
                module = __import__("mapmuts.%s" % ext)
                del module
                sys.stderr.write(" import succeeded.")
            except ImportError:
                msg = 'FAILED to import %s. The package will still'\
                        % ext + ' run in without this C extension, but'\
                        + ' it will be slower since it will have to use'\
                        + ' pure Python code.\n'
                sys.stderr.write(msg)
                self.assertTrue(False, msg=msg)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
