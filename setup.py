"""Setup script for mapmuts.

This script uses distutils, the standard python mechanism for installing
packages. To build, test, and install the package, use the following
commands::

    python setup.py build
    python setup.py test
    python setup.py install

If the user does not have permissions to write to the install directory,
the last command may need to be replaced by::

    sudo python setup.py install

The test command runs a variety of tests to check that the program
appears to working properly. This is optional and will take some time
to run, but is probably a good idea -- particularly if you are doing
further development / changes to the code.

Written by Jesse Bloom.
"""


import sys
import os
from distutils.core import setup
from distutils.core import Extension
from distutils.core import Command


# create a class to handle the 'test' command
class TestCommand(Command):
    """Run all of the tests for this package.

    This is an automatic test run class to make distutils perform the
    package testing. To run these tests, type:

    python setup.py build
    python setup.py test
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        """Run the test script tests/run_tests.py"""
        currdir = os.getcwd()
        os.chdir('tests')
        sys.path.insert(0, '')
        import run_tests
        run_tests.main()
        os.chdir(currdir)

# list of C extensions
calign = Extension('mapmuts.calign', sources=['src/calign.c'])
csequtils = Extension('mapmuts.csequtils', sources=['src/csequtils.c'])

# main setup command
setup(
    name = 'mapmuts', 
    version = '1.1', 
    author = 'Jesse D. Bloom', 
    author_email = 'jbloom@fhcrc.org', 
    url = 'https://github.com/jbloom/mapmuts', 
    description = 'Analyzes mutant library deep-sequencing results.',
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: Free for non-commercial use',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    platforms = 'Tested on Mac OS X.',
    packages = ['mapmuts'],
    package_dir = {'mapmuts':'src'},
    ext_modules = [calign, csequtils],
    cmdclass = {'test':TestCommand},
    scripts = [
            'scripts/mapmuts_makealignments.py',
            'scripts/mapmuts_alignmentsummaryplot.py',
            'scripts/mapmuts_parsecounts.py',
            'scripts/mapmuts_parsesummaryplots.py',
            'scripts/mapmuts_countparsedmuts.py',
            'scripts/mapmuts_aafracsplots.py',
            'scripts/mapmuts_inferpreferences.py',
            'scripts/mapmuts_inferdifferentialpreferences.py',
            'scripts/mapmuts_preferencescorrelate.py',
            'scripts/mapmuts_preferencemeans.py',
            'scripts/mapmuts_inferenrichment.py',
            'scripts/mapmuts_enrichmentcorrelate.py',
            'scripts/mapmuts_enrichmentgeomeans.py',
            'scripts/mapmuts_siteprofileplots.py',
            'scripts/mapmuts_entropycomparison.py',
            ],
    package_data = {'mapmuts':['template.eps']}, # template from weblogo version 3.4
)
