#!python

"""Makes summary plot after several runs of mapmuts_makealignments.py.

This mapmuts_alignmentsummaryplot.py script is an analysis script for the 
mapmuts package.

This script creates a stacked bar graph that summarizes the overall
number of reads for each sample and how many reads were retained
and removed during the alignment process.

To run this script from the prompt, first create a text infile of the
format described below. Then simply type mapmuts_alignmentsummaryplot.py
followed by the infile name. For example, if the name is infile.txt,
type::

    mapmuts_alignmentsummaryplot.py infile.txt

If the script is not executable on your system, you can instead type::

    python mapmuts_alignmentsummaryplot.py infile.txt

This script will only work if pylab / matplotlib are available.
"""


import sys
import os
import time
import mapmuts.plot
import mapmuts.io


def main():
    """Main body of script."""
    print "Beginning execution of mapmuts_alignmentsummaryplot.py..."
    mapmuts.io.PrintVersions(sys.stdout)
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument"\
                + ' specifying the name of the input file.')
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    lines = [line for line in open(infilename).readlines() if \
            line[0] != '#' and not line.isspace()]
    plotfile = None
    files = []
    samples = []
    for line in lines:
        entries = line.split(None, 1)
        if len(entries) != 2:
            raise ValueError("line must contain two entries:\n%s" % line)
        if entries[0].strip() == 'plotfile':
            plotfile = entries[1].strip()
        else:
            file = entries[0].strip()
            sample = entries[1].strip()
            if not os.path.isfile(file):
                raise ValueError("file %s does not exist, but is specified"\
                        % file + ' by:\n%s' % line)
            files.append(file)
            samples.append(sample)
    assert len(files) == len(samples)
    if len(files) < 2:
        raise ValueError("infile specifies < 2 samples")
    if not mapmuts.plot.PylabAvailable():
        raise OSError("pylab / matplotlib is not available, so cannot run this script.")
    if not plotfile:
        raise ValueError("infile did not specify plotfile")
    mapmuts.plot.PlotAlignmentStatistics(files, samples, plotfile)
    print "Script complete"


if __name__ == '__main__':
    main() # run the script
