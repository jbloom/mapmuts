#!python

"""Makes plots for each residue summarizing fraction of mutant amino acids.

This mapmuts_aafracsplots.py script is an analysis script for the 
mapmuts package.

For each residue in a protein, plots the fraction of non-wildtype
amino acids at each position for a set of different samples.

To run this script from the prompt, first create a text infile of the
format described below. Then simply type mapmuts_aafracsplots.py
followed by the infile name. For example, if the name is infile.txt,
type::

    mapmuts_aafracsplots.py infile.txt

If the script is not executable on your system, you can instead type::

    python mapmuts_aafracsplots.py infile.txt

This script is designed to be run after you have already run 
mapmuts_parsecounts.py on several samples. Each run of that program
generates a *_codoncounts.txt file, which is used as input for this
script.

This script will only work if pylab / matplotlib are available.
"""


import sys
import os
import time
import mapmuts.plot
import mapmuts.io
import mapmuts.sequtils


def main():
    """Main body of script."""
    print "Beginning execution of mapmuts_aafracsplots.py..."
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
    printprogress = plotfileprefix = None
    infiles = []
    names = []
    for line in lines:
        entries = line.split(None, 1)
        if len(entries) != 2:
            raise ValueError("line must contain at least two entries:\n%s" % line)
        if entries[0].strip() == 'plotfileprefix':
            if plotfileprefix != None:
                raise ValueError("Duplicate plotfileprefix keys")
            plotfileprefix = entries[1].strip()
        elif entries[0].strip() == 'printprogress':
            if printprogress != None:
                raise IOError("Duplicate entries for printprogress")
            printprogress = entries[1].strip()
            if printprogress == 'False':
                printprogress = False
            elif printprogress == 'True':
                printprogress = True
            else:
                raise ValueError("printprogress must be True or False")
        else:
            infile = entries[0].strip()
            names.append(entries[1].strip())
            if not os.path.isfile(infile):
                raise IOError("Cannot find specified file %s" % infile)
            infiles.append(infile)
    assert len(infiles) == len(names)
    if not infiles:
        raise ValueError("infile fails to specify any samples")
    if not mapmuts.plot.PylabAvailable():
        raise OSError("pylab / matplotlib is not available, so cannot run this script.")
    if not plotfileprefix:
        raise ValueError("infile did not specify plotfileprefix")
    if printprogress == None:
        raise ValueError("infile did not specify printprogress")
    dir = os.path.dirname(plotfileprefix)
    if dir:
        if not os.path.isdir(dir):
            raise IOError("plotfileprefix specifies directory of %s" % dir\
                    + ' but this directory does not exist. You must create'\
                    + ' it manually.')
    counts = []
    maxcodon = None
    for infile in infiles:
        codon_counts = mapmuts.io.ReadCodonCounts(open(infile))
        mapmuts.sequtils.ClassifyCodonCounts(codon_counts)
        counts.append(codon_counts)
        imaxcodon = max([i for i in codon_counts.keys() if isinstance(i, int)])
        if maxcodon != None:
            if imaxcodon != maxcodon:
                raise ValueError("codoncounts files specify different lengths")
        else:
            maxcodon = imaxcodon
    if not maxcodon:
        raise ValueError("codoncounts files specify empty sequence")
    for icodon in range(1, maxcodon + 1):
        wtaas = {}
        for codon_counts in counts:
            wtcodon = codon_counts[icodon]['WT']
            wtaa = mapmuts.sequtils.Translate([('head', wtcodon)])[0][1]
            if not wtaa:
                wtaa = 'STOP'
            wtaas[wtaa] = True
        for wtaa in wtaas.keys():
            aas = [aa for aa in mapmuts.sequtils.AminoAcids() if aa != wtaa]
            if wtaa != 'STOP':
                aas.append('STOP')
            title = "%s%d" % (wtaa, icodon)
            plotfile = "%s_%s.pdf" % (plotfileprefix, title)
            if printprogress:
                print "Making plot %s..." % plotfile
                sys.stdout.flush()
            mapmuts.plot.PlotAAFracs(counts, names, icodon, plotfile, aas, title)
    if printprogress:
        print "Plotting completed."
    print "Script complete."


if __name__ == '__main__':
    main() # run the script
