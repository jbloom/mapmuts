#!python

"""Makes summary plots after several runs of mapmuts_makeparsecounts.py.

This mapmuts_parsesummaryplots.py script is an analysis script for the 
mapmuts package.

It creates stacked bar graphs giving the fraction of all sites that
mutation types. This allows you to compare mutation frequences among samples.

To run this script from the prompt, first create a text infile of the
format described below. Then simply type mapmuts_parsesummaryplots.py
followed by the infile name. For example, if the name is infile.txt,
type::

    mapmuts_parsesummaryplots.py infile.txt

If the script is not executable on your system, you can instead type::

    python mapmuts_parsesummaryplots.py infile.txt

This script is designed to be run after you have already run 
mapmuts_parsecounts.py on several samples. Each run of that program
generates a *_codoncounts.txt and a *_ntcounts.txt file. This script reads
those files and creates a stacked bar graph that summarizes the output
for the several samples.

This script will only work if pylab / matplotlib are available.
"""


import sys
import os
import time
import mapmuts.plot
import mapmuts.io


def main():
    """Main body of script."""
    print "Beginning execution of mapmuts_parsesummaryplots.py..."
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
    textwritefracs_codontypes = textwritefracs_codonnmuts = writefracs = plotfileprefix = None
    formatting = {}
    codonfiles = []
    codonnmutfiles = []
    ntfiles = []
    codonsamples = []
    codonnmutsamples = []
    ntsamples = []
    pairedcodonplot = False
    for line in lines:
        entries = line.split(None, 1)
        if len(entries) != 2:
            raise ValueError("line must contain at least two entries:\n%s" % line)
        if entries[0].strip() == 'plotfileprefix':
            if plotfileprefix != None:
                raise ValueError("Duplicate plotfileprefix keys")
            plotfileprefix = entries[1].strip()
        elif entries[0].strip() == 'figwidth':
            formatting['figwidth'] = float(entries[1])
        elif entries[0].strip() == 'writefracs':
            if writefracs != None:
                raise ValueError("Duplicate writefracs key")
            if entries[1].strip() == 'True':
                writefracs = True
            elif entries[1].strip() == 'False':
                writefracs = False
            else:
                raise ValueError("writefracs must be True or False on line:\n%s" % line)
        elif entries[0].strip() == 'manual_ntfracs':
            fracs = entries[1].split(None, 6)
            if len(fracs) != 7:
                raise ValueError("Invalid number of entries for:\n%s" % line)
            ntsamples.append(fracs[-1].strip())
            d = {}
            for (key, x) in zip(mapmuts.sequtils.NTMutTypes(), fracs[ : -1]):
                try:
                    d[key] = float(x)
                except TypeError:
                    raise ValueError("%s specifies invalid fraction on:\n%s" % (x, line))
            ntfiles.append(d)
        elif entries[0].strip() == 'manual_codontypes':
            fracs = entries[1].split(None, 3)
            if len(fracs) != 4:
                raise ValueError("Invalid number of entries for:\n%s" % line)
            codonsamples.append(fracs[-1].strip())
            d = {}
            for (key, x) in zip(['synonymous', 'nonsynonymous', 'stop codon'], fracs[ : -1]):
                try:
                    d[key] = float(x)
                except TypeError:
                    raise ValueError("%s specifies invalid fraction on:\n%s" % (x, line))
            codonfiles.append(d)
        elif entries[0].strip() == 'manual_codonnmuts':
            fracs = entries[1].split(None, 3)
            if len(fracs) != 4:
                raise ValueError("Invalid number of entries for:\n%s" % line)
            codonnmutsamples.append(fracs[-1].strip())
            d = {}
            for (key, x) in zip(['1 nucleotide mutation', '2 nucleotide mutations', '3 nucleotide mutations'], fracs[ : -1]):
                try:
                    d[key] = float(x)
                except TypeError:
                    raise ValueError("%s specifies invalid fraction on:\n%s" % (x, line))
            codonnmutfiles.append(d)
        elif entries[0].strip() == 'textwritefracs':
            prefix = entries[1].strip()
            if prefix.upper() not in ['NONE', 'FALSE']:
                textwritefracs_codontypes = "%s_codontypes.txt" % prefix
                textwritefracs_codonnmuts = "%s_codonnmuts.txt" % prefix
        elif entries[0].strip() == 'pairedcodonplot':
            value = entries[1].strip().upper()
            if value in ['NONE', 'FALSE']:
                pairedcodonplot = False
            elif value == 'TRUE':
                pairedcodonplot = True
            else:
                raise ValueError("Invalid value for pairedcodonplot")
        else:
            prefix = entries[0].strip()
            sample = entries[1].strip()
            codonfile = '%s_codoncounts.txt' % prefix
            ntfile = '%s_ntcounts.txt' % prefix
            if not os.path.isfile(codonfile):
                raise ValueError("Cannot find file %s" % codonfile)
            if not os.path.isfile(ntfile):
                raise ValueError("Cannot find file %s" % ntfile)
            codonfiles.append(codonfile)
            codonnmutfiles.append(codonfile)
            ntfiles.append(ntfile)
            codonsamples.append(sample)
            codonnmutsamples.append(sample)
            ntsamples.append(sample)
    assert len(codonfiles) == len(codonsamples)
    assert len(ntfiles) == len(ntsamples)
    if len(ntsamples) < 2:
        raise ValueError("infile specifies < 2 ntsamples")
    if len(codonsamples) < 2:
        raise ValueError("infile specifies < 2 codonsamples")
    if not mapmuts.plot.PylabAvailable():
        raise OSError("pylab / matplotlib is not available, so cannot run this script.")
    if not plotfileprefix:
        raise ValueError("infile did not specify plotfileprefix")
    if writefracs == None:
        raise ValueError("infile did not specify writefracs")
    mapmuts.plot.PlotMutFracs(codonfiles, codonsamples, "%s_codontypes.pdf" % plotfileprefix, keys=['synonymous', 'nonsynonymous', 'stop codon'], writefracs=writefracs, formatting=formatting, textwritefracs=textwritefracs_codontypes)
    mapmuts.plot.PlotMutFracs(codonnmutfiles, codonnmutsamples, "%s_codonnmuts.pdf" % plotfileprefix, keys=['1 nucleotide mutation', '2 nucleotide mutations', '3 nucleotide mutations'], writefracs=writefracs, formatting=formatting, textwritefracs=textwritefracs_codonnmuts)
    mapmuts.plot.PlotNTMutFracs(ntfiles, ntsamples, "%s_ntfracs.pdf" % plotfileprefix)
    if pairedcodonplot:
        # get all samples specified in both codonnmut and codon
        samples_inboth = [sample for sample in codonsamples if sample in codonnmutsamples]
        infiles_inboth = []
        for sample in samples_inboth:
            itypes = codonsamples.index(sample)
            inmuts = codonnmutsamples.index(sample)
            if codonfiles[itypes] == codonnmutfiles[inmuts]:
                infiles_inboth.append(codonfiles[itypes])
            elif isinstance(codonfiles[itypes], dict) and isinstance(codonnmutfiles[inmuts], dict):
                d = dict([(key, value) for (key, value) in codonfiles[itypes].items() + codonnmutfiles[inmuts].items()])
                infiles_inboth.append(d)
            else:
                raise ValueError("Mismatch on files for %s. Can't use pairedcodonplot option" % sample)
        mapmuts.plot.PlotPairedMutFracs(infiles_inboth, samples_inboth, '%s_codon_types_and_nmuts.pdf' % plotfileprefix)
    print "Script complete."


if __name__ == '__main__':
    main() # run the script
