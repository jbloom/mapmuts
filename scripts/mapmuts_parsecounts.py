#!python

"""Parses counts from alignments built by mapmuts_makealignments.py

This mapmuts_parsecounts.py script is one of the core scripts of the 
mapmuts package.

To run this script from the prompt, first create a text infile of the
format described below. Then simply type mapmuts_parsecounts.py
followed by the infile name. For example, if the name is infile.txt,
type::

    mapmuts_parsecounts.py infile.txt

If the script is not executable on your system, you can instead type::

    python mapmuts_parsecounts.py infile.txt

This function parses information from the *_alignments.txt.gz file
created by MakeAlignments. The counts of nucleotide and codon variants
at each position is parsed, and written to two output text files.

The overlapping paired-end read alignments are used to call nucleotide
identities. Nucleotides are called as follows:

* If both reads agree on the identity and neither relevant read position
  is in the exclusion lists (r1exclude or r2exclude), then the
  nucleotide is called at that identity using upper-case nucleotides,
  for example 'A' or 'T'.

* If one read specifies an identity and the other read gives an 'N'
  or is at an excluded position, the identity is called to the 
  lower-case value of the read the specifies an identity. For
  example, a position might be called to 'a' or 't'.

* If the two reads specify different identities, or if both specify
  an N or have the position excluded, a position is called to 'N'

Codons are called to the constituent nucleotides using the 
aforementioned scheme. Codons for which only partial nucleotide
coverage is achieved are not called.
"""


import sys
import os
import time
import mapmuts
import mapmuts.io
import mapmuts.main
import mapmuts.plot
import mapmuts.latex


def main():
    """Main body of script."""
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument"\
                + ' specifying the name of the input file.')
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = mapmuts.io.ParseInfile(open(infilename))
    outfileprefix = mapmuts.io.ParseStringValue(d, 'outfileprefix')
    logfile = "%s_parsecounts_log.txt" % outfileprefix
    log = open(logfile, 'w')
    try:
        log.write("Beginning execution of mapmuts_parsecounts.py"\
                " in directory %s" % (os.getcwd()))
        mapmuts.io.PrintVersions(log)
        log.write("Input data being read from infile %s\n\n" % infilename)
        log.write("Progress being logged to this file, %s\n\n" % logfile)
        log.write("Read the following key/value pairs from infile %s:"\
                % (infilename))
        for (key, value) in d.iteritems():
            log.write("\n%s %s" % (key, value))
        del d['outfileprefix']
        alignmentfile = mapmuts.io.ParseFileList(d, 'alignmentfile')
        if len(alignmentfile) != 1:
            raise IOError("Failed to find exactly one alignmentfile")
        alignmentfile = alignmentfile[0]
        del d['alignmentfile']
        fullgene = mapmuts.io.ParseSeqValue(d, 'fullgenefile')
        del d['fullgenefile']
        generange = mapmuts.io.ParseStringValue(d, 'generange')
        try:
            # subtract 1 from first to deal with 0 versus 1 start indexing
            generange = (int(generange.split()[0]) - 1, int(generange.split()[1]))
        except:
            raise ValueError("Invalid generange: %s" % ' '.join(generange))
        gene = fullgene[generange[0] : generange[1]].upper()
        del d['generange']
        if 'test' == mapmuts.io.ParseStringValue(d, 'upcase'):
            upcase = 'test'
        else:
            upcase = mapmuts.io.ParseBoolValue(d, 'upcase')
        del d['upcase']
        r1exclude = mapmuts.io.ParseStringValue(d, 'r1exclude')
        if r1exclude.strip() == 'None':
            r1exclude = set([])
        else:
            try:
                r1exclude = set([int(x) for x in r1exclude.split()])
            except:
                raise ValueError("Invalid r1exclude: %s" % r1exclude)
        del d['r1exclude']
        r2exclude = mapmuts.io.ParseStringValue(d, 'r2exclude')
        if r2exclude == 'None':
            r2exclude = set([])
        else:
            try:
                r2exclude = set([int(x) for x in r2exclude.split()])
            except:
                raise ValueError("Invalid r2exclude: %s" % r2exclude)
        del d['r2exclude']
        samplename = mapmuts.io.ParseStringValue(d, 'samplename')
        del d['samplename']
        if d != {}:
            raise IOError('Input file of %s specifies extra keys: %s'\
                    % (infilename, str(d.keys())))
        log.write('\n\n')
        mapmuts.main.ParseNTCodonCounts(alignmentfile, outfileprefix, gene,\
                r1exclude, r2exclude, upcase=upcase, log=log)
        mapmuts.io.WriteAAsFromCodonCounts("%s_codoncounts.txt" % outfileprefix,\
                "%s_aacounts.txt" % outfileprefix)
        if mapmuts.plot.PylabAvailable():
            log.write('\nPylab / matplotlib appear to be available, so'\
                    + ' we will make the summary plots.\n')
            latexsummary = bool(mapmuts.latex.PdflatexAvailable())
            if latexsummary:
                log.write('\npdflatex is available, so we will also make'\
                        + ' an overall summary PDF.\n')
            else:
                log.write('\npdflatex is unavailable, so we cannot make'\
                        + ' an overall summary PDF.\n')
            mapmuts.main.ParseNTCodonCountsPlots(outfileprefix, r1exclude,\
                    r2exclude, latexsummary, samplename)
        else:
            log.write('\nPylab / matplotlib are not available, so we '\
                    + 'cannot make any summary plots.\n')
    except:
        for x in sys.exc_info():
            log.write("\n\n%s" % str(x))
        log.write("\n\nPrematurely closing log due to execution error.")
        raise
    finally:
        log.write("\n\nExecution completed at %s." % time.ctime())
        log.close()



if __name__ == '__main__':
    main() # run the script
