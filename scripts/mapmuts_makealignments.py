#!python

"""Aligns overlapping paired-end reads to a target gene.

This mapmuts_makealignments.py script is one of the core scripts of the 
mapmuts package.

To run this script from the prompt, first create a text infile of the
format described below. Then simply type mapmuts_makealignments.py
followed by the infile name. For example, if the name is infile.txt,
type::
    mapmuts_makealignments.py infile.txt

If the script is not executable on your system, you can instead type::

    python mapmuts_makealignments.py infile.txt

This script aligns overlapping paired-end sequencing reads to a target
gene sequence. It has been tested on the alignment of 50 nt paired-end
Illumina reads to a roughly 1.5 kb gene, although it should work more
generally. The result of running this script is a variety of output
files containing the actual alignments and summary information about the
alignment process. If pylab / matplotlib are available, summary plots
are also generated. If pdflatex is available, these plots are merged
into a summary PDF document as well.

This script was written by Jesse Bloom.
"""


import sys
import os
import time
import glob
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
    logfile = "%s_makealignments_log.txt" % outfileprefix
    log = open(logfile, 'w')
    try:
        log.write("Beginning execution of mapmuts_makealignments.py"\
                " in directory %s" % (os.getcwd()))
        mapmuts.io.PrintVersions(log)
        log.write("Input data being read from infile %s\n\n" % infilename)
        log.write("Progress being logged to this file, %s\n\n" % logfile)
        log.write("Read the following key/value pairs from infile %s:"\
                % (infilename))
        for (key, value) in d.iteritems():
            log.write("\n%s %s" % (key, value))
        del d['outfileprefix']
        try:
            r1files = mapmuts.io.ParseFileList(d, 'r1files')
        except:
            r1files = mapmuts.io.ParseStringValue(d, 'r1files')
            r1files = glob.glob(r1files)
            if not r1files:
                raise IOError("Can't find any r1files.")
            r1files.sort()
        del d['r1files']
        try:
            r2files = mapmuts.io.ParseFileList(d, 'r2files')
        except:
            r2files = mapmuts.io.ParseStringValue(d, 'r2files')
            r2files = glob.glob(r2files)
            if not r2files:
                raise IOError("Can't find any r2files.")
            r2files.sort()
        del d['r2files']
        gzipped = mapmuts.io.ParseBoolValue(d, 'gzipped')
        del d['gzipped']
        applyfilter = mapmuts.io.ParseBoolValue(d, 'applyfilter')
        del d['applyfilter']
        minq = mapmuts.io.ParseFloatValue(d, 'minq')
        del d['minq']
        fullgene = mapmuts.io.ParseSeqValue(d, 'fullgenefile')
        del d['fullgenefile']
        generange = mapmuts.io.ParseStringValue(d, 'generange')
        try:
            # subtract 1 from first to deal with 0 versus 1 start indexing
            generange = (int(generange.split()[0]) - 1, int(generange.split()[1]))
        except:
            raise ValueError("Invalid generange: %s" % ' '.join(generange))
        del d['generange']
        a1 = mapmuts.io.ParseSeqValue(d, 'a1file')
        del d['a1file']
        a2 = mapmuts.io.ParseSeqValue(d, 'a2file')
        del d['a2file']
        maxn = mapmuts.io.ParseIntValue(d, 'maxn')
        del d['maxn']
        minoverlap = mapmuts.io.ParseIntValue(d, 'minoverlap')
        del d['minoverlap']
        maxrm = mapmuts.io.ParseIntValue(d, 'maxrm')
        del d['maxrm']
        maxa1m = mapmuts.io.ParseIntValue(d, 'maxa1m')
        del d['maxa1m']
        maxa2m = mapmuts.io.ParseIntValue(d, 'maxa2m')
        del d['maxa2m']
        maxgenem = mapmuts.io.ParseIntValue(d, 'maxgenem')
        del d['maxgenem']
        if 'test' == mapmuts.io.ParseStringValue(d, 'upcase'):
            upcase = 'test'
        else:
            upcase = mapmuts.io.ParseBoolValue(d, 'upcase')
            if upcase != True:
                raise ValueError("Invalid value of upcase, must be test or True")
        del d['upcase']
        samplename = mapmuts.io.ParseStringValue(d, 'samplename')
        del d['samplename']
        write_unaligned = mapmuts.io.ParseBoolValue(d, 'write_unaligned')
        del d['write_unaligned']
        if d != {}:
            raise IOError('Input file of %s specifies extra keys: %s'\
                    % (infilename, str(d.keys())))
        log.write('\n\n')
        mapmuts.main.MakeAlignments(r1files, r2files, gzipped, applyfilter,\
                minq, fullgene, generange, a1, a2, maxn, minoverlap, maxrm,\
                maxa1m, maxa2m, maxgenem, outfileprefix, log=log,\
                write_unaligned=write_unaligned)
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
            log.write('\nNow running mapmuts.main.MakeAlignmentsPlots...\n')
            mapmuts.main.MakeAlignmentsPlots(outfileprefix, minq, maxn,\
                    minoverlap, maxrm, maxa1m, maxa2m, maxgenem, latexsummary,\
                    summarytitle=samplename)
            log.write('Completed running mapmuts.main.MakeAlignmentsPlots.\n')
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
