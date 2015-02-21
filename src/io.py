"""Performs input/output operations for ``mapmuts`` package.

Written by Jesse Bloom

List of functions
-------------------
`IteratePairedFASTQ` : iterates over FASTQ files for paired-end reads.

`IterateAlignmentFile` : iterates over files written by `main.MakeAlignments`.

`WriteInsertLengths` : writes distribution of insert lengths. 

`ReadInsertLengths` : reads data written by `WriteInsertLengths`.

`WriteAlignmentStatistics` : writes summary of alignment statistics.

`ReadAlignmentStatistics` : reads data written by `WriteAlignmentStatistics`.

`WriteReadMismatches` : writes distribution of mismatches along read.

`ReadReadMismatches` : reads data written by `WriteReadMismatches`.

`WriteNTCounts` : writes counts of nucleotide identities.

`ReadNTCounts` : reads output of `WriteNTCounts`.

`WriteCodonCounts` : writes counts of codon identities.

`ReadCodonCounts` : reads output of `WriteCodonCounts`.

`ReadEnrichmentRatios` : reads enrichment ratios from ``mapmuts_inferenrichment.py``

`ReadEntropyAndEquilFreqs` : reads entropies / equilibrium freqs from 
``mapmuts_inferenrichment.py`` or ``mapmuts_inferpreferences.py``.

`ReadDifferentialPreferences` : reads differential preferences from ``mapmuts_inferdifferentialpreferences.py``

`WriteAAsFromCodonCounts` : reads files created by `WriteCodonCounts`
and writes the amino acid frequencies.

`ParseInfile` : reads key / value pairs from an input file.

`ParseBoolValue` : reads Boolean values from input dictionary.

`ParseIntValue` : reads integer values from input dictionary.

`ParseFloatValue` : reads float values from input dictionary.

`ParseStringValue` : reads string values from input dictionary.

`ParseSeqValue` : reads sequence from file specified by input dictionary.

`ParseFileList` : reads file list specified by input dictionary.

`AllNTs` : lists nucleotide codes accepted by `WriteNTCounts` / `ReadNTCounts`.

`AllCodons` : list codon codes accepted by `WriteCodonCounts` / `ReadCodonCounts`.

`PrintVersions` : prints version numbers of *mapmuts* and associated packages.


Documentation for individual functions
----------------------------------------
Documentation for individual functions is provided in their definitions below.
"""


import os
import sys
import re
import gzip
import tempfile
import time
import platform
import cStringIO
import subprocess
import warnings
import mapmuts
import mapmuts.sequtils
import mapmuts.latex
import mapmuts.weblogo


def PrintVersions(f):
    """Prints information about the versions of ``mapmuts`` and utilized packages.

    *f* is a writable file-like object (for example, could be *sys.stdout* if
    you want to write to standard output).

    This function writes to *f* a summary of the software versions. It
    also writes the current time and platform. For the optional packages,
    if any of the packages are not available then that is indicated. Prints
    information about the following:

        * The current time / date (accessed via *time.asctime()*).

        * The current platform (accessed via *platform.platform()*).

        * The version of Python (accessed via *sys.version*).

        * The version of ``mapmuts`` (accessed via *mapmuts.__version__*).

        * The version of ``numpy`` (accessed via *numpy.__version__*).

        * The version of ``pymc`` (accessed via *pymc.__version__*).

        * The version of ``scipy`` (accessed via *scipy.__version__*).

        * The version of ``matplotlib`` (accessed via *matplotlib.__version__*).

        * The version of ``pyPdf``: unfortunately the ``pyPdf`` package does
          not contain a built-in version string, so we only print whether
          ``pyPdf`` is installed at all.

        * The version of ``pdflatex`` (accessed via *mapmuts.latex.PdflatexAvailable()*).

        * The version of ``weblogo`` (acccessed via *mapmuts.latex.WebLogoAvailable()*).

    Here is an example of the output of this function::

        ****************************************************
        Version information for mapmuts and associated programs.

        Time and date: Thu Dec  5 16:40:30 2013

        Platform: Darwin-10.8.0-x86_64-i386-64bit

        Python version: 2.6.8 (unknown, Apr 14 2012, 04:15:37) 
        [GCC 4.2.1 (Apple Inc. build 5666) (dot 3)]

        mapmuts version: 1.0

        numpy version: 1.7.1

        pymc version: 2.3

        scipy version: 0.12.0

        matplotlib version: 1.2.1

        pyPdf version: pyPdf is available, but no version string accessible

        pdflatex version: pdfTeX 3.1415926-1.40.11-2.2 (TeX Live 2010)
        kpathsea version 6.0.0
        Copyright 2010 Peter Breitenlohner (eTeX)/Han The Thanh (pdfTeX).
        There is NO warranty.  Redistribution of this software is
        covered by the terms of both the pdfTeX copyright and
        the Lesser GNU General Public License.
        For more information about these matters, see the file
        named COPYING and the pdfTeX source.
        Primary author of pdfTeX: Peter Breitenlohner (eTeX)/Han The Thanh (pdfTeX).
        Compiled with libpng 1.2.40; using libpng 1.2.40
        Compiled with zlib 1.2.3; using zlib 1.2.3
        Compiled with xpdf version 3.02pl4

        weblogo version: WebLogo 3.3 (2012-07-02)
        ****************************************************

    """
    f.write('\n\n****************************************************')
    f.write('\nVersion information for mapmuts and associated programs.')
    f.write('\n\nTime and date: %s' % time.asctime())
    f.write('\n\nPlatform: %s' % platform.platform())
    f.write('\n\nPython version: %s' % sys.version)
    f.write('\n\nmapmuts version: %s' % mapmuts.__version__)
    try:
        import numpy
        f.write('\n\nnumpy version: %s' % numpy.__version__)
    except ImportError:
        f.write('\n\nnumpy version: numpy is NOT available')
    try:
        import pymc
        f.write('\n\npymc version: %s' % pymc.__version__)
    except ImportError:
        f.write('\n\npymc version: pymc is NOT available')
    try:
        import scipy
        f.write('\n\nscipy version: %s' % scipy.__version__)
    except ImportError:
        f.write('\n\nscipy version: scipy is NOT available')
    try:
        import matplotlib
        f.write('\n\nmatplotlib version: %s' % matplotlib.__version__)
    except ImportError:
        f.write('\n\nmatplotlib version: matplotlib is NOT available')
    try:
        import pyPdf
        f.write('\n\npyPdf version: pyPdf is available, but no version string accessible')
    except ImportError:
        f.write('\n\npyPdf version: pyPdf is NOT available')
    pdflatexversion = mapmuts.latex.PdflatexAvailable()
    if not pdflatexversion:
        pdflatexversion = 'pdflatex is NOT available'
    f.write('\n\npdflatex version: %s' % pdflatexversion.strip())
    weblogoversion = mapmuts.weblogo.WebLogoAvailable()
    if not weblogoversion:
        weblogoversion = 'weblogo is NOT available'
    f.write('\n\nweblogo version: %s' % weblogoversion.strip())
    f.write('\n****************************************************\n\n')


def IteratePairedFASTQ(r1files, r2files, gzipped, applyfilter, usegzip=False):
    """Iterates over FASTQ file pairs for paired-end sequencing reads.

    This function has been tested for its ability to read FASTQ files
    generated by the Illumina Casava 1.8 pipeline for paired-end reads.

    CALLING VARIABLES:

    `r1files` : list of file name(s) containing the R1 reads.

    `r2files` : list of file name(s) containing the R2 reads. The files
    specified here must specify the exact same number of reads as
    those in `r1files`, and they must be in the same sequence (i.e.
    the first R1 read found going through `r1files` in order must
    match the first R2 read found going through `r2files` in order,
    etc).

    `gzipped` : Boolean switch specifying whether the FASTQ files
    specified by `r1files` and `r2files` are gzipped. If set to `True`,
    they are gzipped. This iterator will then read them without
    unzipping the files.

    `applyfilter` : Boolean switch specifying whether we remove read
    pairs in which one or more of the reads failed the Illumina
    chastity filter. If `True`, the iterator simply returns `False`
    each time it encounters a read pair where the chastity filter is
    failed. The chastity filter is indicated by Y or N in the
    sequence header -- values of Y indicate a failed filter.

    `usegzip` : Boolean switch that is meaningful only if `gzipped`
    is `True`. This specifies that we use the Python ``gzip`` module to 
    unzip the reads. Otherwise we use the system ``gzip -cd command``,
    which may be faster. You should definitely specify
    `usegzip` to `True` if you are only reading a few lines, as with
    `usegzip` set to `False` the entire file is unzipped before anything happens.
    `False` by default.

    RESULT OF THIS FUNCTION:

    The function will iterate over all read pairs in the indicated
    files. For each iteration, it will return:

    * The return variable is `False` if ``applyfilter == True`` and one of the
      reads in the pair failed the Illumina chastity filter. 
    
    * Otherwise the
      return variable is the following tuple:
      `(name, r1, r2, q1avg, q2avg)` where:

      `name` : name of the read pair (a string)

      `r1` : sequence of the first read R1 (a string)

      `r2` : sequence of the second read R2 (a string)

      `q1avg` : the mean Q-score for the R1 read (a float)

      `q2avg` : the mean Q-score for the R2 read (a float)
    """
    hm = re.compile('^\@(?P<name>[\-\w]+\:\d+\:[\w\-]+\:\d+\:\d+\:\d+\:\d+)'\
            + ' (?P<direction>[1,2])\:(?P<filtered>[Y,N])\:\d+\:'\
            + '\w+\\n$')
    if len(r1files) != len(r2files):
        raise IOError('r1files and r2files are of different length.')
    nfiles = len(r1files)
    ifile = 0
    while ifile < nfiles:
        if not os.path.isfile(r1files[ifile]):
            raise IOError("Can't find file %s" % r1files[ifile])
        if not os.path.isfile(r2files[ifile]):
            raise IOError("Can't find file %s" % r2files[ifile])
        if gzipped:
            # rather than using gzip module, we try using gzip -cd and then
            # read output. This is supposed to give better performance:
            # http://codebright.wordpress.com/2011/03/25/139/
            # unfortunately piping to subprocess.PIPE does not work
            # for the reason described in http://www.macaronikazoo.com/?p=607
            # so we have to use the following less efficient procedure of
            # first writing everything to a temporary file.
            try:
                if usegzip:
                    raise OSError # to make us go to except statement
                f1 = tempfile.TemporaryFile()
                p1 = subprocess.Popen(['gzip', '-cd', r1files[ifile]],\
                    stdout=f1.fileno())
                p1.wait()
                f1.seek(0)
                f2 = tempfile.TemporaryFile()
                p2 = subprocess.Popen(['gzip', '-cd', r2files[ifile]],\
                    stdout=f2.fileno())
                p2.wait()
                f2.seek(0)
            except OSError:
                # cannot use gzip -cd, try using gzip module
                if not usegzip:
                    warnings.warn("Failed to use gzip -cd, will have to use "+\
                            "the slower Python gzip module instead.",\
                            RuntimeWarning)
                f1 = gzip.open(r1files[ifile])
                f2 = gzip.open(r2files[ifile])
        else:
            f1 = open(r1files[ifile])
            f2 = open(r2files[ifile])
        try:
            while True:
                (h1, s1, x1, q1) = (f1.next(), f1.next(), f1.next(),\
                        f1.next())
                (h2, s2, x2, q2) = (f2.next(), f2.next(), f2.next(),\
                        f2.next())
                h1m = hm.search(h1)
                h2m = hm.search(h2)
                if not h1m:
                    raise IOError("Failed to match head:\n%s" % h1)
                if not h2m:
                    raise IOError("Failed to match head:\n%s" % h2)
                (n1, d1, cf1) = (h1m.group('name'),\
                        h1m.group('direction'), h1m.group('filtered'))
                (n2, d2, cf2) = (h2m.group('name'),\
                        h2m.group('direction'), h2m.group('filtered'))
                if n1 != n2:
                    raise IOError("Name mismatch:\n%s\n%s" % (n1, n2))
                if not (d1 == '1' and d2 == '2'):
                    raise IOError("Direction mismatch:\n%s\n%s" % (n1, n2))
                if applyfilter and (cf1 == 'Y' or cf2 == 'Y'):
                    yield False # failed chastity filter
                s1 = s1.strip()
                s2 = s2.strip()
                q1avg = mapmuts.sequtils.MeanQValueFromString(q1.strip())
                q2avg = mapmuts.sequtils.MeanQValueFromString(q2.strip())
                yield (n1, s1, s2, q1avg, q2avg)
        except StopIteration:
            f1.close()
            f2.close()
            ifile += 1


def IterateAlignmentFile(filename, gzipped=True, usegzip=False):
    """Iterates over alignment file written by main.MakeAlignments

    `filename` : name of a file of the type written by `main.MakeAlignments`
    as the ``*_alignments.txt.gz`` file. Each entry
    in these alignment files consists of five lines, giving the
    read name, gene, read 1, read 2, blank line, as in::

            @DH1DQQN1:241:C1433ACXX:2:1101:1623:2365
            14 ATAGATACTAGAT 26
             1 .....C.....T. 13
            13 .....C.......  1

    `gzipped` : Boolean switch specifying whether the alignment
    file is gzipped. In general, this will be `True` since
    `main.MakeAlignments` writes gzipped files. Default
    value is `True`.

    `usegzip` : a Boolean switch that is meaningful only if `gzipped` 
    is `True`. This specifies that we use the Python ``gzip`` module to 
    unzip the reads. Otherwise we use the system ``gzip -cd`` command,
    which may be faster. You should definitely specify
    `usegzip` to `True` if you are only reading a few lines, as with
    `usegzip` as `False` the entire file is unzipped before anything happens.
    `False` by default.

    The function will iterate over all alignment entries in
    filename. Each iteration returns the following tuple:
    `(gstart, gend, gseq, r1start, r1end, r1seq, r2start, r2end, r2seq)`
    The tuple entries have the following meanings:

        `gstart` : index of first gene nt; 1, 2, ... numbering.

        `gend` : index last gene nt; 1, 2, ... numbering.

        `gseq` : aligned gene nucleotides.

        `r1start` -> index of first R1 nt; 1, 2, ... numbering.

        `r1end` : index of last R1 nt; 1, 2, ... numbering.

        `r1seq` : aligned R1 nucleotides, with dots at identities.
    
        `r2start` : index of first R2 nt; 1, 2, ... numbering.

        `r2end` : index of last R2 nt; 1, 2, ... numbering.
    
        `r2seq` : aligned R2 nucleotides, with dots at identities.

    Note that the `r1` and `r2` indices are for the orientation in
    which the read aligns, so either `r1start <= r1end` or
    `r2start <= r2end`.

    For example, the alignment shown above under the filename
    option description would return::

        (14, 26, 'ATAGATACTAGAT', 1, 13, '.....C.....T.', 13, 1, '.....C.......')
    """
    if not os.path.isfile(filename):
        raise IOError("Can't find file %s" % filename)
    if gzipped:
        # rather than using gzip module, we try using gzip -cd and then
        # read output. This is supposed to give better performance:
        # http://codebright.wordpress.com/2011/03/25/139/
        # unfortunately piping to subprocess.PIPE does not work
        # for the reason described in http://www.macaronikazoo.com/?p=607
        # so we have to use the following less efficient procedure of
        # first writing everything to a temporary file.
        try:
            if usegzip:
                raise OSError # to make us go to except statement
            f = tempfile.TemporaryFile()
            p = subprocess.Popen(['gzip', '-cd', filename],\
                    stdout=f.fileno())
            p.wait()
            f.seek(0)
        except OSError:
            # cannot use gzip -cd, try using gzip module
            if not usegzip:
                warnings.warn("Failed to use gzip -cd, will have to use "+\
                        "the slower Python gzip module instead.",\
                        RuntimeWarning)
            f = gzip.open(filename)
    else:
        f = open(filename)
    try:
        while True:
            (n, g, r1, r2, b) = (f.next(), f.next(), f.next(),\
                        f.next(), f.next())
            try:
                entries = g.split()
                (gstart, gseq, gend) = (int(entries[0]), entries[1],\
                        int(entries[2]))
            except ValueError:
                raise IOError("Can't parse gene alignment:\n%s" % g)
            try:
                entries = r1.split()
                (r1start, r1seq, r1end) = (int(entries[0]), entries[1],\
                        int(entries[2]))
            except ValueError:
                raise IOError("Can't parse R1 alignment:\n%s" % r1)
            try:
                entries = r2.split()
                (r2start, r2seq, r2end) = (int(entries[0]), entries[1],\
                        int(entries[2]))
            except ValueError:
                raise IOError("Can't parse R2 alignment:\n%s" % r2)
            yield (gstart, gend, gseq, r1start, r1end, r1seq, r2start,\
                    r2end, r2seq)
    except StopIteration:
        f.close()


def WriteAlignmentStatistics(s, f):
    """Writes summary of statistics from `MakeAlignment` attempt.

    This function writes a summary of the attempted / successful read
    alignment statistics to a writeable file-like object.

    `s` : a dictionary holding the statistics in key / value format, with
    string keys and integer values. The specific keys are described
    below.

    `f` : a writeable file-like object to which we write the summary.

    `s` should have the following keys. These constitute the data written
    to `f`. Their meanings are explained in the example below:

        'nread'

        'nfiltered'

        'nlowq'

        'nexcessn'

        'nattempted'

        'npaired'

        'nalignedfullgene'

        'nalignedgene'

    Here is example output, comment lines beginning with # introduce
    each variable, which is then named and followed by its value.

    >>> s = {'nread':110, 'nfiltered':10, 'nlowq':5, 'nexcessn':5, 'nattempted':90,
    ...    'npaired':65, 'nalignedfullgene':55, 'nalignedgene':50}
    >>> WriteAlignmentStatistics(s, sys.stdout)
    # nread -> total number of read pairs read from the FASTQ files.
    nread 110
    # nfiltered -> number of read pairs removed by the Illumina filter.
    nfiltered 10
    # nlowq -> number of read pairs removed due to low mean Q-score
    nlowq 5
    # nexcessn -> number of read pairs removed due to excess N/n nucleotides
    nexcessn 5
    # nattempted -> number of read pairs we attempted to align
    nattempted 90
    # npaired -> number of read pairs where R1 and R2 were successfully
    # overlapped
    npaired 65
    # nalignedfullgene -> number of read pairs where the overlapped read
    # pairs were successfully aligned to fullgene (the template
    # amplicon).
    nalignedfullgene 55
    # nalignedgene -> number of read pairs which overlapped with the
    # gene of interest (region of fullgene specified by generange).
    nalignedgene 50

    """
    f.write('# nread -> total number of read pairs read from'\
            + ' the FASTQ files.\nnread %d\n' % s['nread'])
    f.write('# nfiltered -> number of read pairs removed by the Illumina'\
            + ' filter.\nnfiltered %d\n' % s['nfiltered'])
    f.write('# nlowq -> number of read pairs removed due to low mean '\
            + 'Q-score\nnlowq %d\n' % s['nlowq'])
    f.write('# nexcessn -> number of read pairs removed due to excess N/n'\
            + ' nucleotides\nnexcessn %d\n' % s['nexcessn'])
    f.write('# nattempted -> number of read pairs we attempted to'\
            + ' align\nnattempted %d\n' % s['nattempted'])
    f.write('# npaired -> number of read pairs where R1 and R2 were'\
            + ' successfully\n# overlapped\nnpaired %d\n' %\
            s['npaired'])
    f.write('# nalignedfullgene -> number of read pairs where the'\
            + ' overlapped read\n# pairs were successfully aligned'\
            + ' to fullgene (the template\n#'\
            + ' amplicon).\nnalignedfullgene %d\n' % s['nalignedfullgene'])
    f.write('# nalignedgene -> number of read pairs which overlapped'\
            + ' with the\n# gene of interest (region of fullgene'\
            + ' specified by generange).\nnalignedgene %d' %\
            s['nalignedgene'])


def ReadAlignmentStatistics(f):
    """Reads alignment statistics written by `WriteAlignmentStatistics`.

    `f` : a readable file-like object, opened to the position
    where output was written by `WriteAlignmentStatistics`.

    This function returns a dictionary `s` containing the statistics used
    to call `WriteAlignmentStatistics` to generate the output in `f`.

    Here is an example showing how `ReadAlignmentStatistics` recovers the
    output written by `WriteAlignmentStatistics`.

    >>> s = {'nread':110, 'nexcessn':5, 'nfiltered':10, 'nlowq':5, 'nattempted':90, 
    ...    'npaired':65, 'nalignedfullgene':55, 'nalignedgene':50}
    >>> f = cStringIO.StringIO()
    >>> WriteAlignmentStatistics(s, f)
    >>> f.seek(0)
    >>> s == ReadAlignmentStatistics(f)
    True
    """
    if not hasattr(f, 'read'):
        raise ValueError("f is not readable file-like object")
    lines = [line for line in f.readlines() if line[0] != '#']
    s = {}
    for line in lines:
        entries = line.split()
        if len(entries) != 2:
            raise IOError("Invalid line:\n%s" % line)
        s[entries[0].strip()] = int(entries[1])
    return s


def ReadDifferentialPreferences(infile):
    """"Reads differential preferences written by ``mapmuts_inferdifferentialpreferences.py``.

    The single calling argument *infile* should specify a ``differentialpreferences_selection_*.txt`` file
    written by ``mapmuts_inferdifferentialpreferences.py``.

    These files have the following format::

        #SITE   WT_AA   RMS_dPI dPI_A   dPI_C   dPI_D   dPI_E   dPI_F   dPI_G   dPI_H   dPI_I   dPI_K   dPI_L   dPI_M   dPI_N   dPI_P   dPI_Q   dPI_R   dPI_S   dPI_T   dPI_V   dPI_W   dPI_Y   dPI_*
        1   M   0.033508    -0.0041802  -0.00160924 -0.00102518 -0.00293463 -0.00112235 -0.000467141    -0.00262952 -0.0032969  -0.00679166 -0.00304567 0.0300719   0.00925138  -0.000442459    -8.11477e-05    -0.00129343 -0.00128639 0.000560046 -0.00243062 -0.00313798 -0.002528   -0.00158079
        2   K   0.0630635   -0.0159228  0.0128652   -0.00367557 -0.0167789  0.005023    -0.000661766    0.0260142   -0.0040351  0.0025629   0.00523521  0.00449474  0.00192782  0.0189381   0.0171668   -0.00881324 -0.00971018 0.00748405  0.0136674   -0.00289523 -0.0288229  -0.024063

    The first column gives the residue number, the second column gives the wildtype amino acid, the third
    column gives the root-mean-square (RMS) differential preference, and the remaining columns give the values
    for each of the 20 amino acids. The columns can be delimited by a single space or a tab. There **may** be
    a final column ``dPI_*`` giving the differential preference for a stop codon. That column can be either
    present or absent -- either one is OK.

    The returned variable is a dictionary keyed by residue numbers, and with
    the values for each residue number being dictionaries keyed by all of the
    other 22 or 23 column labels ('WT_AA', 'RMS_dPI', 'dPI_A', 'dPI_C', ...)
    and possibly 'dPI_*' if that is present in the header.
    """
    if not os.path.isfile(infile):
        raise IOError("Cannot find infile %s" % infile)
    lines = open(infile).readlines()
    headmatch = re.compile('^#SITE\sWT_AA\sRMS_dPI\sdPI_A\sdPI_C\sdPI_D\sdPI_E\sdPI_F\sdPI_G\sdPI_H\sdPI_I\sdPI_K\sdPI_L\sdPI_M\sdPI_N\sdPI_P\sdPI_Q\sdPI_R\sdPI_S\sdPI_T\sdPI_V\sdPI_W\sdPI_Y(?P<stop>(\sdPI_\*){0,1})\n$')
    m = headmatch.search(lines[0])
    if not m:
        raise ValueError("Invalid header of:\n%s" % lines[0])
    has_stop = bool(m.group('stop'))
    keys = lines[0].split()[1 : ]
    nkeys = len(keys)
    assert (not has_stop and nkeys == 22) or (has_stop and nkeys == 23), "Invalid header length, not the expected number of keys (found %d)." % nkeys
    d = {}
    for line in lines[1 : ]:
        if line and not line.isspace():
            entries = line.split()
            if len(entries) != nkeys + 1:
                raise ValueError("Invalid line of:\n%s" % line)
            site = int(entries[0])
            if site in d:
                raise ValueError("Duplicate entry for site %d" % site)
            site_d = {}
            for (key, value) in zip(keys, entries[1 : ]):
                if key == 'WT_AA':
                    site_d[key] = value
                else:
                    site_d[key] = float(value)
            d[site] = site_d
    return d


def ReadEntropyAndEquilFreqs(infile):
    """Reads entropies and preferences written by ``mapmuts_inferenrichment.py`` and ``mapmuts_inferpreferences.py``.

    The single calling argument *infile* should specify a 
    ``*_equilibriumfreqs.txt`` file written by ``mapmuts_inferenrichment.py``
    or an ``*_equilibriumpreferences.txt`` written by ``mapmuts_inferpreferences.py``.

    These files have the following format::

        #SITE WT_AA SITE_ENTROPY PI_A PI_C PI_D PI_E PI_F PI_G PI_H PI_I PI_K PI_L PI_M PI_N PI_P PI_Q PI_R PI_S PI_T PI_V PI_W PI_Y
        1   M   3.952438    0.027276    0.057787    0.061477    0.0589700.058843    0.023888    0.050938    0.017797    0.062260    0.019259    0.216150    0.070343    0.029413    0.052597    0.018734    0.013009    0.025890    0.020677    0.058900    0.055793
        2   A   2.362402    0.537118    0.012475    0.013213    0.0258260.014365    0.008897    0.005926    0.009136    0.013902    0.003019    0.191576    0.010007    0.001631    0.007791    0.003013    0.004130    0.005614    0.099395    0.020214    0.012752
        3   S   2.745295    0.024799    0.002901    0.002167    0.0071370.311824    0.002504    0.022821    0.008849    0.004478    0.059771    0.015584    0.022394    0.005062    0.004390    0.006270    0.091925    0.014399    0.003761    0.048111    0.340854

    The first column gives the residue number, the second gives the wildtype
    residue, the third the site entropy, and the remainder the equilibrium
    frequencies of each of the 20 amino acids. The columns can be delimited
    by a single space or single tab.

    There **may** be a final column ``PI_*`` giving the frequency for stop codons.
    That column can either be present or absent -- either one is OK.

    The returned variable is a dictionary keyed by residue numbers, and with
    the values for each residue number being dictionaries keyed by all of the
    other 22 or 23 column labels ('WT_AA', 'SITE_ENTROPY', 'PI_A', 'PI_C', ...)
    and possibly 'PI_*' if that is present in the header.
    """
    if not os.path.isfile(infile):
        raise IOError("Cannot find infile %s" % infile)
    lines = open(infile).readlines()
    headmatch = re.compile('^#\s{0,1}(SITE|POSITION)\sWT(_AA){0,1}\sSITE_ENTROPY\sPI_A\sPI_C\sPI_D\sPI_E\sPI_F\sPI_G\sPI_H\sPI_I\sPI_K\sPI_L\sPI_M\sPI_N\sPI_P\sPI_Q\sPI_R\sPI_S\sPI_T\sPI_V\sPI_W\sPI_Y(?P<stop>(\sPI_\*){0,1})\s')
    m = headmatch.search(lines[0])
    if not m:
        raise ValueError("Invalid header of:\n%s" % lines[0])
    has_stop = bool(m.group('stop'))
    keys = lines[0][1 : ].split()[1 : ]
    nkeys = len(keys)
    assert (not has_stop and nkeys >= 22) or (has_stop and nkeys >= 23), "Invalid header length, not the expected number of keys (found %d)." % nkeys
    d = {}
    for line in lines[1 : ]:
        if line and not line.isspace():
            entries = line.split()
            if len(entries) != nkeys + 1:
                raise ValueError("Invalid line of:\n%s\nExpected %d entries, but got %d" % (line, nkeys + 1, len(entries)))
            site = int(entries[0])
            if site in d:
                raise ValueError("Duplicate entry for site %d" % site)
            site_d = {}
            for (key, value) in zip(keys, entries[1 : ]):
                if key == 'WT_AA' or key == 'WT':
                    site_d[key] = value
                elif '_95' != key[-3 : ]:
                    site_d[key] = float(value)
            d[site] = site_d
    return d


def ReadEnrichmentRatios(infile):
    """Reads enrichment ratios written by ``mapmuts_inferenrichment.py``.

    The single calling argument *infile* should specify a 
    ``*_enrichmentratios.txt`` file written by ``mapmuts_inferenrichment.py``.
    Specifically, these files have the following format::

        #MUTATION PHI PHI_HPD95_LOW PHI_HPD95_HIGH DIRECT_RATIO MUTDNA_COUNTS
        A22N    0.030913    0.013112    0.049075    0.021807    364
        A22Q    0.070328    0.043265    0.097744    0.062814    408
        A22P    0.133431    0.103709    0.160012    0.203831    1128
        A22S    1.936427    1.787140    2.098920    1.740187    1065
        A22R    0.725581    0.653297    0.802220    0.747375    827
        A22T    0.175550    0.141060    0.210771    0.243729    1140
        A22W    0.141486    0.077679    0.211875    0.120190    121

    where the first column gives the mutation, the second gives the enrichment
    ratio (always a number > 0),
    the third and fourth give the 95% highest posterior density for the
    enrichment ratio, the fifth column gives the directly estimated
    enrichment ratio denotes the enrichment ratio (a number that can range
    from minus infinity to infinity inclusive), and the fifth gives the number
    of counts in the mutant library (a number >= 0).

    The returned variable is a dictionary keyed by mutation name. For each
    key, the value is itself a dictionary with the string keys *PHI*,
    *PHI_HPD95_LOW*, *PHI_HPD95_HIGH*, *DIRECT_RATIO*, and *MUTDNA_COUNTS*
    with the values matching the numbers specified for that mutation in
    *infile*. If a mutation is found multiple times in *infile*, and exception
    is raised. For *DIRECT_RATIO*, if the values are ``-inf`` or ``inf`` then
    the value is *-float('inf')* or *float('inf')*.
    """
    if not os.path.isfile(infile):
        raise IOError("Failed to find infile of %s" % infile)
    lines = open(infile).readlines()
    if lines[0].replace('\t', ' ').strip() != '#MUTATION PHI PHI_HPD95_LOW PHI_HPD95_HIGH DIRECT_RATIO MUTDNA_COUNTS':
        raise ValueError("The first line does not match the expected header. This line is:\n%s" % lines[0])
    ratio_d = {}
    for line in lines[1 : ]:
        if line and not line.isspace():
            entries = line.split()
            if len(entries) != 6:
                raise ValueError("Line does not have six entries:\n%s" % line)
            mut = entries[0].strip()
            if mut in ratio_d:
                raise ValueError("Duplicate mutation name:\n%s" % mut)
            ratio_d[mut] = {}
            for (x, key) in zip(entries[1 : ], ['PHI', 'PHI_HPD95_LOW', 'PHI_HPD95_HIGH', 'DIRECT_RATIO', 'MUTDNA_COUNTS']):
                if key == 'MUTDNA_COUNTS':
                    try:
                        x = int(x)
                    except ValueError:
                        raise ValueError("Last entry not an integer in:\n%s" % line)
                else:
                    try:
                        x = float(x)
                    except ValueError:
                        raise ValueError("Didn't find a number for %s in line:\n%s" % (key, line))
                ratio_d[mut][key] = x
    assert len(ratio_d) == len(lines) - 1
    return ratio_d


def WriteInsertLengths(d, f):
    """Writes distribution of paired read insert lengths.
    
    This function writes a text summary of the distribution of different
    insert lengths found between successfully overlapped paired-end
    reads. The insert length is the size of the template fragment
    between the two read adaptors.

    `d` : dictionary holding the distribution of insert lengths. It is
    keyed by integers giving the insert lengths. For each length,
    the value is an integer giving the number of times inserts of
    this length are found. This method will raise a `ValueError` if `d`
    is empty or records no statistics.

    `f` : writeable file-like object to which we write the summary.

    Here is example output, comment lines beginning with # introduce
    head each column.

    >>> d = {20:3, 21:5, 22:7, 23:10, 24:5, 25:2, 27:1}
    >>> WriteInsertLengths(d, sys.stdout)
    # Distribution of insert lengths
    # Length Number Fraction
    20 3 0.0909
    21 5 0.1515
    22 7 0.2121
    23 10 0.3030
    24 5 0.1515
    25 2 0.0606
    26 0 0.0000
    27 1 0.0303
    """
    tot = 0
    minl = maxl = None
    for (x, y) in d.iteritems():
        if minl == None:
            minl = maxl = x
        else:
            minl = min(x, minl)
            maxl = max(x, maxl)
        tot += y
    if not tot:
        raise ValueError("No insert lengths specified:\n%s" % str(d))
    f.write('# Distribution of insert lengths\n')
    f.write('# Length Number Fraction\n')
    for x in range(minl, maxl + 1):
        try:
            y = d[x]
        except KeyError:
            y = 0
        f.write("%d %d %.4f\n" % (x, y, y / float(tot)))


def ReadInsertLengths(f):
    """Reads alignment statistics written by `WriteInsertLengths`.

    `f` : a readable file-like object, opened to the position
    where output was written by `WriteInsertLengths`.

    This function returns a dictionary `d` containing the statistics used
    to call `WriteInsertLengths` to generate the output in `f`.

    Here is an example showing how `ReadInsertLengths` recovers the
    output written by `WriteInsertLengths`:

    >>> d = {20:3, 21:5, 22:7, 23:10, 24:5, 25:2, 27:1}
    >>> f = cStringIO.StringIO()
    >>> WriteInsertLengths(d, f)
    >>> f.seek(0)
    >>> d == ReadInsertLengths(f)
    True
    """
    lines = [line for line in f.readlines() if line[0] != '#']
    d = {}
    for line in lines:
        entries = line.split()
        if len(entries) != 3:
            raise ValueError("Invalid number of entries:\n%s" % line)
        (x, y) = (int(entries[0]), int(entries[1]))
        if y:
            d[x] = y
    return d


def WriteReadMismatches(rlengths, rmismatches, f):
    """Writes distribution of mismatches along the read length.
    
    This function writes a text summary of the distribution of
    mismatches with a target gene as a function of position along the
    read. 

    `rlengths` : a dictionary keyed by the length of the reads that are
    being analyzed for mismatches, and with values being the number
    of reads of this length. For example, ``{25:2, 26:3}`` indicates
    that we have 2 reads of length 25 and 3 reads of length 26. All
    reads are assumed to begin at the first position and extend up
    to the indicated length.

    `rmismatches` : a dictionary giving the position of mismatches along
    the read. The key is the read position in 1, 2, ... numbering.
    The values are the number of mismatches found at that position
    in the read.

    `f` : a writeable file-like object to which we write the summary.

    Here is example output, comment lines beginning with # introduce
    head each column:

    >>> rlengths = {10:2, 11:3, 12:2}
    >>> rmismatches = {1:2, 3:1, 5:1, 8:2, 11:1, 12:1}
    >>> WriteReadMismatches(rlengths, rmismatches, sys.stdout)
    # Distribution of mismatches along the read length.
    # Read positions numbered 1, 2, ... from read start.
    # ReadPosition NReads NMismatches FracMismatches
    1 7 2 0.285714
    2 7 0 0.000000
    3 7 1 0.142857
    4 7 0 0.000000
    5 7 1 0.142857
    6 7 0 0.000000
    7 7 0 0.000000
    8 7 2 0.285714
    9 7 0 0.000000
    10 7 0 0.000000
    11 5 1 0.200000
    12 2 1 0.500000

    """
    maxl = max(rlengths.keys())
    nreads = dict([(x, 0) for x in range(1, maxl + 1)])
    for (x, y) in rlengths.iteritems():
        for i in range(1, x + 1):
            nreads[i] += y
    f.write("# Distribution of mismatches along the read length.\n")
    f.write("# Read positions numbered 1, 2, ... from read start.\n")
    f.write("# ReadPosition NReads NMismatches FracMismatches\n")
    for i in range(1, maxl + 1):
        try:
            y = rmismatches[i]
        except KeyError:
            y = 0
        x = nreads[i]
        if x and y:
            f.write("%d %d %d %.6f\n" % (i, x, y, y / float(x)))
        elif x and not y:
            f.write("%d %d %d %.6f\n" % (i, x, y, 0))
        else:
            raise ValueError("Mismatch at position %d, no reads." % i)


def ReadReadMismatches(f):
    """Reads read mismatches written by `WriteReadMismatches`.

    `f` : should be a readable file-like object, opened to the position
    where output was written by `WriteReadMismatches`.

    This function returns a dictionary `d` containing information about
    the distribution of read mismatches along the read from the data
    written to `f` by `WriteReadMismatches`. The dictionary is keyed by read
    positions (in 1, 2, ... numbering) and the values are the fraction
    of mismatches at that position among all reads. 
   
    Here is an example:

    >>> rlengths = {10:2, 11:3, 12:2}
    >>> rmismatches = {1:2, 3:1, 5:1, 8:2, 11:1, 12:1}
    >>> f = cStringIO.StringIO()
    >>> WriteReadMismatches(rlengths, rmismatches, f)
    >>> f.seek(0)
    >>> d = {1:0.285714, 2:0.0, 3:0.142857, 4:0.0, 5:0.142857,
    ...      6:0.0, 7:0.0, 8:0.285714, 9:0.0, 10:0.0, 11:0.2, 12:0.5}
    >>> print d == ReadReadMismatches(f)
    True
    """
    if not hasattr(f, 'read'):
        raise ValueError("f is not readable file-like object")
    lines = [line for line in f.readlines() if line[0] != '#']
    d = {}
    for line in lines:
        entries = line.split()
        if len(entries) != 4:
            raise ValueError("Invalid line:\n%s" % line)
        d[int(entries[0])] = float(entries[3])
    return d


def WriteNTCounts(nt_counts, f):
    """Writes counts of nucleotide identities.
   
    This function takes a dictionary that specifies the number of
    times each nucleotide identity was observed at a position. This
    is used to write a text output file.

    `nt_counts` : a dictionary keyed by nucleotide number. These integer
    keys should range from 1 to `max(nt_counts.keys())` inclusive, and
    there must be a key for every number within this range. For key
    `nti`, `nt_counts[nti]` is another dictionary specifying the counts
    for this specific nucleotide. This nucleotide dictionary must
    have the following string keys:

        'WT' : value should be a letter indicating wildtype nucleotide
        identity (such as 'A', 'T', 'C', or 'G')

        'COUNTS' -> integer total number of counts for this nucleotide

        All nucleotide codes in `AllNTs()`. For these codes, upper-case
        letters (A, T, C, G) specify identities called by both reads.
        Lower-case letter (a, t, c, g) specify identities called by
        just one read, with the other read ambiguous or excluded at
        the relevant position. N specifies identities where the 
        two reads disagree, both are ambiguous, or both are excluded.
        If any of the required keys are missing, raises a `ValueError`.

    `f` : writable file like object to which we write the counts.

    Here is example output, comment lines beginning with # introduce
    each variable, which is then named and followed by its value. The
    order of the nucleotide entries matches the order returned by `AllNTs()`:

    >>> nt_counts = {1:{'WT':'A', 'COUNTS':7, 'A':3, 'T':1, 'C':0, 'G':0,
    ...      'a':1, 't':1, 'c':0, 'g':0, 'N':1}, 2:{'WT':'C', 'COUNTS':8,
    ...       'A':0, 'T':0, 'C':6, 'G':0, 'a':1, 't':0, 'c':1, 'g':0, 'N':0}}
    >>> WriteNTCounts(nt_counts, sys.stdout)
    # nucleotide counts
    # upper case indicates nucleotide called by both reads
    # lower case indicates nucleotide called by one read with other ambiguous
    # N indicates disagreement between reads, or both ambiguous
    # POSITION WT COUNTS A T C G a t c g N
    1 A 7 3 1 0 0 1 1 0 0 1
    2 C 8 0 0 6 0 1 0 1 0 0

    """
    if not nt_counts:
        raise ValueError("empty nt_counts")
    maxnt = max(nt_counts.keys())
    keys = ['WT', 'COUNTS'] + AllNTs()
    f.write('# nucleotide counts\n')
    f.write('# upper case indicates nucleotide called by both reads\n')
    f.write('# lower case indicates nucleotide called by one read with other ambiguous\n')
    f.write('# N indicates disagreement between reads, or both ambiguous\n')
    f.write('# POSITION')
    for key in keys:
        f.write(' %s' % key)
    for i in range(1, maxnt + 1):
        f.write("\n%d" % i)
        if i not in nt_counts:
            raise ValueError("Didn't find key for %d in nt_counts" % i)
        if 'WT' not in nt_counts[i]:
            raise ValueError("Didn't find key for WT in nt_counts[%d]" % i)
        f.write(" %s" % nt_counts[i]['WT'])
        for key in keys[1 : ]:
            f.write(" %d" % nt_counts[i][key])


def ReadNTCounts(f):
    """Reads counts of nucleotide identities written by `WriteNTCounts`.
  
    This function reads the files written by `WriteNTCounts` to
    reconstruct the calling dictionary specifying the number
    of time each nucleotide identity is observed at a position.

    `f` : a readable file-like object of the type written by 
    `WriteNTCounts`, open to the position where the nucleotide
    counts begin.

    The returned variable is a dictionary that should be identical in
    its contents to the `nt_counts` dictionary used to call `WriteNTCounts`.

    Here is an example:

    >>> nt_counts = {1:{'WT':'A', 'COUNTS':7, 'A':3, 'T':1, 'C':0, 'G':0,
    ...      'a':1, 't':1, 'c':0, 'g':0, 'N':1}, 2:{'WT':'C', 'COUNTS':8,
    ...       'A':0, 'T':0, 'C':6, 'G':0, 'a':1, 't':0, 'c':1, 'g':0, 'N':0}}
    >>> f = cStringIO.StringIO()
    >>> WriteNTCounts(nt_counts, f)
    >>> f.seek(0)
    >>> d = ReadNTCounts(f)
    >>> d == nt_counts
    True

    """
    if not hasattr(f, 'read'):
        raise ValueError("f is not readable file-like object")
    nt_counts = {}
    allnts = AllNTs()
    keys = ['WT', 'COUNTS'] + allnts
    for line in f:
        if line[0] == '#' or line.isspace() or not line:
            continue # comment line
        try:
            entries = line.split()
            if len(entries) != len(keys) + 1:
                raise ValueError
        except ValueError:
            raise ValueError("Couldn't parse line:\n%s" % line)
        if entries[1] not in allnts:
            raise ValueError("Invalid WT of %s" % entries[1])
        position = int(entries[0])
        if position in nt_counts:
            raise ValueError("Duplicate positions of %d" % position)
        id = {keys[0]:entries[1]}
        for (key, value) in zip(keys[1 : ], entries[2 : ]):
            id[key] = int(value)
        nt_counts[position] = id
    if not nt_counts:
        raise IOError("Failed to read any content.")
    maxi = max(nt_counts.keys())
    if maxi < 1:
        raise ValueError("Maximum position < 1")
    for i in range(1, maxi + 1):
        if i not in nt_counts:
            raise ValueError("No information for position %d" % i)
    return nt_counts


def WriteCodonCounts(codon_counts, f):
    """Writes counts of codon identities.
   
    This function takes a dictionary that specifies the number of
    times each codon identity was observed at a position. This
    is used to write a text output file.

    This function operates exactly like `WriteNTCounts`, except
    that it operates on codon rather than nucleotide counts,
    so `codon_counts` should be keyed by codons rather than
    nucleotides.

    So the required keys for `codon_counts[i]` for each codon `i` (1, 2, ...
    numbering) are:

        'WT' : value should be an upper-case string indicating the
        wildtype identity, such as 'ATG'.

        'COUNTS' : integer total number of counts for this codon.

        All codon codes in `AllCodons()`. For these codes, upper-case
        letters (A, T, C, G) specify identities called by both reads.
        Lower-case letter (a, t, c, g) specify identities called by
        just one read, with the other read ambiguous or excluded at
        the relevant position. N specifies identities where the 
        two reads disagree, both are ambiguous, or both are excluded.
        If any of the required keys are missing, raises a `ValueError`.

    So upper-case letters indicate positions called by both
    reads, lower-case letters indicate positions called by
    one read and the other ambiguous, and N indicates disagreement
    between reads or both ambiguous.
    
    The codons are listed in the order returned by `AllCodons()`.
    """
    if not codon_counts:
        raise ValueError("empty codon_counts")
    maxcodon = max(codon_counts.keys())
    allcodons = AllCodons()
    keys = ['WT', 'COUNTS'] + allcodons
    f.write('# codon counts\n')
    f.write('# upper case indicates nucleotide called by both reads\n')
    f.write('# lower case indicates nucleotide called by one read with other ambiguous\n')
    f.write('# N indicates disagreement between reads, or both ambiguous\n')
    f.write('# POSITION')
    for key in keys:
        f.write(' %s' % key)
    for i in range(1, maxcodon + 1):
        f.write("\n%d" % i)
        if i not in codon_counts:
            raise ValueError("Didn't find key for %d in codon_counts" % i)
        if 'WT' not in codon_counts[i]:
            raise ValueError("Didn't find key for WT in codon_counts[%d]" % i)
        if codon_counts[i]['WT'] not in allcodons:
            raise ValueError("Invalid WT codon of %s" % codon_counts[i]['WT'])
        if codon_counts[i]['WT'] != codon_counts[i]['WT'].upper():
            raise ValueError("WT codon not upper case %s" % codon_counts[i]['WT'])
        f.write(" %s" % codon_counts[i]['WT'])
        for key in keys[1 : ]:
            f.write(" %d" % codon_counts[i][key])


def ReadCodonCounts(f):
    """Reads counts of codon identities written by `WriteCodonCounts`.
  
    This function reads the files written by `WriteCodonCounts` to
    reconstruct the calling dictionary specifying the number
    of time each codon identity is observed at a position.

    `f` : a readable file-like object of the type written by 
    `WriteNTCounts`, open to the position where the codon
    counts begin.

    The returned variable is a dictionary that should be similar
    to the `codon_counts` dictionary used to call `WriteCodonCounts`.
    """
    if not hasattr(f, 'read'):
        raise ValueError("f is not readable file-like object")
    codon_counts = {}
    allcodons = AllCodons()
    keys = ['WT', 'COUNTS'] + allcodons
    for line in f:
        if line[0] == '#' or line.isspace() or not line:
            continue # comment line
        try:
            entries = line.split()
            if len(entries) != len(keys) + 1:
                raise ValueError
        except ValueError:
            raise ValueError("Couldn't parse line:\n%s" % line)
        if entries[1] not in allcodons:
            raise ValueError("Invalid WT of %s" % entries[1])
        position = int(entries[0])
        if position in codon_counts:
            raise ValueError("Duplicate positions of %d" % position)
        id = {keys[0]:entries[1]}
        for (key, value) in zip(keys[1 : ], entries[2 : ]):
            id[key] = int(value)
        codon_counts[position] = id
    if not codon_counts:
        raise IOError("Failed to read any content.")
    maxi = max(codon_counts.keys())
    if maxi < 1:
        raise ValueError("Maximum position < 1")
    for i in range(1, maxi + 1):
        if i not in codon_counts:
            raise ValueError("No information for position %d" % i)
    return codon_counts


def WriteAAsFromCodonCounts(infile, outfile):
    """Reads in codon counts, writes out amino acid counts.
    
    `infile` : name of an existing file of the type created by 
    `WriteCodonCounts`.

    `outfile` : name of the resulting output file, which is created.

    Reads in codon counts from `infile`. For all definitively called
    codons (all upper case nucleotides) writes the number of
    amino acids of that identity. Also writes the total number
    of definitively called codons at that position, and the wildtype
    amino acid. The format is self-explanatory.
    """
    if not os.path.isfile(infile):
        raise ValueError("Cannot find codoncounts infile %s" % infile)
    codon_counts = ReadCodonCounts(open(infile))
    mapmuts.sequtils.ClassifyCodonCounts(codon_counts)
    f = open(outfile, 'w')
    f.write('# amino-acid counts among definitively called codons\n')
    f.write("# POSITION TOTAL_COUNTS WT")
    aas = mapmuts.sequtils.AminoAcids()
    for aa in aas:
        f.write(" N_%s" % aa)
    f.write(" N_STOP\n")
    maxcodon = max([i for i in codon_counts.keys() if isinstance(i, int)])
    for i in range(1, maxcodon + 1):
        wtcodon = codon_counts[i]['WT']
        wtaa = mapmuts.sequtils.Translate([('head', wtcodon)])[0][1]
        if not wtaa:
            wtaa = 'STOP'
        elif wtaa not in aas:
            raise ValueError("Unrecognized aa of %s" % wtaa)
        f.write("%d %d %s" % (i, codon_counts[i]['PAIRED_COUNTS'], wtaa))
        for aa in aas:
            f.write(" %d" % codon_counts[i]['N_%s' % aa])
        f.write(" %d\n" % codon_counts[i]['N_STOP'])
    f.close()


def ParseInfile(f):
    """Reads key / value pairs from an input file.

    `f` : should be a readable file-like object.

    Starting at the current position in `f`, reads all remaining lines
    until the end of the file-like object. Any blank line or line
    beginning with a # character is disregarded (# indicates a comment
    line). All other lines should contain exactly two entries, the first
    being the key and the second being the value. The key is construed
    to be the first word, and cannot have any spaces. The value is all
    of the text that follows the key up to the newline. 

    The key / value pairs are returned in a
    dictionary, as illustrated in the example below. If there is a
    duplicate key, and exception is raised.

    Example of successfully reading two key/value pairs:

    >>> f = cStringIO.StringIO()
    >>> f.write('# comment line followed by blank line\\n\\n')
    >>> f.write('key1 first_value\\n')
    >>> f.write('# now another key with two values\\nkey2 2 value2')
    >>> f.seek(0)
    >>> ParseInfile(f) == {'key1':'first_value', 'key2':'2 value2'}
    True

    Example of duplicate key name:

    >>> f = cStringIO.StringIO()
    >>> f.write('# comment line followed by blank line\\n\\n')
    >>> f.write('key1 first_value\\n')
    >>> f.write('# now another key\\nkey2 2')
    >>> f.write('\\nkey1 1')
    >>> f.seek(0)
    >>> ParseInfile(f) == {'key1':'first_value', 'key2':'2'}
    Traceback (most recent call last):
        ...
    ValueError: duplicate key: key1

    """
    d = {}
    for line in f:
        if not (not line or line.isspace() or line[0] == '#'):
            entries = line.split(None, 1)
            if len(entries) != 2:
                raise ValueError("invalid value/key, not two entries:"\
                        + ' %s' % line.strip())
            key = entries[0].strip()
            value = entries[1].strip()
            if key in d:
                raise ValueError("duplicate key: %s" % key)
            d[key] = value
    return d


def ParseFileList(d, key):
    """Gets list of files from input dictionary.

    `d` : a dictionary of key/value string pairs as returned by
    `ParseInfile`.

    `key` : a string in this dictionary that specifies a value that
    should give a list of one or more file names.

    Returns a list of the file names specified by `key`. If one or more of
    these files does not exist, raises and `IOError`. If `key` does not
    exist in `d`, raises an `KeyError`.
    """
    if key not in d:
        raise KeyError("Did not find key of %s" % key)
    files = []
    for f in d[key].split():
        f = f.strip()
        if not os.path.isfile(f):
            raise IOError("Cannot find specified file %s" % f)
        files.append(f)
    return files


def ParseBoolValue(d, key):
    """Gets Boolean argument from input dictionary.

    `d` : a dictionary of key/value string pairs as returned by
    `ParseInfile`.
    
    `key` : a string in this dictionary that specifies a value that
     should be 'True' or 'False'.

    Returns the Boolean truth value specified by `key`. Raises a `ValueError`
    if `key` does not exist in `d`, and a `ValueError` if it specifies
    something other than `True` or `False`.

    Here are some examples:

    >>> d = {'gzipped':'True', 'applyfilter':'False', 'a1':'a1.txt'}
    >>> ParseBoolValue(d, 'gzipped')
    True

    >>> ParseBoolValue(d, 'applyfilter')
    False

    >>> ParseBoolValue(d, 'a1')
    Traceback (most recent call last):
       ...
    ValueError: value a1.txt for a1 is not True/False

    >>> ParseBoolValue(d, 'otherkey')
    Traceback (most recent call last):
       ...
    ValueError: did not find a key for otherkey

    """
    if key not in d:
        raise ValueError("did not find a key for %s" % key)
    elif d[key] == 'True':
        return True
    elif d[key] == 'False':
        return False
    else:
        raise ValueError("value %s for %s is not True/False" % (d[key],\
                key))


def ParseIntValue(d, key):
    """Gets integer argument from input dictionary.

    `d` : a dictionary of key/value string pairs as returned by
    `ParseInfile`.

    `key` : a string in this dictionary that specifies a value that
    should be an integer.

    Returns the integer specified by `key`. Raises a `ValueError`
    if `key` does not exist in `d`, and a `ValueError` if it specifies
    something other than an integer.

    Here are some examples:

    >>> d = {'maxn':'2', 'minq':'20.5'}
    >>> ParseIntValue(d, 'maxn')
    2

    >>> ParseIntValue(d, 'minq')
    Traceback (most recent call last):
       ...
    ValueError: value 20.5 for minq is not integer

    >>> ParseIntValue(d, 'otherkey')
    Traceback (most recent call last):
       ...
    ValueError: did not find a key for otherkey

    """
    if key not in d:
        raise ValueError("did not find a key for %s" % key)
    try:
        return int(d[key])
    except ValueError:
        raise ValueError('value %s for %s is not integer' % \
                (d[key], key))


def ParseFloatValue(d, key):
    """Gets float argument from input dictionary.

    `d` : a dictionary of key/value string pairs as returned by
    `ParseInfile`.

    `key` : a string in this dictionary that specifies a value that
    should be convertible to a float.

    Returns the float specified by `key`. Raises a `ValueError`
    if `key` does not exist in `d`, and a `ValueError` if it specifies
    something other than an integer.

    Here are some examples:

    >>> d = {'maxn':'2', 'minq':'20.5', 'string':'hello'}
    >>> print "%.1f" % ParseFloatValue(d, 'maxn')
    2.0

    >>> print "%.1f" % ParseFloatValue(d, 'minq')
    20.5

    >>> ParseFloatValue(d, 'string')
    Traceback (most recent call last):
       ...
    ValueError: value hello for string is not float

    >>> ParseFloatValue(d, 'otherkey')
    Traceback (most recent call last):
       ...
    ValueError: did not find a key for otherkey

    """
    if key not in d:
        raise ValueError('did not find a key for %s' % key)
    try:
        return float(d[key])
    except ValueError:
        raise ValueError('value %s for %s is not float' % \
                (d[key], key))


def ParseSeqValue(d, key):
    """Reads sequence from FASTA file specified by input dictionary.

    `d` : a dictionary of key/value string pairs as returned by
    `ParseInfile`.

    `key` : a string in this dictionary that specifies a value that is a
    filename of a readable file that contains exactly one sequence
    in FASTA format.

    Returns a string corresponding to the sequence specified in the
    FASTA file. Raises a `ValueError` if the filename does not exist, or
    does not specify exactly one sequence.
    """
    if key not in d:
        raise ValueError("did not find a key for %s" % key)
    filename = d[key]
    if not os.path.isfile(filename):
        raise ValueError("Cannot find file %s specified by %s" % \
                (filename, key))
    seqs = mapmuts.sequtils.ReadFASTA(filename)
    if len(seqs) != 1:
        raise ValueError("file %s does not contain exactly one"\
                % filename + ' sequence')
    return seqs[0][1]


def ParseStringValue(d, key):
    """Reads string argument specified by input dictionary.

    `d` : a dictionary of the key/value string pairs as returned by
    `ParseInfile`.

    `key` : a string in this dictionary that specifies a value that
    is a string.

    Returns string corresponding to `key`. Raises a `ValueError` if 
    `key` does not exist as a key in `d`.

    Here are some examples:

    >>> d = {'outfileprefix':'test_example', 'samplename':'test example'}
    >>> ParseStringValue(d, 'outfileprefix')
    'test_example'

    >>> ParseStringValue(d, 'samplename')
    'test example'

    >>> ParseStringValue(d, 'otherkey')
    Traceback (most recent call last):
       ...
    ValueError: did not find a key for otherkey

    """
    if key not in d:
        raise ValueError("did not find a key for %s" % key)
    return d[key].strip()


def AllNTs():
    """List of all nucleotide codes for `ReadNTCounts` / `WriteNTCounts`.

    Returns a list of all nucleotide codes in the scheme used by
    `ReadNTCounts` / `WriteNTCounts`. The order of entries in this list
    is important, and should be preserved to maintain compatibility
    for reading/writing files.

    Example:

    >>> AllNTs()
    ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'N']

    """
    return ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'N']


def AllCodons():
    """List of all codon codes for `ReadCodonCounts` / `WriteCodonCounts`.

    Returns a list of all codon codes in the scheme used by
    `ReadCodonCounts` / `WriteCodonCounts`. The order of entries in this list
    is important, and should be preserved to maintain compatibility
    for reading/writing files.

    Example:

    >>> AllCodons()[ : 9]
    ['AAA', 'AAT', 'AAC', 'AAG', 'AAa', 'AAt', 'AAc', 'AAg', 'AAN']

    >>> len(AllCodons()) == 9**3
    True

    """
    all = []
    for nt1 in AllNTs():
        for nt2 in AllNTs():
            for nt3 in AllNTs():
                all.append("%s%s%s" % (nt1, nt2, nt3))
    return all


if __name__ == '__main__':
    import doctest
    doctest.testmod()
