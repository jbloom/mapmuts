"""Module for performing main script functions for `mapmuts` package.

This module contains higher-level functions that are designed to be
called fairly directly by scripts to perform complete operations in the
mapping of sequencing reads to genes.

List of functions
-------------------

`MakeAlignments` : Takes input FASTQ files for paired reads
and creates alignments of these reads to a target gene. Also creates
output files summarizing the alignments and reads.

`MakeAlignmentsPlots` : creates plots summarizing output of
MakeAlignments.

`ParseNTCodonCounts` : parses frequency of different nucleotide and codon
variants from the ``*_alignments.txt.gz`` file built by `MakeAlignments`.

`ParseNTCodonCountsPlots` : creates plots summarizing the output of
`ParseNTCodonCounts`.

Documentation for functions
-----------------------------

Documentation for the individual functions is provided in their 
definitions below.

"""


import os
import sys
import gzip
import time
import mapmuts.sequtils
import mapmuts.io
import mapmuts.align
import mapmuts.plot
import mapmuts.latex


def MakeAlignments(r1files, r2files, gzipped, applyfilter, minq, fullgene,\
        generange, a1, a2, maxn, minoverlap, maxrm, maxa1m, maxa2m,\
        maxgenem, outfileprefix, upcase='test', logevery=int(1e5),\
        log=None, write_unaligned=False):
    """Aligns overlapping paired-end reads to a gene.

    This function aligns overlapping paired-end sequencing reads to a
    gene sequence. It has been tested on the alignment of 50 nt
    paired-end Illumina reads to a roughly 1.5 kb gene, although it
    should work more generally. The result of running this function is a
    variety of output files containing the actual alignments and summary
    information about the alignment process.

    Some important notes:

    * The alignments are not actually true alignments, but rather direct
      searches for substrings allowing some mismatches. This works
      efficiently for short Illumina reads (which tend not to have gap
      mutations) being aligned to relatively short genes (a few kb). It
      will become increasingly less computationally efficient for longer
      reads or gene, and is not a good option if the main mode of
      sequencing errors is insertions/deletions.

    * N / n nucleotides do not count as mismatches, regardless of what they
      are aligned with. Instead, there is simply a limit on how many of
      these can occur in a read as specified by `maxn`.

    * The alignments are case-sensitive. In general, you will want to
      guarantee that the sequences are all of the same case. See the 
      `upcase` option.

    This function has been tested on FASTQ files generated from the
    Illumina Casava 1.8 pipeline for paired-end reads. It should work
    more generally on FASTQ files, but that has not been tested.

    The general procedure is as follows:

    1) Paired Illumina reads (R1 and R2 reads) are read from a FASTQ
       file. If the applyfilter option is True, any reads that fail the
       Illumina chastity filter are discarded. In addition, any reads
       that fail the Q-score cutoff specified by `minq` are discarded.

    2) The R1 and R2 reads are aligned to each other assuming that reads
       have the terminal adaptor sequences specified by `a1` and `a2`. The
       alignment must meet the standards specified by `maxn`, `minoverlap`,
       `maxrm`, `maxa1m`, and `maxa2m`. If the paired reads cannot be aligned
       to each other, they are discarded.

    3) The overlapped region of the R1 and R2 reads is aligned to the
       full template sequence given by `fullgene`. If the region cannot
       be aligned with <= `maxgenem` mismatches, the reads are discarded.

    4) The portion of the overlap to the gene range of interest
       (specified by *generange*) is written to an output file, assuming
       that some of the overlap is to the gene range of interest.

    5) Other output files are generated with statistics from the
       alignment process.

    CALLING PARAMETERS: 

    `r1files` : list of file name(s) containing the R1 reads.

    `r2files` : list of file name(s) containing the R2 reads. The files
    specified here must specify the exact same number of reads as
    those in `r1files`, and they must be in the same sequence (i.e.
    the first R1 read found going through `r1files` in order must
    match the first R2 read found going through `r2files` in order).

    `gzipped` : a Boolean switch specifying whether the FASTQ files
    specified by `r1files` and `r2files` are gzipped. If set to `True`,
    they are gzipped -- they will be read without being unzipped.

    `applyfilter` : a Boolean switch specifying whether we remove read
    pairs in which one or more of the reads failed the Illumina
    chastity filter. If `True`, all reads that failed the Illumina
    chastity filter are discarded. You probably want this option to
    be `True` unless you have a good reason to do otherwise.

    `minq` : specifies the minimum acceptable average Q-score (with the
    average taken over all positions in the read). Any read pair in
    which either read has an average Q-score < `minq` is discarded.

    `fullgene` : a string giving the full sequence of the gene that is
    being sequenced (i.e. the template amplicon). Note that this
    full sequence will typically be longer than the actual coding
    gene range of interest, for example due to terminal PCR
    sequences. The substring of `fullgene` of interest is specified by
    `generange`, and is denoted as `gene` in the remaining
    documentation.

    `generange` : a 2-tuple `(gstart, gend)` that specifies the portion of
    `fullgene` that contains the gene of interest. Specifically, ``gene
    = fullgene[generange[0] : generange[1]]``.

    `a1` : a string giving the adaptor sequence that will be found at the
    3' end of any R1 reads that read past the insert and into the
    adaptor.

    `a2` : a string giving the adaptor sequence that will be found at the
    3' end of any R2 reads that read past the insert and into the
    adaptor.

    `maxn` : the maximum number of N or n nucleotides that are allowed
    in either read of a read pair. If one of the reads has > `maxn`
    N / n nucleotides, the pair is discarded.

    `minoverlap` : the minimum length of the overlap between R1 and R2,
    representing their joint coverage of the target sequence. The
    overlap length must be >= this number.

    `maxrm` : the maximum number of mismatches allowed in the overlap of
    R1 and R2. The number of mismatches must be <= this number.

    `maxa1m` : the maximum number of mismatches allowed in the overlap of
    R1 with adaptor `a1`. The number of mismatches must be <= this
    number.

    `maxa2m` : the maximum number of mismatches allowed in the overlap of
    R2 with adaptor `a2`. The number of mismatches must be <= this
    number.

    `maxgenem` : the maximum number of mismatches that is allowed for
    either read individually when it is aligned with `fullgene` after
    removing the adaptor sequences (`a1` or `a2`).
    
    `outfileprefix` : all output files generated by this method have this
    prefix appended to the suffix described below.

    `upcase` : specifies how we handle possible upper / lower case
    differences in nucleotide codes. The alignments handle upper and
    lower case letters as different. If `upcase` has its default value
    of 'test', then `fullgene`, `a1`, and `a2` are all converted to upper
    case. For the reads, the very first read pair is tested to make
    sure that it is upper case -- if it is not, an exception is
    raised. You can also set `upcase` to `False` -- in this case, no
    case conversion is done for any of the sequences. Finally, you
    can set `upcase` to `True` -- in this case, all sequences are
    converted to upper case. This will guarantee no case mismatches,
    but will also be slower due to the extra time of converting all
    sequences to upper case.

    `logevery` : an integer specifying that we print an update to the log
    output file each time we have finished this many read pairs.

    `log` : an optional argument specifying that we write the log of
    this function's progress to a specified file-like object rather
    than to the file ``"%s_log.txt" % outfileprefix``. By default,
    `log` is `None`, meaning that the output is logged to the aforementioned
    newly created file. If `log` is instead set to a writable file-like
    object, we write the output to this file-like object.

    `write_unaligned` : a Boolean switch specifying that we write a
    file containing reads that cannot be aligned. `False` by default,
    meaning that no such file is written. If set to `True`, then a 
    file ``"%s_unaligned.fasta.gz" % outfileprefix`` is written
    containing unaligned reads. If `applyfilter` is `True`, this
    file does NOT contain reads that fail the filter.

    OUTPUT FILES:

    The output is a series of files. All of these files begin with
    `outfileprefix`, conforming to ``"%s_suffix" % outfileprefix``. Any
    existing files with the same names are overwritten. We have
    the following files with the following suffixes:

    "alignments.txt.gz" : A gzipped text file containing the read pair
    to gene alignments for all read pairs that pass the various
    filters and cutoffs. Each alignment consists of five lines.
    These lines give: 
    
    1) the read name,

    2) the alignment region of the gene with indices, 

    3) the alignment region of R1 with indices, 

    4) the alignment region of R2 with indices

    5) a blank line. 
    
    The indices indicate the number of the first and last nucleotide shown
    in 1, 2, ... numbering. For the reads, a '.' is used to indicated
    identities with the gene at that position. The alignment region
    is only shown for regions of gene covered by both reads. For
    example, here is a possible entry::

        @DH1DQQN1:241:C1433ACXX:2:1101:1623:2365
        14 ATAGATACTAGAT 26
         1 .....C.....T. 13
        13 .....C.......  1

    "log.txt" : A text file logging progress of the function. This can
    be tracked in real time to see how we are proceeding. Note
    that if the `log` option is set a file-like object rather than
    its default value of `None`, this file is not created and written.

    "alignmentstatistics.txt" : A text file giving a summary of the
    alignment, in the format written by `io.WriteAlignmentStatistics`.

    "insertlengths.txt" : A text file giving the distribution of insert
    lengths for all R1-R2 read pairs that align with each other
    according to the given specifications. The shortest possible
    value of the insert is `minoverlap`, while the longest possible
    value is the sum of the two read lengths minus `minoverlap`. In
    the format written by `io.WriteInsertLengths`.

    "R1mismatches.txt" : A text file giving the fraction and number of
    mismatches between the R1 read and the `fullgene` template for all
    positions in the read that align with `fullgene`, for reads that
    pass the pairwise alignment with each other. In the format
    written by `io.WriteReadMismatches`.

    "R2mismatches.txt" : Like "R1mismatches.txt" but for the R2 reads.

    "unaligned.fasta.gz" : Only created if `write_unaligned` is `True`,
    contains unaligned reads that pass `applyfilter` option.

    C-EXTENSIONS:

    This function uses all available C-extensions to speed up 
    the code as much as possible.
    """
    if log:
        closelog = False
        if not isinstance(log, file):
            raise ValueError("log is not file-like object")
        try:
            log.write('')
        except IOError:
            raise IOError("log file-like object is not writable")
    else:
        closelog = True
        log = open("%s_log.txt" % outfileprefix, 'w')
    try:
        start = time.clock()
        log.write("Beginning execution of " +\
                "mapmuts.main.MakeAlignments.\n")
        log.write("Current time is %s.\n" % time.ctime())
        log.write("Current directory is %s.\n" % os.getcwd())
        for (r, files) in [('R1', r1files), ('R2', r2files)]:
            log.write("\nThe %s reads will come from the following" % r\
                    + " FASTQ files:\n")
            for f in files:
                if not os.path.isfile(f):
                    raise IOError("Cannot find file %s" % f)
                log.write("%s\n" % f)
        if applyfilter:
            log.write('\napplyfilter = True: all reads flagged Y '\
                    + ' by the Illumina filter will be removed.\n')
        else:
            log.write('\napplyfilter = False: even reads flagged Y '\
                    + ' by the Illumina filter will be retained.\n')
        log.write('\nminq = %.2f: any read pair where the average Q-'\
                % minq + 'score of either read is < this will be removed.\n')
        if upcase == 'test':
            log.write("\nupcase = 'test': converting fullgene, a1, and"\
                    + " a2 to upper case.\n")
            fullgene = fullgene.upper()
            a1 = a1.upper()
            a2 = a2.upper()
            log.write("Testing case of first R1 and R2 reads...")
            for (n, r1, r2, q1, q2) in mapmuts.io.IteratePairedFASTQ(\
                    r1files, r2files, gzipped, applyfilter=False, usegzip=True):
                break
            if (r1 == r1.upper()) and (r2 == r2.upper()):
                log.write(" test passed, first reads upper case.\n")
            else:
                raise ValueError("The first R1 and R2 reads are lower"\
                        + " case, so failing upcase='test' test.")
            upcase = False # no further upcase conversions will be done
        elif upcase == True:
            log.write("\nupcase = True: all sequences will be converted "\
                    + "to upper case.\n")
            fullgene = fullgene.upper()
            a1 = a1.upper()
            a2 = a2.upper()
        elif upcase == False:
            log.write("\nupcase = False: No conversion to upper case "\
                    + "will be done for any sequences.\n")
        else:
            raise ValueError("Invalid value of upcase: %r" % upcase)
        log.write("\nThe value of fullgene is:\n%s\n" % fullgene)
        if not (isinstance(generange, tuple) and 0 <= generange[0] <\
                generange[1] <= len(fullgene)):
            raise ValueError("Invalid value of generange: %s" % \
                    str(generange))
        log.write("\nThe value of generange is:\n%s\n" % str(generange))
        log.write("\nThis means that the value of gene (the region of"\
                + " fullgene specified by generange) is:\n%s\n" % \
                fullgene[generange[0] : generange[1]])
        log.write("\nThe value of a1 (the adaptor at the 3' end of R1"\
                + " reads) is:\n%s\n" % a1)
        log.write("\nThe value of a2 (the adaptor at the 3' end of R2"\
                + " reads) is:\n%s\n" % a2)
        log.write('\nThe value of maxn (the maximum number of N / n '\
                + 'nucleotides allowed in a read) is %d\n' % maxn)
        log.write("\nThe value of minoverlap (minimum acceptable "+\
                "overlap between R1 and R2) is %d.\n" % minoverlap)
        log.write("\nThe value of maxrm (maximum allowed mismatches "+\
                "between R1 and R2 in overlap) is %d.\n" % maxrm)
        log.write("\nThe value of maxa1m (maximum allowed mismatches "+\
                "between R1 and its adaptor a1) is %d.\n" % maxa1m)
        log.write("\nThe value of maxa2m (maximum allowed mismatches "+\
                "between R2 and its adaptor a2) is %d.\n" % maxa2m)
        log.write("\nThe value of maxgenem (maximum allowed mismatches"\
                + " of either read with fullgene after removing read "\
                + "adaptors) is %d.\n" % maxgenem)
        if write_unaligned:
            unalignedfilename = '%s_unaligned.fasta.gz' % outfileprefix
            unalignedfile = gzip.open(unalignedfilename, 'w')
            log.write('\nUnaligned reads that passed applyfilter will'\
                    + ' be written to %s\n' % unalignedfilename)
        aligner = mapmuts.align.ReadPairToGeneAligner(fullgene,\
                generange, a1, mapmuts.sequtils.ReverseComplement(a2), maxn,\
                minoverlap, maxrm, maxa1m, maxa2m, maxgenem, upcase,\
                use_cext=True)
        nread = nfiltered = nlowq = 0
        alignfilename = "%s_alignments.txt.gz" % outfileprefix
        alignfile = gzip.open(alignfilename, 'w')
        log.write("\nSuccessfully aligned reads will be written to"\
                + " %s\n" % alignfilename)
        log.write("\nBeginning reading and aligning reads...\n")
        log.flush()
        for tup in mapmuts.io.IteratePairedFASTQ(r1files, r2files,\
                gzipped, applyfilter):
            nread += 1
            if tup:
                (n, r1, r2, q1, q2) = tup 
                if q1 < minq or q2 < minq:
                    nlowq += 1
                    if write_unaligned:
                        unalignedfile.write(">%s R1 avgQ=%.1f REASON_UNALIGNED:lowQ\n%s\n" \
                                % (n, q1, r1) + '>%s R2 avgq=%1.f REASON_UNALIGNED:lowQ\n%s\n'\
                                % (n, q2, r2))
                else:
                    a = aligner.Align(r1,\
                            mapmuts.sequtils.ReverseComplement(r2))
                    if isinstance(a, str):
                        alignfile.write("%s\n%s\n\n" % (n, a))
                    elif isinstance(a, tuple):
                        reason = a[1]
                        if write_unaligned:
                            unalignedfile.write(">%s R1 avgQ=%.1f REASON_UNALIGNED:%s\n%s\n" \
                                    % (n, q1, reason, r1) + '>%s R2 avgq=%1.f REASON_UNALIGNED:%s\n%s\n'\
                                    % (n, q2, reason, r2))
                    else:
                        raise ValueError("Invalid return for aligner.Align")
            else:
                nfiltered += 1
            if not (nread % logevery):
                t = time.clock() - start
                log.write("Completed %d reads in %.3f seconds.\n" %\
                        (nread, t))
                log.flush()
        alignfile.close()
        log.write("\nNow writing statistics to output files.\n")
        s = aligner.Statistics()
        s['nfiltered'] = nfiltered
        s['nread'] = nread
        s['nlowq'] = nlowq
        summaryfilename = '%s_alignmentstatistics.txt' % outfileprefix
        log.write('Writing summary statistics to %s.\n' % summaryfilename)
        summaryfile = open(summaryfilename, 'w')
        mapmuts.io.WriteAlignmentStatistics(s, summaryfile)
        summaryfile.close()
        insertlenfilename = '%s_insertlengths.txt' % outfileprefix
        log.write('Writing insert length distribution to %s.\n' %\
                insertlenfilename)
        insertlenfile = open(insertlenfilename, 'w')
        mapmuts.io.WriteInsertLengths(s['insertlength'], insertlenfile)
        insertlenfile.close()
        for r in [1, 2]:
            if s['naligned'] == 0:
                log.write("Failed to align any reads, so not writing mismatch distribution.\n")
                continue
            rmfilename = '%s_R%dmismatches.txt' % (outfileprefix, r)
            log.write('Writing R%d mismatch distribution to %s.\n' %\
                    (r, rmfilename))
            rmfile = open(rmfilename, 'w')
            mapmuts.io.WriteReadMismatches(s['r%dlengths' % r],\
                    s['r%dmismatches' % r], rmfile)
            rmfile.close()
        log.write("\nCompleted execution of mapmuts.main.MakeAlignments at %s."\
                % time.ctime())
        if closelog:
            log.close()
        if write_unaligned:
            unalignedfile.close()
    except:
        for x in sys.exc_info():
            log.write("\n%s" % str(x))
        if closelog:
            log.write("\nPrematurely closing log due to execution error!")
            log.close()
        if write_unaligned:
            if 'unalignedfile' in vars():
                unalignedfile.close()
        raise


def MakeAlignmentsPlots(outfileprefix, nlowq, maxn, minoverlap,\
        maxrm, maxa1m, maxa2m, maxgenem, latexsummary, summarytitle):
    """Makes plots summarizing the output of `MakeAlignments`.

    This function creates plots summarizing the results output by the
    `mapmuts.main.MakeAlignments` function. These plots are created using
    ``pylab`` / ``matplotlib``, so this function will raise an exception if
    these packages are not available (if ``mapmuts.plot.PylabAvailable() ==
    False``). If ``pdflatex`` is available (if 
    *mapmuts.latex.PdflatexAvailable()*) this function can also optionally
    create a summary PDF file containing all of these plots.
    All output files from `MakeAlignments` are expected to
    exist with the names created using the value `outfileprefix` used
    to call the `MakeAlignments` function.

    CALLING VARIABLES:

    `outfileprefix` : prefix of output files, value used to call
    `MakeAlignments`.

    `nlowq`, `maxn`, `minoverlap`, `maxrm`, `maxa1m`, `maxa2m`, 
    `maxgenem` : the values of these variables
    used when calling `MakeAlignments`. These values are
    used to output information about the alignment parameters.

    `latexsummary` : Boolean switch specifying that we create the PDF
    summary file using `mapmuts.latex.MakeAlignmentsSummary`. We do
    this if `latexsummary` is `True`.

    `summarytitle` : string giving the title printed at the top of the
    latex summary PDF if it is created.

    OUTPUT PLOTS:

    Plots are created summarizing each of the four output files ``*.txt``
    files created by `MakeAlignments`. These files have the same names as
    the output file, but with the suffix ``.pdf``. So these files are:

    ``"%s_alignmentstatistics.pdf" % outfile_prefix``

    ``"%s_insertlengths.pdf" % outfile_prefix``

    ``"%s_R1mismatches.pdf" % outfile_prefix``

    ``"%s_R2mismatches.pdf" % outfile_prefix``

    The `latexsummary` file (if created) has the name:

    ``"%s_makealignments_summary.pdf" % outfileprefix``
    """
    if not mapmuts.plot.PylabAvailable():
        raise ImportError("pylab / matplotlib are not available")
    mapmuts.plot.PlotAlignmentStatistics(["%s_alignmentstatistics.txt"\
            % outfileprefix], [''], "%s_alignmentstatistics.pdf" %\
            outfileprefix)
    mapmuts.plot.PlotInsertLengths("%s_insertlengths.txt" %\
            outfileprefix, "%s_insertlengths.pdf" % outfileprefix,\
            normed=True)
    for r in ['R1', 'R2']:
        mapmuts.plot.PlotReadMismatches("%s_%smismatches.txt" %\
                (outfileprefix, r), "%s_%smismatches.pdf" %\
                (outfileprefix, r),\
                logy=True, title='%s mismatches' % r)
    if latexsummary:
        mapmuts.latex.MakeAlignmentsSummary(outfileprefix,\
                nlowq, maxn, minoverlap, maxrm, maxa1m, maxa2m, maxgenem,\
                summarytitle)


def ParseNTCodonCounts(alignmentfile, outfileprefix, gene, r1exclude,\
        r2exclude, upcase='test', logevery=int(1e5), log=sys.stdout):
    """Parses nucleotide and codon counts from alignment from `MakeAlignments`.

    This function parses information from the ``*_alignments.txt.gz`` file
    created by `MakeAlignments`. The counts of nucleotide and codon variants
    at each position is parsed, and written to two output files.

    The overlapping paired-end read alignments are used to call nucleotide
    identities. Nucleotides are called as follows:

    * If both reads agree on the identity and neither relevant read position
      is in the exclusion lists (`r1exclude` or `r2exclude`), then the
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

    CALLING VARIABLES:

    `alignmentfile` : the name of a ``*_alignments.txt.gz`` file built
    by `MakeAlignments`. It is crucial that the alignments be written
    in all upper case nucleotides for this method to work correctly.
    See the `upcase` option.

    `outfileprefix` : the prefix used at the beginning of the output files. The
    two output files created by this function have the following two names::
    
        "%s_ntcounts.txt" % outfileprefix
        "%s_codoncounts.txt" % outfileprefix

    `gene` : a string giving the gene sequence for the gene used to build
    alignmentfile. If `fullgene` and `generange` are the calling variables
    to `MakeAlignments`, ``gene = fullgene[generange[0] : generange[1]]``.
    The gene is converted to upper case (see the `upcase` option).
    In addition, it must be a protein-coding
    sequence with length that is a multiple of three and no stop codons
    except for possibly the last codon. This is checked, and there will
    be an exception raised if gene is not a valid protein-coding gene.

    `r1exclude` : a list, dictionary, or set that specifies the integer
    positions of any sites in the R1 reads that are excluded from
    being considered, in 1, 2, ... numbering. `r1exclude` must support
    the ``in`` operator, such that we exclude any position i for which 
    ``i in r1exclude == True``. You would want to specify read positions
    that have high error rates, and so for which you do not want
    to include when identifying variants.

    `r2exclude` : like `r1exclude`, but for the R2 reads.

    `upcase` : specifies how we handle possible upper / lower case
    differences in nucleotide codes. Note that because of the coding
    used in calling the identities, this method will not work if some
    of the nucleotides are lower case. The alignments handle upper and
    lower case letters as different. If upcase has its default value
    of 'test', then gene is converted to upper case. For the alignments,
    the very first alignment is tested to make sure that it is upper case
    -- if it is not, an exception is raised. You can also set upcase to
    `True`, in which case all sequences are converted to upper case. This
    is the safest option, but is slwo due to all the case conversion.

    `logevery` : if `log` (the next parameter) evaluates `True`, write updated progress
    every time we have processed this many read alignments. Is set to
    ``int(1e5)`` by default.

    `log` : a file-like object to which we write progress of this function.
    Is set to ``sys.stdout`` by default. You can also set it to any other
    writable file-like object, or to `None` if you do not want to log
    any progress.

    OUTPUT FILES:

    ``"%s_ntcounts.txt" % outfileprefix`` : a file giving the nucleotide
    identity counts at each position, in the format written by
    `mapmuts.io.WriteNTCounts`.

    ``"%s_codoncountx.txt" % outfileprefix`` : a file giving the
    codon identity counts at each position, in the format written
    by `mapmuts.io.WriteCodonCounts`.
    """
    r1exclude = set(r1exclude) # convert two sets, necessary for calign
    r2exclude = set(r2exclude)
    ntcountsfile = "%s_ntcounts.txt" % outfileprefix
    codoncountsfile = "%s_codoncounts.txt" % outfileprefix
    if log:
        start = time.clock()
        log.write("Beginning execution of mapmuts.main.ParseNTCodonCounts.\n")
        log.write("Current time is %s.\n" % time.ctime())
        log.write("Current directory is %s.\n" % os.getcwd())
        log.write("Alignments will be read from %s\n" % alignmentfile)
        log.write("The nucleotide counts output will be written to %s"\
                % ntcountsfile)
        log.write("\nThe codon counts output will be written to %s"\
                % codoncountsfile)
        if r1exclude:
            log.write("\nThe following positions will be excluded in R1:\n")
            for i in r1exclude:
                log.write("%d\n" % i)
        else:
            log.write("\nNo positions will be excluded in R1.\n")
        if r2exclude:
            log.write("The following positions will be excluded in R2:\n")
            for i in r2exclude:
                log.write("%d\n" % i)
        else:
            log.write("\nNo positions will be excluded in R2.\n")
    if upcase == 'test':
        if log:
            log.write("\nupcase = 'test': ")
            log.write("Testing case of first R1 and R2 reads...")
        for aligntup in mapmuts.io.IterateAlignmentFile(alignmentfile,\
                gzipped=True, usegzip=True):
            break
        (gseq, r1seq, r2seq) = (aligntup[2], aligntup[5], aligntup[8])
        if (gseq == gseq.upper()) and (r1seq == r1seq.upper()) and\
                (r2seq == r2seq.upper()):
            if log:
                log.write(" test passed, first alignment upper case.\n")
        else:
            raise ValueError("The alignment contains lower"\
                    + " case, so failing upcase='test' test.")
    gene = gene.upper()
    if len(gene) % 3:
        raise ValueError("gene length of %d is not a multiple of 3" % len(gene))
    try:
        prot = mapmuts.sequtils.Translate([('gene', gene)])[0][1]
    except:
        sys.stderr.write("Failed to translate gene. Must be a coding sequence"\
                + ' and stop codon is only allowed at the terminus.')
        raise
    nnts = len(gene)
    ncodons = nnts // 3
    if len(prot) == ncodons:
        pass
    elif len(prot) == ncodons + 1:
        # last codon is stop codon
        if not mapmuts.sequtils.Translate([('lastcodon', gene[len(gene) - 3 : ])])[0][1]:
            # indeed a stop codon, indicate with *
            prot = "%s*" % prot
        else:
            raise ValueError("this should not happen, problem with last codon")
    else:
        raise ValueError("this should not happen, problem with translation?")
    if log:
        log.write("\nAligning to the following gene sequence (length %d):\n%s\n"\
                % (nnts, gene))
        log.write("\nThis gene has the following translation (%d codons):\n%s\n"\
                % (ncodons, prot))
    allnts = mapmuts.io.AllNTs()
    nt_counts = {}
    for i in range(1, nnts + 1):
        d = {'WT':gene[i - 1].upper(), 'COUNTS':0}
        for nt in allnts:
            d[nt] = 0
        nt_counts[i] = d
    allcodons = mapmuts.io.AllCodons()
    codon_counts = {}
    for i in range(1, ncodons + 1):
        d = {'WT':gene[3 * i - 3 : 3 * i].upper(), 'COUNTS':0}
        for codon in allcodons:
            d[codon] = 0
        codon_counts[i] = d
    i = 0
    if not os.path.isfile(alignmentfile):
        raise IOError("Cannot find alignmentfile %s" % alignmentfile)
    log.write("\nNow reading alignments from %s...\n" % alignmentfile)
    log.flush()
    for aligntup in mapmuts.io.IterateAlignmentFile(alignmentfile, gzipped=True):
        (gstart, gend, gseq) = (aligntup[0], aligntup[1], aligntup[2])
        if upcase == True:
            aligntup = (aligntup[0], aligntup[1], aligntup[2].upper(),\
                    aligntup[3], aligntup[4], aligntup[5].upper(),\
                    aligntup[6], aligntup[7], aligntup[8].upper())
        assert 1 <= gstart <= gend <= nnts, "Invalid gene start/end"
        assert gseq == gene[gstart - 1 : gend], "gene mismatch - did you remember to make upper case?"
        (ntstartindex, ntidentities) = mapmuts.align.ParsePairedAlignment(\
                aligntup, r1exclude, r2exclude)
        lena = len(ntidentities)
        for j in range(lena):
            di = nt_counts[ntstartindex + j]
            di['COUNTS'] += 1
            di[ntidentities[j]] += 1
        (icodon, ipos) = divmod(ntstartindex - 1, 3)
        if ipos:
            startcodon = icodon + 2
        else:
            startcodon = icodon + 1
        if ipos == 0:
            endcodon = startcodon + lena // 3 - 1
            codonshift = 0
        elif ipos == 1:
            endcodon = startcodon + (lena - 2) // 3 - 1
            codonshift = 2
        elif ipos == 2:
            endcodon = startcodon + (lena - 1) // 3 - 1
            codonshift = 1
        for j in range(endcodon - startcodon + 1):
            di = codon_counts[startcodon + j]
            k = j * 3 + codonshift
            di['COUNTS'] += 1
            di[ntidentities[k : k + 3]] += 1
        i += 1
        if not (i % logevery):
            t = time.clock() - start
            log.write("Read %d alignments in %.3f seconds...\n" % (i, t))
            log.flush()
    log.write("Finished reading alignments.\n")
    log.write("\nNow writing nucleotide counts to %s\n" % ntcountsfile)
    f = open(ntcountsfile, 'w')
    mapmuts.io.WriteNTCounts(nt_counts, f)
    f.close()
    log.write("\nNow writing codon counts to %s\n" % codoncountsfile)
    f = open(codoncountsfile, 'w')
    mapmuts.io.WriteCodonCounts(codon_counts, f)
    f.close()
    log.write('\nFinished executation of mapmuts.main.ParseNTCodonCounts '\
            + 'at %s.' % time.ctime())


def ParseNTCodonCountsPlots(outfileprefix, r1exclude, r2exclude,\
        latexsummary, summarytitle):
    """Makes plots summarizing the output of `ParseNTCodonCounts`.

    This function creates plots summarizing the results output by the
    `mapmuts.main.ParseNTCodonCounts` function. These plots are created using
    ``pylab`` / ``matplotlib``, so this function will raise an exception if
    these packages are not available (if ``mapmuts.plot.PylabAvailable() ==
    False``). If ``pdflatex`` is available (if 
    *mapmuts.latex.PdflatexAvailable()*) this function can also optionally
    create a summary PDF file containing all of these plots.
    All output files from `ParseNTCodonCounts` are expected to
    exist with the names created using `outfileprefix` as the calling variable
    to that function. 

    CALLING VARIABLES:

    `outfileprefix` : prefix of output files, used to call `MakeAlignments`.

    `r1exclude`, `r2exclude` : the values of the variables used when
    calling `ParseNTCodonCounts`. These values are used
    to output information about the parsing parameters.

    `latexsummary` : Boolean switch specifying that we create the PDF
    summary file using `mapmuts.latex.ParseNTCodonCountsSummary`. We do
    this if `latexsummary` is `True`.

    `summarytitle` : string giving the title printed at the top of the
    latex summary PDF if it is created.

    OUTPUT FILES

    ``"%s_codondepth.pdf" % outfile_prefix`` : summarizes codon depth
    across the gene.

    ``"%s_syn-ns-dist.pdf" % outfile_prefix`` : summarizes distribution
    of synonymous / nonsynonymous mutations.

    ``"%s_nmutspercodon-dist.pdf" % outfile_prefix`` : summarizes
    distribution of number of nucleotide mutations per codon mutation.

    ``"%s_parsesummary.pdf" % outfileprefix`` : LaTex summary file if
    created (if ``latexsummary == True``).
    """
    if not mapmuts.plot.PylabAvailable():
        raise ImportError("pylab / matplotlib are not available")
    codonfile = "%s_codoncounts.txt" % outfileprefix
    if not os.path.isfile(codonfile):
        raise IOError("Cannot find expected codonfile %s" % codonfile)
    mapmuts.plot.PlotCodonDepth(codonfile, "%s_codondepth.pdf" % outfileprefix)
    mapmuts.plot.PlotCodonDist(codonfile, "%s_syn-ns-dist.pdf"\
            % outfileprefix, dist_type='syn_ns_stop')
    mapmuts.plot.PlotCodonDist(codonfile, "%s_nmutspercodon-dist.pdf"\
            % outfileprefix, dist_type='codon_nmuts')
    if latexsummary:
        mapmuts.latex.MakeParseSummary(outfileprefix, r1exclude,\
                r2exclude, summarytitle)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
