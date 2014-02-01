"""Module for performing sequence and read alignment operations.

This module is part of the ``mapmuts`` package.

Written by Jesse Bloom.


List of functions
-------------------
`AlignReads` : Aligns overlapping ungapped paired-end reads with known adaptors.

`AlignReadToGene` : aligns read to a gene, does not handle gaps.

`Needle` : runs the EMBOSS needle program.

`AddDots` : adds dots at identities for two sequences.

`RemoveDots` : inverse of `AddDots`.

`ReadToGeneMismatches` : finds mismatches between aligned read and gene.

`ParsePairedAlignment` : parses alignment returned by `io.IterateAlignmentFile`.


List of classes
-------------------
`ReadPairToGeneAligner` : aligns overlapping paired-end reads with known adaptors to a gene.


List of C-extensions
--------------------------
Some of the functions in this module are implemented in faster format
the C extension module `calign`. Calling these functions cause them to 
in turn call the faster C implementations, assuming that `calign` can
be imported. This is done for the following functions unless they are
called with `use_calign set` to `False`:

* `AlignReads`

* `AlignReadToGene`

* `AddDots`

* `ReadToGeneMismatches`

* `ParsePairedAlignment`


Documentation for functions and classes
-------------------------------------------
Documentation for individual functions and classes is provided
in their definitions below.
"""

import math
import os
import sys
import tempfile
import subprocess
import warnings
import mapmuts.sequtils
try:
    import mapmuts.calign
    _calignimported = True # global variable, is calign available?
except ImportError:
    _calignimported = False
    warnings.warn("Cannot import calign. Will have to use pure python"\
            + ' implementations.', RuntimeWarning)


def AddDots(s1, s2, use_calign=True):
    """Adds dots at identities for two aligned sequences.
   
    Take as input two string `s1` and `s2` that represent aligned
    sequences. Adds a dot to every position of `s2` that is identical to
    the corresponding position of `s1`. Upper and lower case letter are
    treated differently. Then returns the copy of `s2` with the dots added.

    `s1` : a string giving the first sequence aligned to `s2`

    `s2` : a string giving the second sequence aligned to `s1`
    
    `use_calign` : specifies that we actually perform the calculations
    using the implementation of this function in `calign`. This 
    dramatically speeds the computation. The `calign` extension
    is used assuming that it is available and that this
    parameter is set to its default value of `True`. If `use_calign`
    is set to `False`, or if `calign` is unavailable, then the 
    computations instead use the pure Python implementation.

    >>> s1 = 'ATGACATAGACA'
    >>> s2 = 'ATGACTTACACA'
    >>> AddDots(s1, s2)
    '.....T..C...'

    >>> s1 = 'ATGNATAGACA'
    >>> s2 = 'ATGCANAGACA'
    >>> AddDots(s1, s2)
    '...C.N.....'
    
    >>> s1 = 'ATGACATAGACA'
    >>> s2 = 'aTGACTTACACA'
    >>> AddDots(s1, s2)
    'a....T..C...'

    >>> s1 = 'ATGACATA-ACA'
    >>> s2 = '-TGACTTACACA'
    >>> AddDots(s1, s2)
    '-....T..C...'
    """
    if _calignimported and use_calign:
        return mapmuts.calign.AddDots(s1, s2)
    assert len(s1) == len(s2)
    s2dotted = []
    for (x, y) in zip(s1, s2):
        if x == y:
            s2dotted.append('.')
        else:
            s2dotted.append(y)
    return ''.join(s2dotted)


def RemoveDots(s1, s2):
    """Removes dots from identities for two aligned sequences.

    This function is the inverse of `AddDots`. Takes as input
    two strings `s1` and `s2` of the same length that represent
    aligned sequences with a dot at each position of `s2` that
    matches `s1`. Returns a copy of `s2` where these dots have
    been replaced with the actual identity at that position.

    `s1` : a string giving the first sequence.

    `s2` : a string giving the second sequence with dots added,
    as would be done by calling ``s2 = AddDots(s1, s2a)`` where
    `s2a` is the original second sequence with all of the 
    nucleotide identities.

    >>> s1 = 'ATGACATAGACA'
    >>> s2 = 'ATGACTTACACA'
    >>> s2 == RemoveDots(s1, AddDots(s1, s2))
    True
    
    >>> s1 = 'ATGACATAGACA'
    >>> s2 = 'aTGACTTACACA'
    >>> s2 == RemoveDots(s1, AddDots(s1, s2))
    True

    >>> s1 = 'ATGACATA-ACA'
    >>> s2 = '-TGACTTACACA'
    >>> s2 == RemoveDots(s1, AddDots(s1, s2))
    True
    """
    assert len(s1) == len(s2)
    s2undotted = []
    for (x, y) in zip(s1, s2):
        if y == '.':
            s2undotted.append(x)
        else:
            s2undotted.append(y)
    return ''.join(s2undotted)


def ReadToGeneMismatches(g, r, gi, maxmm, indexback=False,\
        n_is_mismatch=False, use_calign=True):
    """Gets a list of read positions that mismatch with a target gene.

    This function is designed to get a list of all positions where a read
    mismatches with a target gene, presuming that you already know the 
    position in the gene to which the read aligns.

    `r` : a string giving the read sequence.

    `g` : a string giving the gene sequence.

    `gi` : the position in g where r aligns, in 0, 1, ... numbering. We 
    expect `r` to align to ``g[gi : gi + len(r)]``.

    `maxmm` : integer giving the maximum number of allowed mismatches.

    `indexback` : a Boolean switch specifying that we count backwards
    from the end of the read when we label the mismatches, 
    rather than forward from the start. 

    `n_is_mismatch` : a Boolean switch specifying that we consider 'N'
    or 'n' nucleotides to be mismatches if they align with something
    else. By default this is `False`, meaning that an 'N' or 'n'
    nucleotide is not considered to be a mismatch regardless of
    what it aligns with. If set to `True`, then an 'N' or 'n' is 
    considered a mismatch unless it aligns with the same character.

    `use_calign` : specifies that we actually perform the calculations
    using the implementation of this function in `calign`. This 
    dramatically speeds the computation. The `calign` extension
    is used assuming that it is available and that this
    parameter is set to its default value of `True`. If `use_calign`
    is set to `False`, or if `calign` is unavailable, then the 
    computations instead use the pure Python implementation.

    The read and the gene are assumed to be in the same orientation, so
    no checking for reverse-complement matching is done.

    If the read extends beyond the gene (``gi + len(r) > len(g)``) then 
    returns `False`.  Otherwise, returns list of each position in the read 
    (in 1, 2, ... numbering) where `r` differs from the aligned position
    in `g`. Note that the matching is case-sensitive, and the treatment of
    'N' or 'n' is determined by `n_is_mismatch`. If the number of mismatches is
    > `maxmm`, returns False. If indexback is `True`, then the last position of
    `r` is considered 1, the second to last position is 2, etc.

    Example with perfect match, then too many mismatches as read misaligned.

    >>> g = 'ATGGGCATGACCTGA'
    >>> r = 'CATGACC'
    >>> gi = 5
    >>> maxmm = 2
    >>> ReadToGeneMismatches(g, r, gi, maxmm)
    []
    >>> ReadToGeneMismatches(g, r, gi, maxmm, indexback=True)
    []
    >>> gi = 4
    >>> ReadToGeneMismatches(g, r, gi, maxmm)
    False
    >>> ReadToGeneMismatches(g, r, gi, maxmm, indexback=True)
    False

    Example with mismatches

    >>> g = 'ATGGGCATGACCTGA'
    >>> r = 'CATCACC'
    >>> gi = 5
    >>> maxmm = 2
    >>> ReadToGeneMismatches(g, r, gi, maxmm)
    [4]
    >>> ReadToGeneMismatches(g, r, gi, maxmm, indexback=True)
    [4]
    >>> r = 'CATCACCA'
    >>> ReadToGeneMismatches(g, r, gi, maxmm)
    [4, 8]
    >>> ReadToGeneMismatches(g, r, gi, maxmm, indexback=True)
    [5, 1]
    >>> r = 'cATCACC'
    >>> ReadToGeneMismatches(g, r, gi, maxmm)
    [1, 4]
    >>> ReadToGeneMismatches(g, r, gi, maxmm, indexback=True)
    [7, 4]

    Example where read extends beyond gene

    >>> g = 'ATGGGCATGACCTGA'
    >>> r = 'ATCACCTGA'
    >>> gi = 6
    >>> maxmm = 2
    >>> ReadToGeneMismatches(g, r, gi, maxmm)
    [3]
    >>> r = 'ATCACCTGAC'
    >>> ReadToGeneMismatches(g, r, gi, maxmm)
    False

    Example usage of `n_is_mismatch`:

    >>> g = 'ATGGGCATGACCTGA'
    >>> r = 'ATCACCTNA'
    >>> gi = 6
    >>> maxmm = 1
    >>> ReadToGeneMismatches(g, r, gi, maxmm)
    [3]
    >>> ReadToGeneMismatches(g, r, gi, maxmm, n_is_mismatch=True)
    False
    >>> maxmm = 2
    >>> ReadToGeneMismatches(g, r, gi, maxmm, n_is_mismatch=True)
    [3, 8]
    >>> g = 'ATGGGCNTGACCTGA'
    >>> maxmm = 3
    >>> ReadToGeneMismatches(g, r, gi, maxmm)
    [3]
    >>> ReadToGeneMismatches(g, r, gi, maxmm, n_is_mismatch=True)
    [1, 3, 8]
    """
    if _calignimported and use_calign:
        return mapmuts.calign.ReadToGeneMismatches(g, r, gi, maxmm, \
                indexback, n_is_mismatch)
    lenr = len(r)
    imax = gi + lenr
    if imax > len(g):
        return False
    mm = []
    for i in range(gi, imax):
        if r[i - gi] != g[i]:
            if n_is_mismatch or (r[i - gi].upper() != 'N' and\
                    g[i].upper() != 'N'):
                if indexback:
                    mm.append(lenr - i + gi)
                else:
                    mm.append(i + 1 - gi)
                if len(mm) > maxmm:
                    return False
    return mm


def ParsePairedAlignment(aligntup, r1exclude, r2exclude, use_calign=True):
    """Parses nucleotide identities from paired-read alignment.

    `aligntup` : a tuple specifying the alignment of a gene
    with the overlapping region of two paired-end reads. In the
    format returned by `io.IterateAlignmentFile`. Specifically, this
    is the tuple
    `(gstart, gend, gseq, r1start, r1end, r1seq, r2start, r2end, r2seq)`
    This method assumes that the strings specified in the alignment
    are upper case, but does not check this (in the interest of time).
    So you want to make sure that such tests are somehow built in
    elsewhere. Acceptable nucleotide codes are 'A', 'T', 'C', 'G',
    and 'N' (for ambiguous).

    `r1exclude` -> a set that specifies the integer
    positions of any sites in the R1 reads that are excluded from
    being considered, in 1, 2, ... numbering. Note that this must be
    a set for the C-extension in `calign` to work.

    `r2exclude` : like `r1exclude`, but for the R2 reads.

    `use_calign` : specifies that we use fast C extension implemented
    in `calign`.

    For each position in the alignment, calls the
    nucleotide. A nucleotide is definitively called if both the R1 and
    R2 reads in `aligntup` have the same identity. In this case, the nucleotide
    is called as the upper case value: for example, if both reads have 'A'
    then the call is 'A'. If one read gives 'N' or is excluded (by
    `r1exclude` or `r2exclude`) while the other gives a value, then the
    nucleotide is called as the lower case value. For example,
    'A' and 'N' are called as 'a'. If both are called to different 
    non-'N' values, or both are 'N', or both are excluded,
    or one is 'N' and the other is excluded, then it is called as 'N'.

    All position numbering is done using 1, 2, ... indexing.

    The returned variable is the 2-tuple `(ntstartindex, ntidentities)`
    In this 2-tuple, `ntstartindex` is the index of the first nucleotide position
    in the gene that can be called.
    `ntidentities` is a string giving all of the called nucleotides.
    So `ntindices[0]` gives the value called for `ntstartindex`,
    and `ntindices[j]` gives the value called for `ntstartindex + j`,
    for ``0 <= j < len(ntidentities)``.

    An example:

    >>> aligntup = (14, 26, 'ATAGATACTAGAT', 1, 13, '.....C.....T.', 
    ...       13, 1, '.....C.......')
    >>> (ntstartindex, ntidentities) = ParsePairedAlignment(aligntup, set([]), set([]))
    >>> ntstartindex == 14
    True
    >>> ntidentities == 'ATAGACACTAGNT'
    True

    Examples with excluded nucleotides

    >>> aligntup = (14, 26, 'ATAGATACTAGAT', 1, 13, '.....C.....T.', 
    ...       13, 1, '.....C.......')
    >>> (ntstartindex, ntidentities) = ParsePairedAlignment(aligntup, set([1, 3, 13]), set([1, 3]))
    >>> ntstartindex == 14
    True
    >>> ntidentities == 'aTaGACACTAgNN'
    True

    Example with N nucleotides 

    >>> aligntup = (14, 26, 'ATAGATACTAGAT', 1, 13, 'N.N..C....NT.', 
    ...       13, 1, 'N....C.......')
    >>> (ntstartindex, ntidentities) = ParsePairedAlignment(aligntup, set([]), set([2, 3]))
    >>> ntstartindex == 14
    True
    >>> ntidentities == 'NTaGACACTANtT'
    True

    """
    if _calignimported and use_calign:
        return mapmuts.calign.ParsePairedAlignment(aligntup, r1exclude, r2exclude)
    (gstart, gend, gseq, r1start, r1end, r1seq, r2start, r2end, r2seq)\
            = aligntup
    lena = gend - gstart + 1
    assert lena == len(gseq) == len(r1seq) == len(r2seq)
    if r1start > r1end:
        r1isrc = True
        assert lena == r1start - r1end + 1 == r2end - r2start + 1
    else:
        r1isrc = False
        assert lena == r1end - r1start + 1 == r2start - r2end + 1
    ntidentities = []
    ig = gstart
    ir1 = r1start
    ir2 = r2start
    for i in range(lena):
        nt1 = r1seq[i]
        nt2 = r2seq[i]
        if (ir1 in r1exclude) or nt1 == 'N':
            if (ir2 in r2exclude) or nt2 == 'N':
                ntidentities.append('N')
            elif nt2 == '.':
                ntidentities.append(gseq[i].lower())
            else:
                ntidentities.append(nt2.lower())
        elif (ir2 in r2exclude) or nt2 == 'N':
            if (ir1 in r1exclude) or nt1 == 'N':
                ntidentities.append('N')
            elif nt1 == '.':
                ntidentities.append(gseq[i].lower())
            else:
                ntidentities.append(nt1.lower())
        elif nt1 == nt2:
            if nt1 == '.':
                ntidentities.append(gseq[i])
            else:
                ntidentities.append(nt1)
        else:
            ntidentities.append('N')
        ig += 1
        if r1isrc:
            ir1 -= 1
            ir2 += 1
        else:
            ir1 += 1
            ir2 -= 1
    return (gstart, ''.join(ntidentities))


class ReadPairToGeneAligner(object):
    """Aligns overlapping paired-end reads with known adaptors to gene.

    This class is designed to process data from the following situation.
    We are interested in accurate sequencing of a gene, as from a PCR
    amplicon. That gene, possibly with some extensions on the termini
    (as from primers), is sheared to create short inserts. These reads
    then have known adaptors ligated onto each end, and are then
    sequenced using a paired-end technology (such as Illumina 50 nt PE).
    If the insert length is similar to the read length, the reads
    overlap, providing double coverage of part of the read. This class
    is designed to align the paired reads to identify the adaptors and
    overlapping read, and then align these to the full gene. Certain
    requirements are imposed for the maximal number of mismatches and
    the minimal amount of overlap. Once an aligner object is created, it
    can be used to process many read pairs. It keeps track of the
    alignment statistics for these reads.

    Note that gaps are not handled by the alignment methods used here.
    Any reads with gaps relative to each other, the adaptors, or the
    gene to which they are being aligned will not be handled, and no
    alignment will be returned.

    Also, N / n nucleotides are not treated as mismatches, rather a
    limit is put on the number of allowed N / n nucleotides by the
    initialization argument maxn.

    Note that this method will only be computationally effective for
    relatively short genes (such as no more than a few thousand
    nucleotides), since it uses the `AlignReadToGene` function, which is
    not efficient for very long alignment templates.

    To use an objects from this class, create an aligner as follows::

        aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2,
            maxn, minoverlap, maxrm, maxa1m, maxa2m, maxgenem, upcase,
            use_cext=True)

    The initialization arguments have the following meanings:

    `fullgene` : a string giving the full sequence of the gene that is
    being sequenced. Note that this full sequence might be longer
    than the actual gene range of interest -- this would be the case
    for example if PCR adds extra terminal sequences to the
    sequenced gene. The substring of fullgene that is of interest
    is defined by generange, and is denoted by gene in the remaining
    documentation.

    `generange` : a 2-tuple of the form `(gstart, gend)`. This tuple
    specifies the portion of `fullgene` that contains the actual
    sequence of interest (for example, it might be the coding
    region). Specifically, the region of interest is given by
    ``gene = fullgene[gstart : gend]``, where we have the constraints
    that ``0 <= gstart < gend <= len(fullgene)``.

    `a1` : the adaptor sequence that is at the 3' end of read 1 (`r1`)
    after the target sequence (a string).

    `a2` : the adaptor sequence that is at the 5' end of the reverse
    complement of read 2 (`r2`) before the target sequence.
    Note that in the actual sequencing read it will be at the 3'
    end of the read, but is specified as the reverse complement 
    since `r2` is assumed to be reverse complemented

    `maxn` : the maximum number of N / n nucleotides that a read
    can contain. Read pairs where either read has > `maxn` 
    N / n nucleotides are discarded.

    `minoverlap` -> the minimum amount of overlap between `r1` and `r2`,
    representing their joint coverage of the target sequence (an
    integer >= 1). All of this joint coverage must also align to
    gene.

    `maxrm` : the number of mismatches allowed in the overlap of `r1` and
    `r2` (an integer >= 0). N / n nucleotides do not count as 
    mismatches.

    `maxa1m` : the maximum number of mismatches allowed in the overlap of
    `r1` with `a1` (an integer >= 0). N / n nucleotides do not count as
    mismatches.

    `maxa2m` : the maximum number of mismatches allowed in the overlap of
    `r2` with `a2` (an integer >= 0). N / n nucleotides do not count as
    mismatches.

    `maxgenem` : the maximum number of mismatches allowed in the
    alignment of the coverage region of either read with gene.
    N / n nucleotides do not count as mismatches.

    `upcase` : argument specifying whether we convert all sequences to
    upper case. If `True`, the results will be upper case and
    that the alignment is case-insensitive. If
    `False`, the alignments will be case-sensitive, with 'a' and 'A'
    considered a mismatch, and in this case, the results will be
    returned with the original case. The functions will run
    faster if `upcase` is `False`, since that avoids converting
    sequences to new cases. So if you are able to guarantee
    that all of the sequences that you are examining are the
    same case, set this to `False`.

    `use_cext` : specifies that we use any C extensions available in
    speed computations. This is `True` by default. If you
    set it to `False`, pure Python code will be used, and the
    functions will be much slower.

    Paired-reads are aligned by calling the method `Align` with
    reads `r1` and `r2`. Recall that `r1` is the read 1, while `r2` is actually
    the reverse complement of the read 2. This method does the following:

    1) If either `r1` or `r2` has > `maxn` N / n nucleotides, the pair is
       discarded, and the method returns *(False, 'excess_N')*.

    2) Uses `AlignReads` function to see if `r1` and `r2` can be aligned to each 
       other using the adaptors `a1` and `a2`, with the parameters specified by
       `minoverlap`, `maxrm`, `maxa1m`, and `maxa2m`. If they cannot, returns
       a 2-tuple *(False, 'unpaired')*. Otherwise...

    3) Uses `AlignReadToGene` function to align the coverage (non-adaptor) 
       region of `r1` to `fullgene` or its reverse complement. If it can be aligned
       with <= `maxgenem` mismatches, proceeds. If it cannot be aligned to `fullgene`
       or its reverse complement, returns a 2-tuple *(False, 'unaligned_R1')*. 
       Otherwise...

    4) Checks to make sure that `r2` also aligns to `fullgene` as would be
       expected based on the alignment of the overlapping region of `r1`.
       If it does not, returns 2-tuple *(False, 'unaligned_R2')*. Otherwise...

    5) Examines where the reads align to gene (the portion of `fullgene`
       specified by `generange`). If the alignment is entirely to a 
       portion of `fullgene` that is not in the range that constitutes
       `gene`, returns 2-tuple *(False, 'outside_gene')*. Otherwise, determines
       the region that overlaps with `gene` and...

    6) Returns information about the alignment to `gene` (the portion of
       `fullgene` specified by `generange`) in the form of an alignment
       string. This alignment string consists of three lines, spanning
       the overlapping portion of `r1` and `r2` that aligns to `gene` (note
       that this is `gene` as specified by `generange`, not `fullgene`). The
       alignment strings look like this when printed::

            14 ATAGATACTAGAT 26
             1 .....C.....T. 13
            13 .....C.......  1

       where the first line is the number of the first nucleotide shown
       in `gene` in 1, 2, ... numbering, a space, the `gene` sequence, a
       space, and the number of the last nucleotide shown in `gene`. The
       second line is the same for `r1`, and the third line is the same
       for `r2` but with the numbering reversed so that it is for the
       true read and not the reverse complement passed to this method.

    The class also has a `Statistics` method which keeps track of statistics
    about how many reads align successfully, the overall length of the 
    inserts being sequenced, the overall length of the overlapping
    coverage region, and the number and position of mismatches along
    the read sequence. These are the accumulated statistics for all reads
    aligned with a given object of this class. The returned variable is a
    dictionary with the following strings as keys, each keying the following
    values:

    'nattempted' : integer number of times the `Align` method has been called.

    'nexcessn' : integer number of read pairs discarded because one of
    the reads had too many N / n nucleotides.

    'npaired' : integer number of the number of times the reads used to
    call `Align` could be aligned to each other using `AlignReads` with the
    specified cutoffs.

    'nalignedfullgene' : integer number of times the paired reads could both
    be aligned to `fullgene` with the specified cutoffs.
    
    'nalignedgene' : integer number of times the paired reads could both
    be aligned to `fullgene` with the specified cutoffs, and this alignment
    possessed at least some overlap with the gene specified by `generange`
    to exist in `fullgene`. This is equivalent to the number of times
    that `Align` returned an alignment string, rather than `False`.

    'insertlength' : a dictionary giving the length of the inserts that
    the two reads are reading for all reads where the read pairs align. 
    This is the region between the adaptors. The shortest possible value
    will be `minoverlap`. The largest possible value will be the sum of
    the two read lengths minus `minoverlap`. This is because shorter or 
    longer inserts cannot be aligned with the specified value of
    `minoverlap`. The dictionary is keyed by integers, with values equal
    to the number of calling read pairs (to `Align` method) that have an 
    insert of that length. Insert lengths that are never found in calling
    reads have no keys. This statistic is for all reads that align
    with each other, regardless of whether they align to `fullgene`.

    'r1lengths' : a dictionary giving the length of the `r1` reads that aligns
    to `fullgene` (precedes the `r1` adaptor). The dictionary is keyed by
    integers equal to the number of calling `r1` reads to `Align`
    that have this length overlap with `gene`. Read lengths that are never
    found in the calling reads have no keys. This statistic is only 
    recorded `r1` / `r2` pairs in which both `r1` and `r2` align to gene with
    the specified mismatch cutoffs.

    'r2lengths' : like 'r1lengths', but for the `r2` reads.

    'r1mismatches' : for each `r1` / `r2` pair that aligns with `fullgene`, gives
    the indices in 1, 2, ... numbering of each position in `r1` that 
    differs from `fullgene`. The dictionary is keyed by integers with
    values equal to the number of times that position in the calling
    `r1` reads to `Align` have a mismatch at this position. Positions
    that never have mismatches have no keys. Note that N / n
    nucleotides do not count as mismatches regardless of what they
    are aligned with.

    'r2mismatches' : like 'r1mismatches, but for mismatches between `r2`
    and `fullgene`.

    The `Statistics` method can also be called with the optional argument `key`
    set to one of the aforementioned string keys, in which case the return
    value is the value for that key in the statistics dictionary, rather than
    full dictionary.

    In the example below, the first set of reads align to the gene. The
    second set don't because of a mismatch between the reads.

    >>> fullgene = 'CATATGGGCATAGCATAGAACT'
    >>> generange = (3, 18)
    >>> a1 = 'TACA'
    >>> a2 = 'CATG'
    >>> aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2, 0, 10,
    ...            0, 0, 0, 0, upcase=False)    
    >>> r1 = 'ATATGGGCATTACA'
    >>> r2 = 'CATGATATGGGCAT'
    >>> print aligner.Align(r1, r2)
    1 ATGGGCAT  8
    3 ........ 10
    8 ........  1
    >>> r1 = 'ATATCGGCATTACA'
    >>> r2 = 'CATGATATGGGCAT'
    >>> print aligner.Align(r1, r2)
    (False, 'unpaired')
    >>> aligner.Statistics('nattempted') == 2
    True
    >>> aligner.Statistics('npaired') == 1
    True
    >>> aligner.Statistics('nalignedfullgene') == aligner.Statistics('nalignedgene') == 1
    True
    >>> aligner.Statistics('insertlength') == {10:1}
    True
    >>> aligner.Statistics('r1lengths') == {10:1}
    True
    >>> aligner.Statistics('r2lengths') == {10:1}
    True
    >>> aligner.Statistics('r1mismatches') == {}
    True
    >>> aligner.Statistics('r2mismatches') == {}
    True

    The same reads that failed to align previously now align because we allow
    one mismatch

    >>> aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2, 0, 8,
    ...            1, 0, 0, 1, upcase=False)    
    >>> print aligner.Align(r1, r2)
    1 ATGGGCAT  8
    3 ..C..... 10
    8 ........  1
    >>> aligner.Statistics('nattempted') == aligner.Statistics('npaired') == 1
    True
    >>> aligner.Statistics('nalignedfullgene') == aligner.Statistics('nalignedgene') == 1
    True
    >>> aligner.Statistics('insertlength') == {10:1}
    True
    >>> aligner.Statistics('r1lengths') == {10:1}
    True
    >>> aligner.Statistics('r2lengths') == {10:1}
    True
    >>> aligner.Statistics('r1mismatches') == {5:1}
    True
    >>> aligner.Statistics('r2mismatches') == {}
    True

    Now alignment fails because the overlap is not long enough

    >>> aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2, 0, 11,
    ...            1, 0, 0, 1, upcase=False)    
    >>> print aligner.Align(r1, r2)
    (False, 'unpaired')
    >>> aligner.Statistics('nattempted') == 1
    True
    >>> aligner.Statistics('npaired') == 0
    True
    >>> aligner.Statistics('nalignedfullgene') == aligner.Statistics('nalignedgene') == 0
    True
    >>> aligner.Statistics('insertlength') == {}
    True
    >>> aligner.Statistics('r1lengths') == {}
    True
    >>> aligner.Statistics('r2lengths') == {}
    True
    >>> aligner.Statistics('r1mismatches') == {}
    True
    >>> aligner.Statistics('r2mismatches') == {}
    True

    In the first example, the alignment succeeds. In the second, the reads
    align to each other, but fail to fully align to the gene, so the alignment
    fails (in the second, the regions before the overlap in read 1 don't match
    the gene).

    >>> fullgene = 'CATATGAGCGGCATAGCATAGAACT'
    >>> generange = (3, 21)
    >>> a1 = 'ATGC'
    >>> a2 = 'AGTC'
    >>> r1 = 'ATGAGCGGCATAG'
    >>> r2 = 'CGGCATAGCATAG'
    >>> aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2, 0, 8,
    ...            0, 0, 0, 0, upcase=False)    
    >>> print aligner.Align(r1, r2)
     6 CGGCATAG 13
     6 ........ 13
    13 ........  6
    >>> r1 = 'GTGAGCGGCATAG'
    >>> r2 = 'CGGCATAGCATAG'
    >>> print aligner.Align(r1, r2)
    (False, 'unaligned_R1')
    >>> aligner.Statistics('nattempted') == aligner.Statistics('npaired') == 2
    True
    >>> aligner.Statistics('nalignedfullgene') == aligner.Statistics('nalignedgene') == 1
    True
    >>> aligner.Statistics('insertlength') == {18:2}
    True
    >>> aligner.Statistics('r1lengths') == {13:1}
    True
    >>> aligner.Statistics('r2lengths') == {13:1}
    True
    >>> aligner.Statistics('r1mismatches') == {}
    True
    >>> aligner.Statistics('r2mismatches') == {}
    True

    Read aligns despite mismatches because of tolerant cutoffs.

    >>> aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2, 0, 8,
    ...            1, 0, 0, 2, upcase=False)    
    >>> r1 = 'ATGAGGGGCTTAG'
    >>> r2 = 'GGGCATAGCATAG'
    >>> print aligner.Align(r1, r2)
     6 CGGCATAG 13
     6 G...T... 13
    13 G.......  6
    >>> aligner.Statistics('nattempted') == aligner.Statistics('npaired') == 1
    True
    >>> aligner.Statistics('nalignedfullgene') == aligner.Statistics('nalignedgene') == 1
    True
    >>> aligner.Statistics('insertlength') == {18:1}
    True
    >>> aligner.Statistics('r1lengths') == {13:1}
    True
    >>> aligner.Statistics('r2lengths') == {13:1}
    True
    >>> aligner.Statistics('r1mismatches') == {6:1, 10:1}
    True
    >>> aligner.Statistics('r2mismatches') == {13:1}
    True

    Reads align to gene reverse complement

    >>> fullgene = 'CATATGAGCGGCATAGCATAGAACT'
    >>> generange = (3, 18)
    >>> a1 = 'CTAG'
    >>> a2 = 'GACC'
    >>> aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2, 0, 8,
    ...            1, 0, 0, 1, upcase=False)    
    >>> r1 = 'CTATGCTATG'
    >>> r2 = 'TATGCTATGC'
    >>> print aligner.Align(r1, r2)
     9 CATAGCA 15
    10 .......  4
     2 .......  8
    >>> r1 = 'TTCTATGCTATG'
    >>> r2 = 'TCTATGCTATGC'
    >>> print aligner.Align(r1, r2)
     9 CATAGCA 15
    12 .......  6
     2 .......  8
    >>> aligner.Statistics('nattempted') == aligner.Statistics('npaired') == 2
    True
    >>> aligner.Statistics('nalignedfullgene') == aligner.Statistics('nalignedgene') == 2
    True
    >>> aligner.Statistics('insertlength') == {11:1, 13:1}
    True
    >>> aligner.Statistics('r1lengths') == {10:1, 12:1}
    True
    >>> aligner.Statistics('r2lengths') == {10:1, 12:1}
    True
    >>> aligner.Statistics('r1mismatches') == {}
    True
    >>> aligner.Statistics('r2mismatches') == {}
    True

    Reads align to gene reverse complement

    >>> fullgene = 'CATATGAGCGGCATAGCATAGAACT'
    >>> generange = (3, 18)
    >>> a1 = 'CTAG'
    >>> a2 = 'GACC'
    >>> r1 = 'CCGCTCATA'
    >>> r2 = 'CGCTCATAT'
    >>> aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2, 0, 8,
    ...            0, 0, 0, 0, upcase=False)    
    >>> print aligner.Align(r1, r2)
    1 ATGAGCG 7
    8 ....... 2
    3 ....... 9
    >>> r1 = 'CGCTCATACTA'
    >>> r2 = 'ACCCGCTCATA'
    >>> print aligner.Align(r1, r2)
    1 ATGAGCG 7
    7 ....... 1
    2 ....... 8
    >>> aligner.Statistics('nattempted') == aligner.Statistics('npaired') == 2
    True
    >>> aligner.Statistics('nalignedfullgene') == aligner.Statistics('nalignedgene') == 2
    True
    >>> aligner.Statistics('insertlength') == {10:1, 8:1}
    True
    >>> aligner.Statistics('r1lengths') == {9:1, 8:1}
    True
    >>> aligner.Statistics('r2lengths') == {9:1, 8:1}
    True
    >>> aligner.Statistics('r1mismatches') == {}
    True
    >>> aligner.Statistics('r2mismatches') == {}
    True

    Here is an example in which the reads align, but it is outside
    of generange, and so False is returned. In the second example,
    we extend the reads one nucleotide into the gene, and then get
    a short returned alignment. 

    >>> fullgene = 'CATGACATAGGCATGGGCATTGCAATGTAGAACCATAGCATA'
    >>> generange = (12, 30)
    >>> a1 = ''
    >>> a2 = ''
    >>> aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2, 0, 8,
    ...            0, 0, 0, 0, upcase=False)    
    >>> r1 = 'ACATAGGC'
    >>> r2 = 'ACATAGGC'
    >>> print aligner.Align(r1, r2)
    (False, 'outside_gene')
    >>> r1 = 'ACATAGGCA'
    >>> r2 = 'ACATAGGCA'
    >>> print aligner.Align(r1, r2)
    1 A 1
    9 . 9
    1 . 1
    >>> aligner.Statistics('nattempted') == aligner.Statistics('npaired') == 2
    True
    >>> aligner.Statistics('nalignedfullgene') == 2
    True
    >>> aligner.Statistics('nalignedgene') == 1
    True
    >>> aligner.Statistics('insertlength') == {8:1, 9:1}
    True
    >>> aligner.Statistics('r1lengths') == {8:1, 9:1}
    True
    >>> aligner.Statistics('r2lengths') == {8:1, 9:1}
    True
    >>> aligner.Statistics('r1mismatches') == {}
    True
    >>> aligner.Statistics('r2mismatches') == {}
    True

    Example of how `maxn` determines treatment of N / n nucleotides

    >>> fullgene = 'CATGACATAGGCATGGGCATTGCAATGTAGAACCATAGCATA'
    >>> generange = (12, 30)
    >>> a1 = ''
    >>> a2 = ''
    >>> aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2, 0, 8,
    ...            0, 0, 0, 0, upcase=False)   
    >>> r1 = 'ATGNGCATTG'
    >>> r2 = 'ATGGGCNTTG'
    >>> print aligner.Align(r1, r2)
    (False, 'excess_N')
    >>> aligner.Statistics('nattempted')
    1
    >>> aligner.Statistics('nexcessn')
    1
    >>> aligner = ReadPairToGeneAligner(fullgene, generange, a1, a2, 1, 8,
    ...            0, 0, 0, 0, upcase=False)  
    >>> print aligner.Align(r1, r2)
     1 ATGGGCATTG 10
     1 ...N...... 10
    10 ......N...  1
    >>> aligner.Statistics('nattempted')
    1
    >>> aligner.Statistics('nexcessn')
    0
    >>> aligner.Statistics('r1mismatches') == {}
    True
    >>> aligner.Statistics('r2mismatches') == {}
    True
    """

    def __init__(self, fullgene, generange, a1, a2, maxn, minoverlap,\
            maxrm, maxa1m, maxa2m, maxgenem, upcase, use_cext=True):
        """Initializes a `ReadPairToGeneAligner` object.
        
        The meaning of the calling parameters is described in the main
        docstring for this class.    
        """
        assert isinstance(maxn, int) and maxn >= 0
        self.maxn = maxn
        self.upcase = upcase
        self.use_cext = use_cext
        assert isinstance(fullgene, str)
        assert isinstance(a1, str)
        assert isinstance(a2, str)
        if self.upcase:
            self.fullgene = fullgene.upper()
            self.a1 = a1.upper()
            self.a2 = a2.upper()
        else:
            self.fullgene = fullgene
            self.a1 = a1
            self.a2 = a2
        assert isinstance(generange, tuple) and len(generange) == 2
        assert 0 <= generange[0] < generange[1] <= len(self.fullgene)
        self.generange = generange
        self.gene = self.fullgene[self.generange[0] : self.generange[1]]
        self.fullgene_rc = mapmuts.sequtils.ReverseComplement(
                self.fullgene, use_csequtils=self.use_cext)
        self.gene_rc = mapmuts.sequtils.ReverseComplement(
                self.gene, use_csequtils=self.use_cext)
        assert isinstance(minoverlap, int) and minoverlap > 0
        self.minoverlap = minoverlap
        assert isinstance(maxrm, int) and maxrm >= 0
        self.maxrm = maxrm
        assert isinstance(maxa1m, int) and maxa1m >= 0
        self.maxa1m = maxa1m
        assert isinstance(maxa2m, int) and maxa2m >= 0
        self.maxa2m = maxa2m
        assert isinstance(maxgenem, int) and maxgenem >= 0
        self.maxgenem = maxgenem
        self.statistics = {'nattempted':0,
                           'nexcessn':0,
                           'npaired':0,
                           'nalignedfullgene':0,
                           'nalignedgene':0,
                           'insertlength':{},
                           'r1lengths':{},
                           'r2lengths':{},
                           'r1mismatches':{},
                           'r2mismatches':{},
                          }

    def Align(self, r1, r2):
        """Aligns paired reads, and then aligns their overlap to gene.

        `r1` : the first read (R1) of a paired-end read.

        `r2` : the reverse complement of the second read (R2) of a
        paired-end read.

        Performs the operations detailed in the main docstring for this
        class, using the mismatch and overlap settings specified at
        object initialization. If the alignment fails, returns
        a 2-tuple *(False, explanation)* where *explanation* is a string
        explaining why the alignment filed (see main docstring for this
        class). If the alignment succeeds, 
        returns the alignment string specified in the main
        docstring for this class.

        Also adds information to this alignment to the aligner's
        statistics, which can be accessed via `Statistics`.
        """
        self.statistics['nattempted'] += 1
        if self.upcase:
            r1 = r1.upper()
            r2 = r2.upper()
        if (r1.count('N') + r1.count('n') > self.maxn) or \
                (r2.count('N') + r2.count('n') > self.maxn):
            self.statistics['nexcessn'] += 1
            return (False, 'excess_N')
        a = AlignReads(r1, r2, self.a1, self.a2, self.minoverlap,
                self.maxrm, self.maxa1m, self.maxa2m, upcase=False,
                use_calign=self.use_cext)
        if not a:
            return (False, 'unpaired') # reads do not align
        (r1overlap, r2overlap, nrm, r1a, r2a) = a
        lenr2 = len(r2)
        insertlength = r1overlap[1] + lenr2 - r2overlap[1]
        try:
            self.statistics['insertlength'][insertlength] += 1
        except KeyError:
            self.statistics['insertlength'][insertlength] = 1
        self.statistics['npaired'] += 1
        lenoverlap = r1overlap[1] - r1overlap[0]
        # r1_na and r2_na are reads without adaptors
        (r1_na, r2_na) = (r1[ : r1a], r2[lenr2 - r2a : ])
        a = AlignReadToGene(r1_na, self.fullgene, self.fullgene_rc,\
                self.maxgenem, upcase=False, use_calign=self.use_cext)
        if not a:
            return (False, 'unaligned_R1') # read 1 does not align to gene
        (r1iposfull, r1isrc, r1nmfull) = a
        # check that r2 aligns OK
        r2iposfull = r1iposfull + r1overlap[0]
        if r1isrc:
            r2_mms = ReadToGeneMismatches(self.fullgene_rc, r2_na,\
                    r2iposfull, self.maxgenem, True, use_calign=self.use_cext)
            if r2_mms == False:
                return (False, 'unaligned_R2') # too many mismatches for read 2
            r1_mms = ReadToGeneMismatches(self.fullgene_rc, r1_na,\
                    r1iposfull, self.maxgenem, use_calign=self.use_cext)
        else:
            r2_mms = ReadToGeneMismatches(self.fullgene, r2_na,\
                    r2iposfull, self.maxgenem, True, use_calign=self.use_cext)
            if r2_mms == False:
                return (False, 'unaligned_R2') # too many mismatches for read 2
            r1_mms = ReadToGeneMismatches(self.fullgene, r1_na,\
                    r1iposfull, self.maxgenem, use_calign=self.use_cext)
        self.statistics['nalignedfullgene'] += 1
        try:
            self.statistics['r1lengths'][len(r1_na)] += 1
        except KeyError:
            self.statistics['r1lengths'][len(r1_na)] = 1
        try:
            self.statistics['r2lengths'][len(r2_na)] += 1
        except KeyError:
            self.statistics['r2lengths'][len(r2_na)] = 1
        for i in r1_mms:
            try:
                self.statistics['r1mismatches'][i] += 1
            except KeyError:
                self.statistics['r1mismatches'][i] = 1
        for i in r2_mms:
            try:
                self.statistics['r2mismatches'][i] += 1
            except KeyError:
                self.statistics['r2mismatches'][i] = 1
        # get the gene overlap sequence
        # fullgstart is index where overlap begins in fullgene
        if r1isrc:
            fullgstart = len(self.fullgene) - r1iposfull - lenoverlap - r1overlap[0]
            gstart = fullgstart - self.generange[0]
            if gstart < 0:
                trim3 = abs(gstart)
                gstart = 0
            else:
                trim3 = 0
            fullgend = fullgstart + lenoverlap
            if fullgstart >= self.generange[1] or fullgend <= self.generange[0]:
                return (False, 'outside_gene') # alignment is entirely outside of generange
            gend = fullgend - self.generange[0]
            genelength = self.generange[1] - self.generange[0]
            if gend > genelength:
                trim5 = gend - genelength
                gend = genelength
            else:
                trim5 = 0 # trim from 5' end of reads before reverse-complementing
            geneoverlap = self.gene[gstart : gend]
            gstartlabel = gstart + 1
            gendlabel = gend
            r1start = r1overlap[0] + trim5
            r1end = r1overlap[1] - trim3
            r1endlabel = r1start + 1
            r1startlabel = r1end
            r2start = r2overlap[0] + trim5
            r2end = r2overlap[1] - trim3
            r2startlabel = lenr2 - r2end + 1
            r2endlabel = lenr2 - r2start
            r1string = mapmuts.sequtils.ReverseComplement(\
                    r1[r1start : r1end], use_csequtils=self.use_cext)
            r2string = mapmuts.sequtils.ReverseComplement(\
                    r2[r2start : r2end], use_csequtils=self.use_cext)
        else:
            fullgstart = r1iposfull + r1overlap[0]
            fullgend = fullgstart + lenoverlap
            if fullgstart >= self.generange[1] or fullgend <= self.generange[0]:
                return (False, 'outside_gene') # alignment is entirely outside of generange
            gstart = fullgstart - self.generange[0]
            gend = fullgend - self.generange[0]
            genelength = self.generange[1] - self.generange[0]
            if gstart < 0:
                trim5 = abs(gstart)
                gstart = 0
            else:
                trim5 = 0
            if gend > genelength:
                trim3 = gend - genelength
                gend = genelength
            else:
                trim3 = 0
            geneoverlap = self.gene[gstart : gend]
            gstartlabel = gstart + 1
            gendlabel = gend
            r1start = r1overlap[0] + trim5
            r1end = r1overlap[1] - trim3
            r1startlabel = r1start + 1
            r1endlabel = r1end
            r2start = r2overlap[0] + trim5
            r2end = r2overlap[1] - trim3
            r2startlabel = lenr2 - r2start
            r2endlabel = lenr2 - r2end + 1
            r1string = r1[r1start : r1end]
            r2string = r2[r2start : r2end]
        self.statistics['nalignedgene'] += 1
        r1string = AddDots(geneoverlap, r1string, use_calign=self.use_cext)
        r2string = AddDots(geneoverlap, r2string, use_calign=self.use_cext)
        startdigits = int(math.ceil(math.log10(max([gstartlabel,\
                r1startlabel, r2startlabel])) + 1e-6))
        enddigits = int(math.ceil(math.log10(max([gendlabel,\
                r1endlabel, r2endlabel])) + 1e-6))
        s = "%*d %s %*d\n%*d %s %*d\n%*d %s %*d" % (startdigits,\
                gstartlabel, geneoverlap, enddigits, gendlabel,\
                startdigits, r1startlabel, r1string, enddigits,\
                r1endlabel, startdigits, r2startlabel, r2string,\
                enddigits, r2endlabel)
        return s

    def Statistics(self, key=None):
        """Get statistics about attempted read alignments.

        This method returns statistics about all attempted alignments
        of reads using the current object. The statistics are returned
        in the form of a dictionary with the keys specified in the main
        documentation string for this class. If the optional argument
        `key` is set to a key string in this dictionary, then just returns
        that argument in the statistics dictionary rather than the full 
        dictionary.
        """
        if isinstance(key, str):
            try:
                return self.statistics[key]
            except KeyError:
                raise ValueError('Invalid key to Statistics method: %s' % key)
        else:
            return self.statistics



def AlignReads(r1, r2, a1, a2, minoverlap, maxrm, maxa1m, maxa2m,
        upcase, n_is_mismatch=False, use_calign=True):
    """Aligns overlapping paired-end reads with known adaptors.

    This function rapidly aligns overlapping paired-end reads `r1`
    and `r2`. These reads cover the same template sequence in opposite
    directions with possible overlap, such as would be generated when
    sequencing a short Illumina insert. Gaps in the read overlap are
    currently not handled, so the presence of any gap will prevent
    alignment. 

    This is NOT an alignment method that will necessarily find the
    optimal alignment. Instead, it searches rapidly for the first match
    that meets the specified criteria, and returns this. So be careful
    if you set lenient values for `minoverlap`, `maxrm`, `maxa1m`, and `maxa2m`.
    However, this method is fast!

    `r1` : the first read (a string).

    `r2` : the second read, should already be reverse complemented (a
    string).

    `a1` : the adaptor sequence that is at the 3' end of `r1` after the
    target sequence (a string).

    `a2` : the adaptor sequence that is at the 5' end of `r2` before the
    target sequence. Note that in the actual sequence read it will
    be at the 3' end of the read, but since `r2` is already reverse
    complemented it is now at the 5' end, and itself should be
    reverse complemented from the actual adaptor (a string).

    `minoverlap` : the minimum amount of overlap between `r1` and `r2`,
    representing their joint coverage of the target sequence (an
    integer >= 1).

    `maxrm` : the number of mismatches allowed in the overlap of `r1` and
    `r2` (an integer >= 0).

    `maxa1m` : the maximum number of mismatches allowed in the overlap of
    `r1` with `a1` (an integer >= 0).

    `maxa2m` : the maximum number of mismatches allowed in the overlap of
    `r2` with `a2` (an integer >= 0).

    `upcase` : argument specifying whether we convert all sequences to
    upper case. If `True`, the results will be upper case and
    that the alignment is case-insensitive. If
    `False`, the alignments will be case-sensitive, with 'a' and 'A'
    considered a mismatch, and in this case, the results will be
    returned with the original case. In general, the method will
    run faster if upcase is `False`, since that avoids converting
    sequences to new cases. So if you are able to guarantee
    that all of the sequences that you are examining are the
    same case, set this to `False`.

    `n_is_mismatch` : a Boolean switch specifying that we consider 'N'
    or 'n' nucleotides to be mismatches if they align with something
    else. By default this is `False`, meaning that an 'N' or 'n'
    nucleotide is not considered to be a mismatch regardless of
    what it aligns with. If set to `True`, then an 'N' or 'n' is 
    considered a mismatch unless it aligns with the same character.

    `use_calign` : specifies that we actually perform the calculations
    using the implementation of this function in `calign`. This 
    dramatically speeds the computation. The `calign` extension
    is used assuming that it is available and that this
    parameter is set to its default value of `True`. If `use_calign`
    is set to `False`, or if calign is unavailable, then the 
    computations instead use the pure Python implementation.

    If the reads cannot be aligned with the specified overlap and
    mismatch parameters, returns `False`. If they can be aligned, returns
    the following tuple:
    `(r1overlap, r2overlap, nmismatches, r1a, r2a)`
    where `r1overlap` and `r2overlap` are 2-tuples of integers such that
    ``r1[r1overlap[0] : r1overlap[1]]`` gives the portion of `r1` that overlaps
    with `r2`, and ``r2[r2overlap[0] : r2overlap[1]]`` gives the portion of `r2`
    that overlaps with `r1`. The length of the overlap is therefore
    ``r1overlap[1] - r1overlap[0] = r2overlap[1] - r2overlap[0]``.
    The number of mismatches in the overlapping regions of `r1` and `r2`
    is given by the integer nmismatches. The integers `r1a` and `r2a`
    specify where in the reads the adaptor sequence begins.
    Specifically, `r1a` gives the length of `r1` before adaptor begins (it
    is just `len(r1)` if there is no adaptor). Similarly, `r2a` gives the
    length of `r2` before the adaptor, keeping in mind that `r2` is reverse
    complemented so the length begins at the end of the string and moves
    backward.

    If the read extends longer than the specified adaptor sequence, it
    does not matter if the remaining read sequence past the adaptor can be 
    aligned anything.

    If the settings are not sufficiently stringent to guarantee a unique
    alignment (i.e. `minoverlap` is too small and `maxrm`, `maxa1m`, and `maxa2m`
    are too large), then simply returns the first alignment found that
    meets the specifications. Note that this may not be the optimal
    alignment.

    Example of alignment of sequences with exact overlap, succeeds when
    when minoverlap is 9 but not when it is 10.

        >>> r1 = 'AGACTGAACGAATCGCA'
        >>> r2 = 'GCATAGTCAGACTGAAC'
        >>> a1 = 'GAATCGCAC'
        >>> a2 = 'CGCATAGTC'
        >>> AlignReads(r1, r2, a1, a2, 9, 0, 0, 0, upcase=False)
        ((0, 9), (8, 17), 0, 9, 9)
        >>> r1[0 : 9] == r2[8 : 17]
        True
        >>> AlignReads(r1, r2, a1, a2, 10, 1, 0, 0, upcase=False)
        False

    Another example:

        >>> r1 = 'GTGAGCGGCATAG'
        >>> r2 = 'CGGCATAGCGGTG'
        >>> a1 = ''
        >>> a2 = ''
        >>> AlignReads(r1, r2, a1, a2, 8, 0, 0, 0, upcase=False)
        ((5, 13), (0, 8), 0, 13, 13)

    Example of alignment of sequences with one mismatch in overlap,
    succeeds when maxrm is 1 but not when it is 0. Also succeeds when
    minoverlap is 10, but not when it is 15.

        >>> r1 = 'AGATCAGCATGCAGACTGA'
        >>> r2 = 'GCATAAGACCAGCATGCAG'
        >>> a1 = 'ACTG'
        >>> a2 = 'CATA'
        >>> AlignReads(r1, r2, a1, a2, 10, 1, 0, 0, upcase=False)
        ((0, 14), (5, 19), 1, 14, 14)
        >>> r1[0 : 14] == 'AGATCAGCATGCAG'
        True
        >>> r2[5 : 19] == 'AGACCAGCATGCAG'
        True
        >>> AlignReads(r1, r2, a1, a2, 10, 0, 0, 0, upcase=False)
        False
        >>> AlignReads(r1, r2, a1, a2, 15, 1, 0, 0, upcase=False)
        False

    Example of alignment when there is mismatch in the adaptor
    sequences.

        >>> r1 = 'AGATCAGCATGCAGACTGA'
        >>> r2 = 'GCATAAGACCAGCATGCAG'
        >>> a1 = 'GCTG'
        >>> a2 = 'CTTA'
        >>> AlignReads(r1, r2, a1, a2, 10, 1, 1, 1, upcase=False)
        ((0, 14), (5, 19), 1, 14, 14)
        >>> AlignReads(r1, r2, a1, a2, 10, 1, 1, 0, upcase=False)
        False
        >>> AlignReads(r1, r2, a1, a2, 10, 1, 0, 1, upcase=False)
        False
        >>> AlignReads(r1, r2, a1, a2, 10, 1, 0, 0, upcase=False)
        False

    Another example with a mismatch, this time as the first character.

        >>> r1 = 'CTGTACTAGCAGTCGACTA'
        >>> r2 = 'ATGTACTAGCAGTCGACTA'
        >>> a1 = ''
        >>> a2 = ''
        >>> AlignReads(r1, r2, a1, a2, 10, 1, 0, 0, upcase=False)
        ((0, 19), (0, 19), 1, 19, 19)
        >>> AlignReads(r1, r2, a1, a2, 10, 0, 0, 0, upcase=False)
        False

    Example of alignment of sequences with a gap. The alignment fails
    even though the sequences are otherwise exact matches. Note that
    this is the case even though the gap is the first nucleotide of `r1`

        >>> r1 = 'GACTGAACGAATCGCA'
        >>> r2 = 'GCATAGTCAGACTGAAC'
        >>> a1 = 'GAATCGCAC'
        >>> a2 = 'CGCATAGTC'
        >>> AlignReads(r1, r2, a1, a2, 5, 1, 1, 1, upcase=False)
        False

    Example of use of upcase. Differences in upper and lower case
    prevent alignment if `upcase` is set to `False`, but not when a mismatch
    is allowed.

        >>> r1 = 'AGACTGAAcgaatcgca'
        >>> r2 = 'gcatagtcAGACTGAAC'
        >>> a1 = 'gaatcgcac'
        >>> a2 = 'cgcatagtc'
        >>> AlignReads(r1, r2, a1, a2, 5, 0, 0, 0, upcase=True)
        ((0, 9), (8, 17), 0, 9, 9)
        >>> AlignReads(r1, r2, a1, a2, 5, 0, 0, 0, upcase=False)
        False
        >>> AlignReads(r1, r2, a1, a2, 5, 1, 0, 0, upcase=False)
        ((0, 9), (8, 17), 1, 9, 9)

    Another example of use of `upcase`. Cases match in both reads, so
    alignment is not prevented.

        >>> r1 = 'AGACTGAAcgaatcgca'
        >>> r2 = 'gcatagtcAGACTGAAc'
        >>> a1 = 'gaatcgcac'
        >>> a2 = 'cgcatagtc'
        >>> AlignReads(r1, r2, a1, a2, 5, 0, 0, 0, upcase=True)
        ((0, 9), (8, 17), 0, 9, 9)
        >>> AlignReads(r1, r2, a1, a2, 5, 0, 0, 0, upcase=False)
        ((0, 9), (8, 17), 0, 9, 9)

    Example of `n_is_mismatch` usage:

        >>> r1 = 'NGACTGAACGNATCGCA'
        >>> r2 = 'GCATAGTCAGACTGAAC'
        >>> a1 = 'GAATCGCAC'
        >>> a2 = 'CGCATAGTC'
        >>> AlignReads(r1, r2, a1, a2, 9, 0, 0, 0, upcase=False)
        ((0, 9), (8, 17), 0, 9, 9)
        >>> AlignReads(r1, r2, a1, a2, 9, 0, 0, 0, upcase=False, n_is_mismatch=True)
        False
        >>> AlignReads(r1, r2, a1, a2, 9, 1, 1, 1, upcase=False)
        ((0, 9), (8, 17), 0, 9, 9)
        >>> AlignReads(r1, r2, a1, a2, 9, 1, 1, 1, upcase=False, n_is_mismatch=True)
        ((0, 9), (8, 17), 1, 9, 9)
    """
    if upcase:
        r1 = r1.upper()
        r2 = r2.upper()
        a1 = a1.upper()
        a2 = a2.upper()
    if use_calign and _calignimported:
        return mapmuts.calign.AlignReads(r1, r2, a1, a2, minoverlap, maxrm,\
                maxa1m, maxa2m, n_is_mismatch)
    lenr1 = len(r1)
    lenr2 = len(r2)
    lena1 = len(a1)
    lena2 = len(a2)
    i2max = lenr2 - minoverlap + 1 # value for i1 = 0, then set to 0
    for i1 in range(lenr1 - minoverlap + 1):
        for i2 in range(i2max):
            nrm = 0 # number of mismatches
            for j in range(min(lenr2 - i2, lenr1 - i1)):
                if r1[i1 + j] != r2[i2 + j]:
                    if n_is_mismatch or (r1[i1 + j].upper() != 'N'\
                            and r2[i2 + j].upper() != 'N'):
                        nrm += 1
                        if nrm > maxrm:
                            break # too many mismatches
            else: 
                r1a = i1 + j + 1 # length of r1 before a1 begins
                r2a = lenr2 - i2  # length of r2 before a2 begins
                # found overlap, now check adaptors
                overlap = j + 1
                nma1 = 0
                for k in range(min(lenr1 - r1a, lena1)):
                    if r1[overlap + k] != a1[k]:
                        if n_is_mismatch or (r1[overlap + k].upper() != 'N'\
                                and a1[k].upper() != 'N'):
                            nma1 += 1
                            if nma1 > maxa1m:
                                break # too many mismatches
                else:
                    nma2 = 0
                    for k in range(min(lenr2 - r2a, lena2)):
                        if r2[lenr2 - overlap - 1 - k] != a2[lena2 - 1 - k]:
                            if n_is_mismatch or (r2[lenr2 - overlap - 1 - k].upper()\
                                    != 'N' and a2[lena2 - 1 - k].upper() != 'N'):
                                nma2 += 1
                                if nma2 > maxa2m:
                                    break # too many mismatches
                    else:
                        # we found an alignment
                        assert nrm <= maxrm
                        return ((i1, i1 + overlap), (i2, i2 + overlap), nrm,
                                r1a, r2a)
        i2max = 1 # if not aligning with i1 = 0, must start at i2 = 0
    return False # found no overlap


def AlignReadToGene(r, gene, gene_rc, maxm, upcase, n_is_mismatch=False,\
        use_calign=True):
    """Aligns a read to a gene.

    This function is designed to align a short read to a gene
    sequence. The read can align with either the gene or its
    reverse complement. This is not a true alignment, but rather a
    search to see if the read occurs in the gene with <= some specified
    number of mismatches. If there are multiple ways that the read
    aligns with the gene with <= the specified mismatches, just 
    returns the first (not necessarily the best) one. Also, gaps are
    not allowed, so if there are any gaps the alignment will not be 
    found. Basically, we are searching for the read being a substring
    allowing some number of mismatches.

    `r` : a string giving the read.

    `gene` : a string giving the gene sequence

    `gene_rc` : a string giving the reverse-complement of gene
    so ``gene_rc = sequtils.ReverseComplement(gene)``

    `maxm` : the maximum number of mismatches that are allowed 
    in the alignment of `r` with `gene` or `gene_rc`. If `maxm`
    is not less than `len(r)`, then obviously any position
    will be an alignment, so in general this function is only
    meaningful if `maxm << len(r)`. We only consider an alignment
    to be valid if the number of mismatches is <= `maxm`.

    `upcase` : argument specifying whether we convert all sequences to
    upper case. If `True`, the results will be upper case and
    that the alignment is case-insensitive. If
    `False`, the alignments will be case-sensitive, with 'a' and 'A'
    considered a mismatch, and in this case, the results will be
    returned with the original case. In general, the method will
    run faster if upcase is `False`, since that avoids converting
    sequences to new cases. So if you are able to guarantee
    that all of the sequences that you are examining are the
    same case, set this to `False`.

    `n_is_mismatch` : a Boolean switch specifying that we consider 'N'
    or 'n' nucleotides to be mismatches if they align with something
    else. By default this is `False`, meaning that an 'N' or 'n'
    nucleotide is not considered to be a mismatch regardless of
    what it aligns with. If set to `True`, then an 'N' or 'n' is 
    considered a mismatch unless it aligns with the same character.

    `use_calign` : specifies that we actually perform the calculations
    using the implementation of this function in `calign`. This 
    dramatically speeds the computation. The `calign` extension
    is used assuming that it is available and that this
    parameter is set to its default value of `True`. If `use_calign`
    is set to `False`, or if `calign` is unavailable, then the 
    computations instead use the pure Python implementation.

    If the function is unable to find an alignment of `r` with either
    `gene` or `gene_rc`, then the return value is `False`. If it does align,
    returns a 3-tuple `(ipos, isrc, nm)` where the elements have
    the following meanings:

    `ipos` : the position in `gene` or `gene_rc` to which `r` aligns in
    0, 1, 2, ... numbering. Specifically, `r` aligns with
    ``gene[ipos : ipos + len(r)]`` (if `isrc` is `False`) or with
    ``gene_rc[ipos : ipos + len(r)]`` (if `isrc` is `True`). Since
    `r` must be a substring, we have the constraint that
    ``0 <= ipos < len(gene) - len(r)``.

    `isrc` : a Boolean argument specifying whether the read is
    the reverse complement of the gene. If `False`, it means
    that `r` aligns with `gene`. If `True`, it means that `r` aligns
    with `gene_rc`.

    `nm` : the number of mismatches in the alignment. This will be
    an integer with ``0 <= nm <= maxm``.

    In this example, `r1` aligns exactly to gene, `r2` aligns with one mismatch,
    `r3` does not align because it is not a substring.

    >>> gene = 'ATGATACTAGATAGATATA'
    >>> gene_rc = mapmuts.sequtils.ReverseComplement(gene)
    >>> r1 = 'ATGATACT'
    >>> r2 = 'ATGATTCT'
    >>> r3 = 'CATGATACT'
    >>> AlignReadToGene(r1, gene, gene_rc, 0, upcase=False)
    (0, False, 0)
    >>> AlignReadToGene(r2, gene, gene_rc, 0, upcase=False)
    False
    >>> AlignReadToGene(r2, gene, gene_rc, 1, upcase=False)
    (0, False, 1)
    >>> AlignReadToGene(r3, gene, gene_rc, 1, upcase=False)
    False

    In this examle, two reads that do not align because they are not
    substrings, but they do align when the last nucleotide is removed to
    make them no longer substrings.

    >>> gene = 'ATGATACTAGATAGATATA'
    >>> gene_rc = mapmuts.sequtils.ReverseComplement(gene)
    >>> r1 = 'GATATAT'
    >>> r2 = 'GTATCATA'
    >>> AlignReadToGene(r1, gene, gene_rc, 0, upcase=False)
    False
    >>> AlignReadToGene(r1[ : -1], gene, gene_rc, 0, upcase=False)
    (13, False, 0)
    >>> AlignReadToGene(r2, gene, gene_rc, 0, upcase=False)
    False
    >>> AlignReadToGene(r2[ : -1], gene, gene_rc, 0, upcase=False)
    (12, True, 0)

    In this example, read aligns in internal position

    >>> gene = 'ATGATACTAGATAGATATA'
    >>> gene_rc = mapmuts.sequtils.ReverseComplement(gene)
    >>> r = 'CTAGATAGATA'
    >>> AlignReadToGene(r, gene, gene_rc, 0, upcase=False)
    (6, False, 0)

    In this example, reads do not align due to gaps

    >>> gene = 'ATGATACTAGATAGATATA'
    >>> gene_rc = mapmuts.sequtils.ReverseComplement(gene)
    >>> r1 = 'ATGTACTA'
    >>> r2 = 'TAGATTAGA'
    >>> r3 = 'GATATAT'
    >>> AlignReadToGene(r1, gene, gene_rc, 0, upcase=False)
    False
    >>> AlignReadToGene(r2, gene, gene_rc, 0, upcase=False)
    False
    >>> AlignReadToGene(r3, gene, gene_rc, 0, upcase=False)
    False

    In this example, read aligns in reverse complement

    >>> gene = 'TAGACTAGATAGATAGCTAT'
    >>> gene_rc = mapmuts.sequtils.ReverseComplement(gene)
    >>> r1 = 'ATAGCTATCT'
    >>> AlignReadToGene(r1, gene, gene_rc, 0, upcase=False)
    (0, True, 0)
    >>> r2 = 'TATCTATCTTG'
    >>> AlignReadToGene(r2, gene, gene_rc, 0, upcase=False)
    False
    >>> AlignReadToGene(r2, gene, gene_rc, 1, upcase=False)
    (5, True, 1)

    Examples of `upcase` usage

    >>> gene = 'TGCATGACTAATAGACATA'
    >>> gene_rc = mapmuts.sequtils.ReverseComplement(gene)
    >>> r1 = 'tgcatgact'
    >>> r2 = 'TGCATGAaT'
    >>> AlignReadToGene(r1, gene, gene_rc, 0, upcase=False)
    False
    >>> AlignReadToGene(r1, gene, gene_rc, 0, upcase=True)
    (0, False, 0)
    >>> AlignReadToGene(r2, gene, gene_rc, 0, upcase=False)
    False
    >>> AlignReadToGene(r2, gene, gene_rc, 1, upcase=False)
    (0, False, 1)

    Example of `n_is_mismatch` usage

    >>> gene = 'ATGATACTAGATAGATATA'
    >>> gene_rc = mapmuts.sequtils.ReverseComplement(gene)
    >>> r1 = 'NTGATACT'
    >>> r2 = 'NTGATTCT'
    >>> AlignReadToGene(r1, gene, gene_rc, 0, upcase=False)
    (0, False, 0)
    >>> AlignReadToGene(r2, gene, gene_rc, 0, upcase=False)
    False
    >>> AlignReadToGene(r2, gene, gene_rc, 1, upcase=False)
    (0, False, 1)
    >>> AlignReadToGene(r1, gene, gene_rc, 0, upcase=False, n_is_mismatch=True)
    False
    >>> AlignReadToGene(r1, gene, gene_rc, 1, upcase=False, n_is_mismatch=True)
    (0, False, 1)
    >>> AlignReadToGene(r2, gene, gene_rc, 1, upcase=False, n_is_mismatch=True)
    False
    >>> AlignReadToGene(r2, gene, gene_rc, 2, upcase=False, n_is_mismatch=True)
    (0, False, 2)
    """
    if upcase:
        r = r.upper()
        gene = gene.upper()
        gene_rc = gene_rc.upper()
    if use_calign and _calignimported:
        return mapmuts.calign.AlignReadToGene(r, gene, gene_rc, maxm,\
                n_is_mismatch)
    lenr = len(r)
    lengene = len(gene)
    assert lengene == len(gene_rc)
    for (isrc, g) in [(False, gene), (True, gene_rc)]:
        for i in range(lengene - lenr + 1): # i ranges over gene indices
            nm = 0
            for j in range(lenr): # j ranges of read indices
                if g[i + j] != r[j]:
                    if n_is_mismatch or (g[i + j].upper() != 'N' and\
                            r[j].upper() != 'N'):
                        nm += 1
                        if nm > maxm:
                            break # not a match
            else: # found a match
                return (i, isrc, nm)
    return False
    
    

def Needle(s1, s2, needlecmd, seqtype, outfile, gapopen=10.0,
           gapextend=0.5, outformat='fasta', overwriteoutfile=False, 
           endweight=False):
    """Needleman-Wunsch pairwise alignment using EMBOSS ``needle``.

    This function uses the Needleman-Wunsch algorithm as implemented in
    EMBOSS by the ``needle`` program to perform a pairwise alignment of one
    or more sequences. One target sequence (`s1`) is pairwise aligned with
    one or more other sequences (`s2`). The results are written to a
    specified output file.

    The return value is `None`, but `outfile` is created.

    This function is currently tested with ``needle`` from EMBOSS 6.5.7.0,
    although it should run with many versions.

    `s1` : string giving the first input sequence, to which all sequences
    in `s2` are aligned.

    `s2` : specifies the second input sequence(s) to which we are
    aligning `s1`. It can be:

        * a string giving the second input sequence to align to `s1` (if we
          are aligning just one pair)

        * a list of strings giving all input sequences to align to `s1`
          (if we are making multiple pairwise alignments).
        
        * a 2-tuple specifying that we read `s2` from an existing file. In
          this case, the 2-tuple should be `(filename, fileformat)`
          where `filename` is a string giving the filename, and
          `fileformat` is a string giving the file format (such as
          'fasta' or 'fastq').

    `needlecmd` : string giving executable path for ``needle``. For example,
    this can be a full pathname (such as ``/Users/jbloom/EMBOSS-6.5.7/emboss/needle``) 
    or if ``needle`` is installed in a directory in the current search path it might 
    just be the name of the program (i.e. 'needle').  
    
    `seqtype` : string specifying the sequence type. Can be either 
    'nucleotide' or 'protein'.
    
    `outfile` : string giving the name of the file to which we write the
    output alignment. The format is specified by `outformat`, and
    whether we overwrite any existing outfile is specified by
    `overwriteoutfile`.

    `gapopen` : the gap opening penalty, 10.0 by default.

    `gapextend` : the gap extension penalty, 0.5 by default.

    `outformat` : string specifying the format of the output alignment.
    Is 'fasta' by default.

    `overwriteoutfile` : Boolean switch specifying whether we overwrite
    any existing output files of this name. If `True`, we overwrite
    existing files. If `False`, we raise an `IOError` if there is
    already an existing file with the name of `outfile`.

    `endweight` : the weight applied against gaps at the end of the
    alignment. Is `False` by default, meaning that no end gap penalty
    is applied. If it is set to some other value, it should be the
    2-tuple `(endopen, endextend)` where `endopen` is the penalty for
    opening an end gap, and endextend is the penalty for extending
    an end gap.
    """
    try:
        arguments = ['-gapopen %f -gapextend %f' % (gapopen, gapextend)]
        if seqtype == 'nucleotide':
            arguments.append('-snucleotide1 -snucleotide2')
        elif seqtype == 'protein':
            arguments.append('-sprotein1 -sprotein2')
        else:
            raise ValueError("Invalid value of seqtype: %r" % seqtype)
        (s1fd, s1file) = tempfile.mkstemp() # temporary file
        mapmuts.sequtils.WriteFASTA([('s1', s1)], s1file)
        arguments.append('-sformat1 fasta')
        if isinstance(s2, tuple) and len(s2) == 2:
            s2file = s2[0]
            if not os.path.isfile(s2file):
                raise IOError("Cannot find s2 file of %s" % s2file)
            arguments.append('-sformat2 %s' % s2[1])
        elif isinstance(s2, str):
            (s2fd, s2file) = tempfile.mkstemp() # temporary file
            mapmuts.sequtils.WriteFASTA([('s2', s2)], s2file)
            arguments.append('-sformat2 fasta')
        elif isinstance(s2, list):
            (s2fd, s2file) = tempfile.mkstemp() # temporary file
            mapmuts.sequtils.WriteFASTA(zip(['s2'] * len(s2), s2),
                    s2file)
            arguments.append('-sformat2 fasta')
        else:
            raise ValueError("Invalid argument for s2: %r" % s2)
        if endweight:
            if isinstance(endweight, tuple) and len(endweight) == 2:
                arguments.append('-endweight Y -endopen %f -endextend\
                    %f' % endweight)
            else:
                raise ValueError("Invalid value for endweight: %r" %
                    endweight)
        else:
            arguments.append('-endweight N')
        if os.path.isfile(outfile) and not overwriteoutfile:
            raise IOError("Needle outfile of %s already exists." %
            outfile)
        arguments.append('-aformat3 %s' % outformat)
        (stderrfd, stderrfile) = tempfile.mkstemp()
        (stdoutfd, stdoutfile) = tempfile.mkstemp()
        stdout = open(stdoutfile, 'w')
        stderr = open(stderrfile, 'w')
        cmdstr = '%s %s -asequence %s -bsequence %s -outfile %s' % \
                (needlecmd, ' '.join(arguments), s1file, s2file, outfile)
        returncode = subprocess.call(cmdstr.split(), stdout=stdout,
                stderr=stderr)
        stdout.close()
        stderr.close()
        if returncode != 0:
            raise IOError("Execution of needle failed. \n\nOutput: \n%s\
                \n\n Errors: \n%s" % ((open(stdoutfile).read()),
                open(stderrfile).read()))
    finally:
        # remove any temporary files created by function
        try:
            if os.path.isfile(s1file):
                os.remove(s1file) # remove the temporary file
                os.close(s1fd)
        except:
            pass
        try:
            if not isinstance(s2, tuple) and os.path.isfile(s2file):
                os.remove(s2file) # remove the temporary file
                os.close(s2fd)
        except:
            pass
        try:
            if os.path.isfile(stderrfile):
                os.remove(stderrfile) # remove the temporary file
                os.close(stderrfd)
        except:
            pass
        try:
            if os.path.isfile(stdoutfile):
                os.remove(stdoutfile) # remove the temporary file
                os.close(stdoutfd)
        except:
            pass


if __name__ == '__main__':
    import doctest
    doctest.testmod()
