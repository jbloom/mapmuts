"""Module for reading / writing / processing sequences.

This module is part of the ``mapmuts`` package.

Written by Jesse Bloom.


List of functions
------------------------

`AminoAcids` : returns a list of upper-case one-letter amino-acid codes.

`Codons` : returns a list of all upper-case codon strings.

`OrderedAminoAcids` : amino acids ordered by properties.

`MutAAsCodons` : returns all mutant amino acids and encoding codons.

`NTMutTypes` : returns a list of all possible nucleotide substitution types.

`ClassifyNTMuts` : classifies nucleotide mutation types.

`ClassifyCodonCounts` : classifies codon mutations.

`TallyCodonCounts` : tallies codon mutations for each site.

`MeanQValueFromString` : gets the mean Q-score from a string.

`QValuesFromString` : gets a list of Q-values from a string.

`WriteFASTA` : writes sequences to a FASTA file.

`ReadFASTA` : reads sequences from a FASTA file.

`GetSequence` : gets a particular sequence from a list based on header.

`Translate` : translates nucleotide sequences to proteins.

`AmbiguousNTCodes` : returns all possible nucleotides for ambiguous code.

`PurgeAmbiguousDNA` : removes sequences with ambiguous nucleotides.

`GetEntries` : gets entries from a potentially very large FASTA file.

`AAThreeToOne` : converts three-letter to one-letter amino acid code.

`ReverseComplement` : returns the reverse-complement of a nucleotide sequence.

`FindMotifs` : finds specific motifs in a nucleotide sequence

`CondenseSeqs` : removes nearly identical sequences.

`KyteDoolittle` : gets the Kyte-Doolittle hydrophobicity scale.


List of C-extensions
-----------------------
Some of the functions in this module are accelerated by using the 
C extensions in `csequtils` if that C extension can be imported. This is
done by default, and can be stopped by calling the relevant
function with `use_csequtils` set to `False`. Note that you do not need to
explicitly call the `csequtils` functions -- they are automatically called
by the functions in this module when available. Currently extensions are
available for:

* `ReverseComplement`

* `MeanQValueFromString`


Documentation for individual functions
---------------------------------------
Documentation for individual functions is provided in their definitions
below.

"""


import os
import re
import warnings
import mapmuts.stats
try:
    import mapmuts.csequtils
    _csequtilsimported = True # global variable, is csequtils available?
except ImportError:
    warnings.warn("Cannot import csequtils. Will have to use pure Python"\
            + ' implementations.', RuntimeWarning)
    _csequtilsimported = False


def AminoAcids(includestop=False):
    """Returns a list of all one-letter amino acid codes, upper case.

    *includestop* is an optional argument that is *False* by default. If set
    to *True*, then we append the character * (for a stop codon) at the end of 
    list.

    EXAMPLES:

    >>> AminoAcids()
    ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    >>> len(AminoAcids()) == 20
    True

    >>> AminoAcids(includestop=True)
    ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*']

    >>> len(AminoAcids(includestop=True)) == 21
    True

    """
    aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    if includestop:
        aas.append('*')
        return aas
    else:
        return aas


def Codons():
    """Returns a list of all codons as upper-case strings in alphabetical order.

    >>> codons = Codons()
    >>> len(codons) == 64
    True
    >>> codons[ : 4]
    ['AAA', 'AAC', 'AAG', 'AAT']
    """
    codons = []
    nts = ['A', 'C', 'G', 'T']
    for nt1 in nts:
        for nt2 in nts:
            for nt3 in nts:
                codons.append("%s%s%s" % (nt1, nt2, nt3))
    return codons


def OrderedAminoAcids(ordering):
    """Returns the amino acid codes ordered according to physical properties.

    This function returns a list of all 20 amino-acid one-letter upper-case
    codes ordered according to physical properties.

    The parameter *ordering* specifies how they are ordered. Possible values
    are the following strings:

        * *icelogo* : the ordering used in the IceLogo heat maps
          as at http://iomics.ugent.be/icelogoserver/manual.pdf

        * *hydrophobicity* : ordering from most to least hydrophobic
          using the Kyte-Doolittle scale.         

    Examples:

    >>> aas = OrderedAminoAcids('icelogo')
    >>> aas == ['G', 'A', 'V', 'S', 'T', 'C', 'M', 'L', 'I', 'K', 'R', 'E', 'D', 'Q','N', 'F', 'Y', 'W', 'P', 'H']
    True
    >>> len(aas) == 20
    True

    >>> aas = OrderedAminoAcids('hydrophobicity')
    >>> aas == ['I', 'V', 'L', 'F', 'C', 'M', 'A', 'W', 'G', 'T', 'S', 'Y', 'P', 'H', 'N', 'D', 'Q', 'E', 'K', 'R']
    True
    >>> len(aas) == 20
    True

    """
    if ordering == 'icelogo':
        return ['G', 'A', 'V', 'S', 'T', 'C', 'M', 'L', 'I', 'K', 'R', 'E', 'D', 'Q', 'N', 'F', 'Y', 'W', 'P', 'H']
    elif ordering == 'hydrophobicity':
        return ['I', 'V', 'L', 'F', 'C', 'M', 'A', 'W', 'G', 'T', 'S', 'Y', 'P', 'H', 'N', 'D', 'Q', 'E', 'K', 'R']
    else:
        raise ValueError("Invalid ordering of %s" % ordering)


def MutAAsCodons(wtcodon):
    """Returns a list of all mutant amino acids and their codons.

    *wtcodon* should be a string giving a codon, such as GGG.

    The return variable is a list of twenty 2-tuples. The first entry
    is one of the 19 mutant amino acids relative to *wtcodon* or the stop
    codon. The second entry is a list of all codons that encode for that amino
    acid. For example, one entry might be *('A', ['GCT', 'GCC', 'GCA', 'GCG'])*.

    The * character is used to indicate stop codons. One-letter codes are used
    for amino acids. All strings are converted to upper case.

    This function is coded very inefficiently, so that might concern you if you
    start using a program that really it a very large number of times.
    """
    genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
         'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
         'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
         'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
         'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
         'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
         'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
         'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
         'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W', 'CGT':'R',
         'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
         'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    assert len(genetic_code) == 64
    aas_to_codons = {}
    for (codon, aa) in genetic_code.iteritems():
        if aa in aas_to_codons:
            aas_to_codons[aa].append(codon)
        else:
            aas_to_codons[aa] = [codon]
    assert len(aas_to_codons) == 21
    wtaa = genetic_code[wtcodon.upper()]
    returnlist = [(aa, codons) for (aa, codons) in aas_to_codons.iteritems() if aa != wtaa]
    assert len(returnlist) == 20
    return returnlist


def NTMutTypes():
    """Returns a list of all possible nucleotide substitution types.

    If mutations are made on both strands of double-stranded DNA, there
    are six unique types of nucleotide mutations. For example, A->G, T->C
    are one type (they represent the same mutation on different strands).
    Returns a list of all such types as 2-tuples with two-letter upper-case
    strings of the wildtype and mutant nucleotide, for example: ('AG', 'TC').

    EXAMPLE:

    >>> NTMutTypes()
    [('AT', 'TA'), ('AC', 'TG'), ('AG', 'TC'), ('GA', 'CT'), ('GC', 'CG'), ('GT', 'CA')]
    """
    return [('AT', 'TA'), ('AC', 'TG'), ('AG', 'TC'), ('GA', 'CT'), ('GC', 'CG'), ('GT', 'CA')]


def ClassifyNTMuts(nt_counts, recompute=False):
    """Classifies nucleotide mutation types.

    There are six unique types of nucleotide mutations, as returned by
    `NTMutTypes()`. At each nucleotide position, and also over the whole gene,
    this function counts the number of occurrences of each of these nucleotide
    mutation types. Keys giving these counts are added to the calling dictionary
    `nt_counts`.

    CALLING VARIABLES:

    `nt_counts` : dictionary of the type returned by `mapmuts.io.ReadNTCounts`.

    `recompute` : Boolean switch specifying that we recompute
    any classification values that are already present. Is
    `False` by default, meaning that values are not recomputed.
    If set to `True`, they are recomputed. When checking if values
    need to be recomputed, we only look for the first key listed
    below ``('AT', 'TA')`` in the first nucleotide value (`nt_counts[1]`)
    and if it is present and `recompute` is `True`, nothing further is done.


    RESULT OF THIS FUNCTION:

    The following keys are added both to the main dictionary `nt_counts` and
    to the dictionary `nt_counts[i]` for all `i`: a key corresponding to each
    mutation type returned by `ClassifyNTMuts()`. These keys give the number of
    occurrences of that particular mutation at that nucleotide (for the 
    `codon_counts[i]` keys) for over the whole gene (for the `codon_counts` keys).
    Only nucleotide identities specified by upper-case nucleotides are 
    considered, since these are the only ones that represent definitively
    called identities. In addition, a key 'PAIRED_COUNTS' is added both
    to `nt_counts` and `nt_counts[i]` that gives the number of definitively
    called (upper case) nucleotide calls at each position and for all
    positions.

    EXAMPLES:

    >>> nt_counts = {1:{'WT':'A', 'COUNTS':7, 'A':3, 'T':1, 'C':0, 'G':0,
    ...      'a':1, 't':1, 'c':0, 'g':0, 'N':1}, 2:{'WT':'C', 'COUNTS':8,
    ...      'A':0, 'T':1, 'C':4, 'G':1, 'a':1, 't':0, 'c':1, 'g':0, 'N':0}}
    >>> ClassifyNTMuts(nt_counts)

    >>> nt_counts['PAIRED_COUNTS'] == 10
    True
    >>> nt_counts[('AT', 'TA')] == 1
    True
    >>> nt_counts[('AC', 'TG')] == 0
    True
    >>> nt_counts[('AG', 'TC')] == 0
    True
    >>> nt_counts[('GA', 'CT')] == 1
    True
    >>> nt_counts[('GT', 'CA')] == 0
    True
    >>> nt_counts[('GC', 'CG')] == 1
    True

    >>> nt_counts[1]['PAIRED_COUNTS'] == 4
    True
    >>> nt_counts[1][('AT', 'TA')] == 1 
    True
    >>> nt_counts[1][('AC', 'TG')] == 0
    True
    >>> nt_counts[1][('AG', 'TC')] == 0 
    True
    >>> nt_counts[1][('GA', 'CT')] == 0
    True
    >>> nt_counts[1][('GT', 'CA')] == 0
    True
    >>> nt_counts[1][('GC', 'CG')] == 0
    True

    >>> nt_counts[2]['PAIRED_COUNTS'] == 6
    True
    >>> nt_counts[2][('AT', 'TA')] == 0
    True
    >>> nt_counts[2][('AC', 'TG')] == 0
    True
    >>> nt_counts[2][('AG', 'TC')] == 0
    True
    >>> nt_counts[2][('GA', 'CT')] == 1
    True
    >>> nt_counts[2][('GT', 'CA')] == 0
    True
    >>> nt_counts[2][('GC', 'CG')] == 1
    True

    """
    if not nt_counts:
        raise ValueError("empty nt_counts")
    if 'PAIRED_COUNTS' in nt_counts[1]:
        if not recompute:
            return # already computed
    maxnt = max([i for i in nt_counts.keys() if isinstance(i, int)])
    ntmuttypes = NTMutTypes()
    keys = ['PAIRED_COUNTS'] + ntmuttypes
    for key in keys:
        nt_counts[key] = 0
        for nt in range(1, maxnt + 1):
            nt_counts[nt][key] = 0
    for nt in range(1, maxnt + 1):
        wt = nt_counts[nt]['WT']
        n = nt_counts[nt][wt]
        nt_counts['PAIRED_COUNTS'] += n
        nt_counts[nt]['PAIRED_COUNTS'] += n
        muts = [x for x in ['A', 'T', 'C', 'G'] if x != wt]
        assert len(muts) == 3, str(muts)
        for mut in muts:
            n = nt_counts[nt][mut]
            if n:
                nt_counts['PAIRED_COUNTS'] += n
                nt_counts[nt]['PAIRED_COUNTS'] += n
                for tup in ntmuttypes:
                    x = "%s%s" % (wt, mut)
                    if x == tup[0] or x == tup[1]:
                        break
                else:
                    raise ValueError("Failed to find mutation type %s" % x)
                nt_counts[nt][tup] += n
                nt_counts[tup] += n


def TallyCodonCounts(codon_counts, sites='all'):
    """Tallies codon mutations at each site.

    *codon_counts* : a list with entries being dictionaries of the type returned by
    *mapmuts.io.ReadCodonCounts*. There must be at least one such dictionary 
    specified. If there are multiple dictionaries, they must all specify the
    same set of sites as their keys.

    *sites* is an optional argument that is the string *all* by default.
    If you set it to some other value, then it should be a list of 
    integer sites. In this case, we only tally sites for the sites in this
    list (which all must be keys in the *codon_counts* dictionaries).

    For all sites *r* represented in *codon_counts* and for which all
    of the entries of *codon_counts* specify the same wildtype residue,
    this function iterates
    through all possible codon mutations (so there are *63 L* of these where *L*
    is the number of codons. It records the number of definitely called
    (all uppercase nucleotides) occurrences of each of these mutations. It then
    returns the following tuple:
    *(all_counts, multi_nt_all_counts, syn_counts, multi_nt_syn_counts)*

    This tuple has the following elements:

        * *all_counts* is a list of *63 L* integers giving the number of
          occurrences of all *63 L* mutations. The order is arbitrary.

        * *multi_nt_all_counts* is like *all_counts* except that it just
          includes counts for the subset of all mutations that involve
          more than one nucleotide change to the codon.

        * *syn_counts* is like *all_counts* except it just includes
          counts for the subset of mutations that are synonymous.

        * *multi_nt_syn_counts* is like *all_counts* except that it just
          includes counts for the subset of mutations that are synonymous
          and involve multiple nucleotide changes to the codon.
    """
    if not (isinstance(codon_counts, list) and len(codon_counts) >= 1):
        raise ValueError("codon_counts must be a list containing at least one entry.")
    if sites == 'all':
        sites = codon_counts[0].keys()
    else:
        for r in sites:
            if r not in codon_counts[0].keys():
                raise ValueError("sites specifies a site that has no information: %s" % r)
    all_counts = []
    multi_nt_all_counts = []
    syn_counts = []
    multi_nt_syn_counts = []
    codons = mapmuts.sequtils.Codons()
    aas = [tup[1] for tup in Translate([('codon', codon) for codon in codons])]
    for r in sites:
        wt = codon_counts[0][r]['WT']
        for i_d in codon_counts[1 : ]:
            if i_d[r]['WT'] != wt:
                continue # not all lists have same wildtype here
        wtaa = Translate([('WT', wt)])[0][1]
        for (codon, aa) in zip(codons, aas):
            if codon != wt:
                ndiffs = len([i for i in range(len(codon)) if codon[i] != wt[i]])
                all_n = multi_nt_all_n = syn_n = multi_nt_syn_n = 0
                for i_d in codon_counts:
                    all_n += i_d[r][codon]
                    if ndiffs > 1:
                        multi_nt_all_n += i_d[r][codon]
                    if aa == wtaa:
                        syn_n += i_d[r][codon]
                        if ndiffs > 1:
                            multi_nt_syn_n += i_d[r][codon]
                all_counts.append(all_n)
                if ndiffs > 1:
                    multi_nt_all_counts.append(multi_nt_all_n)
                if aa == wtaa:
                    syn_counts.append(syn_n)
                    if ndiffs > 1:
                        multi_nt_syn_counts.append(multi_nt_syn_n)
    assert len(all_counts) >= len(syn_counts) >= len(multi_nt_syn_counts)
    assert len(all_counts) >= len(multi_nt_all_counts)
    return (all_counts, multi_nt_all_counts, syn_counts, multi_nt_syn_counts)


def ClassifyCodonCounts(codon_counts, recompute=False):
    """Classifies codon mutations.

    CALLING VARIABLES:

    `codon_counts` : dictionary as returned by `mapmuts.io.ReadCodonCounts`

    `recompute` : a Boolean switch specifying that we recompute
    any classification values that are already present. Is
    `False` by default, meaning that values are not recomputed.
    If set to `True`, they are recomputed. When checking if values
    need to be recomputed, we only look for the first key listed
    below ('PAIRED_COUNTS') in the first codon value (`codon_counts[1]`)
    and if it is present and recompute is `False`, nothing further is done.

    RESULT OF CALLING THIS FUNCTION:

    The following keys are added to `codon_counts[i]` for all `i` if new
    values are being computed. Here are the added keys:

    'PAIRED_COUNTS' : total number of counts of codons that contain
    all nucleotides that were called by both reads (all nucleotides
    upper case, and not N).

    'SINGLE_COUNTS' : total number of counts of codons that contain at
    least one nucleotide that was only called by one read (at least
    one nucleotide lower case, but no N nucleotides).

    'AMBIGUOUS_COUNTS' : total number of counts of codons that contain
    at least one ambiguous nucleotide (containing an 'N').

    'N_WT' : total number of definitively called codons (all
    nucleotides upper case) that match 'WT', the wildtpe codon
    at this position.

    'N_NS' : total number of definitively called codons (all nucleotides
    upper case) that represent non-synonymous mutations from the
    wildtype codon. Mutations from stop codons to non-stop codons
    are classified as nonsynonymous. Mutations from non-stop codons
    to stop codons are not classified as nonsynonymous, instead look
    under 'N_STOP'.

    'N_SYN' : total number of definitively called codons (all nucleotides
    upper case) that represent synonymous mutations from the wildtype
    codon. Mutations from stop codons to other stop codons are
    classified as synonymous.

    'N_1MUT', 'N_2MUT', 'N_3MUT' : the total number of definitively called
    mutant codons (all nucleotides upper case) that contain one, two, or
    three nucleotide mutations relative to the wildtype codon at the site.

    ``'N_%s' % aa`` for all ``aa in AminoAcids()`` : total number of definitively
    called codons (all nucleotides upper case) that encode each amino
    acid in `AminoAcids()`.
    
    'N_STOP' : total number of definitively called codons (all nucleotides
    upper case) that encode stop codons. 

    In addition, the following keys are added directly to `codon_counts`, and
    represent totals over all codon positions:

    'TOTAL_COUNTS' : total number of paired counts for all codons.

    'TOTAL_NS' : total number of 'N_NS' for all codons.

    'TOTAL_SYN' : total number of 'N_SYN' for all codons.

    'TOTAL_STOP' : total number of 'N_STOP' for all codons.

    'TOTAL_MUT' : sum of 'N_NS' + 'N_SYN' + 'N_STOP', total number of mutations.

    'TOTAL_N_1MUT', 'TOTAL_N_2MUT', 'TOTAL_N_3MUT' : total number of
    'N_1MUT', N_2MUT', and 'N_3MUT' for all codons.
    """
    if (not codon_counts) or ('PAIRED_COUNTS' in codon_counts[1] and not recompute):
        return codon_counts # nothing to be added to codon_counts
    codon_counts['TOTAL_COUNTS'] = 0
    codon_counts['TOTAL_NS'] = 0
    codon_counts['TOTAL_SYN'] = 0
    codon_counts['TOTAL_STOP'] = 0
    codon_counts['TOTAL_N_1MUT'] = 0
    codon_counts['TOTAL_N_2MUT'] = 0
    codon_counts['TOTAL_N_3MUT'] = 0
    aminoacids = AminoAcids()
    allcodons = mapmuts.io.AllCodons()
    keys = ['PAIRED_COUNTS', 'SINGLE_COUNTS', 'AMBIGUOUS_COUNTS',\
            'N_WT', 'N_NS', 'N_SYN', 'N_STOP'] + ['N_%s' % aa for \
            aa in aminoacids] + ['N_1MUT', 'N_2MUT', 'N_3MUT']
    maxi = max([i for i in codon_counts.keys() if isinstance(i, int)])
    for i in range(1, maxi + 1):
        di = codon_counts[i]
        for key in keys:
            di[key] = 0 # initialize all counts to zero
        wtcodon = di['WT']
        assert wtcodon == wtcodon.upper() and len(wtcodon) == 3
        wtaa = Translate([('head', wtcodon)])[0][1]
        if not wtaa:
            wtaa = 'STOP'
        else:
            assert wtaa in aminoacids
            pass
        for codon in allcodons:
            n = di[codon]
            if n:
                if 'N' in codon:
                    di['AMBIGUOUS_COUNTS'] += n
                elif codon == codon.upper():
                    di['PAIRED_COUNTS'] += n
                    codon_counts['TOTAL_COUNTS'] += n
                    if wtcodon == codon:
                        di['N_WT'] += n
                        di['N_%s' % wtaa] += n
                    else:
                        aa = Translate([('head', codon)])[0][1]
                        if not aa:
                            aa = 'STOP'
                        di['N_%s' % aa] += n
                        if wtaa == aa:
                            di['N_SYN'] += n
                            codon_counts['TOTAL_SYN'] += n
                        elif aa != 'STOP':
                            di['N_NS'] += n
                            codon_counts['TOTAL_NS'] += n
                        elif aa == 'STOP':
                            codon_counts['TOTAL_STOP'] += n
                        ndiffs = len([j for j in range(3) if wtcodon[j] != codon[j]])
                        di['N_%dMUT' % ndiffs] += n
                        codon_counts['TOTAL_N_%dMUT' % ndiffs] += n
                else:
                    di['SINGLE_COUNTS'] += n
        assert di['COUNTS'] == di['PAIRED_COUNTS'] + di['SINGLE_COUNTS'] + di['AMBIGUOUS_COUNTS']
        if wtaa != 'STOP':
            assert di['PAIRED_COUNTS'] == di['N_WT'] + di['N_NS'] + di['N_SYN'] + di['N_STOP'], "%d %d %d %d %d" % (di['N_WT'], di['PAIRED_COUNTS'], di['N_NS'], di['N_SYN'], di['N_STOP'])
            pass
    codon_counts['TOTAL_MUT'] = codon_counts['TOTAL_SYN'] + codon_counts['TOTAL_NS'] + codon_counts['TOTAL_STOP']


def MeanQValueFromString(qstring, use_csequtils=True):
    """Gets mean Q-value from string.

    CALLING VARIABLES:

    `qstring` : string giving Q-values, encoded by single characters
    according to the Sanger format (this is the same format used by
    Illumina 1.8). A character is converted to its numeric Q-value
    using ``q = ord(qchar) - 33``

    `use_csequtils` : a Boolean switch specifying that we perform the
    calculation using the C-extension, which will be much faster.
    `True` by default. If you do not use this extension, this function
    is not any faster than
    `mapmuts.stats.Mean(QValuesFromString(qstring))`

    RETURN VARIABLE:

    The return variable is the mean of all of the Q-values in the
    string.

    EXAMPLE:

    >>> print "%.3f" % MeanQValueFromString('!"#;@Ih')
    24.429
    """
    if _csequtilsimported and use_csequtils:
        return mapmuts.csequtils.MeanQValueFromString(qstring)
    return mapmuts.stats.Mean(QValuesFromString(qstring))


def QValuesFromString(qstring):
    """Gets list of Q-values from string.

    CALLING VARIABLE:

    `qstring` : string of Q-values encoded by single characters
    according to the Sanger format (this is the same format used by
    Illumina 1.8). A character is converted to its numeric Q-value
    using:
    ``q = ord(qchar) - 33``

    RETURN VARIABLE:

    The return variable is a list of Q-values as integers.

    EXAMPLE:

    >>> QValuesFromString('!"#;@Ih')
    [0, 1, 2, 26, 31, 40, 71]
    """
    return [ord(char) - 33 for char in qstring]


def WriteFASTA(headers_seqs, filename, writable_file=False):
    """Writes sequences to a FASTA file.

    CALLING VARIABLES:

    `headers_seqs` : list of 2-tuples specifying sequences and their
    corresponding headers.  Each entry is the 2-tuple `(header, seq)`
    where `header` is a string giving the header (without the leading ">"),
    and `seq` is the corresponding sequence.

    `filename` : string that specifies the name of the file to which the
     headers and sequences should be written.  If this file already exists,
     it is overwritten. 

    `writable_file` : Boolean switch specifying that rather than `filename`
    giving a string specifying the name of a file to which the sequences
    should be written, it instead specifies a writable file object to which
    the sequences should be written.

    RESULT OF THIS FUNCTION:

    The sequences are written to the file in the same order that they are specified
    in `headers_seqs`.
    """
    assert isinstance(writable_file, bool)
    if writable_file:
        f = filename
    else:
        f = open(filename, 'w')
    for (header, seq) in headers_seqs:
        f.write(">%s\n%s\n" % (header, seq))
    if not writable_file:
        f.close()


def ReadFASTA(fastafile):
    """Reads sequences from a FASTA file.

    CALLING VARIABLE:

    `fastafile` : specify the name of a FASTA file.

    RETURN VARIABLE: 
    
    This function reads all sequences from the FASTA file.  It returns the
    list `headers_seqs`.  This list is composed of a 2-tuple `(header, seq)`
    for every sequence entry in `fastafile`. `header` is the header for
    a sequence, with the leading ">" and any trailing spaces removed. `seq`
    is the corresponding sequence.
    """
    lines = open(fastafile).readlines()
    headers_seqs = []
    header = None
    seq = []
    for line in lines:
        if line[0] == '>':
            if (not header) and (not seq):
                pass # first sequence in file
            elif header and not seq:
                raise ValueError, "Empty sequence for %s" % header
            elif seq and not header:
                raise ValueError, "File does not begin with header."
            else:
                seq = ''.join(seq)
                seq = seq.replace(' ', '')
                headers_seqs.append((header, seq))
            header = line.strip()[1 : ]
            seq = []
        else:
            seq.append(line.strip())
    if (not header) and (not seq):
        pass # first sequence in file
    elif header and not seq:
        raise ValueError, "Empty sequence for %s" % header
    elif seq and not header:
        raise ValueError, "File does not begin with header."
    else:
        seq = ''.join(seq)
        seq = seq.replace(' ', '')
        headers_seqs.append((header, seq))
    return headers_seqs


def GetSequence(header, headers_sequences):
    """Gets a particular sequence based on its header name.

    `header` : the name of a sequence's FASTA header.

    `headers_sequences` : list of tuples `(head, seq)` as would
    be returned by `Read`.

    This function searches through `headers_sequences` and returns the
    sequence corresponding to the first header found that matches
    the calling argument `header`.  If no such header is found,
    raises an exception.
    """
    for (head, seq) in headers_sequences:
        if head == header:
            return seq
    else:
        raise ValueError, "Could not find a header matching %s" % header


def Translate(headers_sequences, readthrough_n=False, readthrough_stop=False, truncate_incomplete=False, translate_gaps=False):
    """Translates a set of nucleotide sequences to amino acid sequences.

    CALLING VARIABLES:

    `headers_sequences` : list of tuples `(header, seq)` as would be returned
    by `Read`.  The sequences should all specify valid coding nucleotide
    sequences.  
    
    The returned variable is a new list in which
    all of the nucleotide sequences have been translated to their
    corresponding protein sequences, given by one letter codes.
    If any of the nucleotide sequences do not translate to
    valid protein sequences, an exception is raised.

    `readthrough_n` : specifies that if any nucleotides
    in the sequence are equal to to an ambiguous nt code and cannot therefore
    be unambiguously translated into an amino acid, we simply translate through these 
    nucleotides by making the corresponding amino acid equal to "X".  By
    default, this option is `False`.  Note that even when this option is `False`,
    certain ambiguous nucleotides may still be translatable if they all lead to 
    the same amino acid.

    `readthrough_stop` : specifies that if we encounter any stop
    we simply translation them to 'X'.  By default,
    this option is `False`, meaning that we instead raise an error
    of an incomplete stop codon.

    `truncate_incomplete` : specifies that if the sequence
    length is not a multiple of three, we simply truncate off the one or two
    final nucleotides to make the length a multiple of three prior to translation.
    By default, this option is `False`, meaning that no such truncation is done.

    `translate_gaps` : specifies that a codon of '---' is translated
    to '-'. Codons with one '-' are also translated to gaps.

    RETURN VARIABLE:

    The returned variable is a new list in which
    all of the nucleotide sequences have been translated to their
    corresponding protein sequences, given by one letter codes.
    If any of the nucleotide sequences do not translate to
    valid protein sequences, an exception is raised.

    EXAMPLES:

    >>> Translate([('seq1', 'ATGTAA'), ('seq2', 'gggtgc')])
    [('seq1', 'M'), ('seq2', 'GC')]

    >>> Translate([('seq2', 'GGNTGC')])
    [('seq2', 'GC')]
    
    >>> Translate([('seq2', 'NGGTGC')])
    Traceback (most recent call last):
       ...
    ValueError: Cannot translate codon NGG
    
    >>> Translate([('seq2', 'NGGTGC')], readthrough_n=True)
    [('seq2', 'XC')]

    >>> Translate([('seq2', 'TAATGC')])
    Traceback (most recent call last):
       ...
    ValueError: Premature stop codon

    >>> Translate([('seq2', 'TAATGC')], readthrough_stop=True)
    [('seq2', 'XC')]

    >>> Translate([('seq2', 'TGCA')])
    Traceback (most recent call last):
       ...
    ValueError: Sequence length is not a multiple of three

    >>> Translate([('seq2', 'TGCA')], truncate_incomplete=True)
    [('seq2', 'C')]

    >>> Translate([('seq2', 'TGC---')])
    Traceback (most recent call last):
       ...
    ValueError: Cannot translate gap.

    >>> Translate([('seq2', 'TGC---')], translate_gaps=True)
    [('seq2', 'C-')]
    """
    genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
            'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
            'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
            'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
            'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
            'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'STOP', 'TAG':'STOP',
            'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
            'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
            'TGT':'C', 'TGC':'C', 'TGA':'STOP', 'TGG':'W', 'CGT':'R',
            'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
            'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    assert isinstance(headers_sequences, list)
    translated_headers_sequences = []
    for (head, seq) in headers_sequences:
        seq = seq.upper()
        if len(seq) % 3:
            if truncate_incomplete:
                seq = seq[ : -(len(seq) % 3)]
            else:
                raise ValueError, "Sequence length is not a multiple of three"
        prot_length = len(seq) // 3
        prot = []
        for i in range(prot_length):
            codon = seq[3 * i : 3 * (i + 1)]
            try:
                aa = genetic_code[codon]
            except KeyError:
                if '-' in codon:
                    if translate_gaps:
                        aa = '-'
                    else:
                        raise ValueError("Cannot translate gap.")
                else:
                    # see if we have an ambiguous nucleotide codon that doesn't matter in the translation
                    possible_nt1 = AmbiguousNTCodes(codon[0])
                    possible_nt2 = AmbiguousNTCodes(codon[1])
                    possible_nt3 = AmbiguousNTCodes(codon[2])
                    possible_codons = []
                    for nt1 in possible_nt1:
                        for nt2 in possible_nt2:
                            for nt3 in possible_nt3:
                                possible_codons.append("%s%s%s" % (nt1, nt2, nt3))
                    try:
                        aa = genetic_code[possible_codons[0]]
                    except KeyError:
                        raise KeyError("Cannot translate codon %s in %s" % (codon, head))
                    for possible_codon in possible_codons:
                        if genetic_code[possible_codon] != aa:
                            if readthrough_n:
                                aa = 'X'
                            else:
                                raise ValueError("Cannot translate codon %s" % codon)
            if aa == 'STOP' and i == prot_length - 1:
                aa = ''
            elif aa == 'STOP':
                if readthrough_stop:
                    aa = 'X'
                else:
                    raise ValueError("Premature stop codon")
            prot.append(aa)
        translated_headers_sequences.append((head, ''.join(prot)))
    return translated_headers_sequences


def AmbiguousNTCodes(nt):
    """Returns all possible nucleotides corresponding to an ambiguous code.

    This method takes as input a single nucleotide character, which is
    assumed to represent a nucleotide as one of the accepted
    codes for an ambiguous character.  Returns a list giving
    all possible codes for which a nucleotide might stand.  Raises
    an exception if `nt` is not a valid nucleotide code.

    EXAMPLES:

    >>> AmbiguousNTCodes('N')
    ['A', 'T', 'G', 'C']

    >>> AmbiguousNTCodes('R')
    ['A', 'G']

    >>> AmbiguousNTCodes('A')
    ['A']

    >>> AmbiguousNTCodes('-')
    ['-']

    >>> AmbiguousNTCodes('F')
    Traceback (most recent call last):
       ...
    ValueError: Invalid nt code of "F"
    """
    if nt in ['A', 'T', 'G', 'C', '-']:
        return [nt]
    elif nt == 'R':
        return ['A', 'G']
    elif nt == 'Y':
        return ['T', 'C']
    elif nt == 'K':
        return ['G', 'T']
    elif nt == 'M':
        return ['A', 'C']
    elif nt == 'S':
        return ['G', 'C']
    elif nt == 'W':
        return ['A', 'T']
    elif nt == 'B':
        return ['C', 'G', 'T']
    elif nt == 'D':
        return ['A', 'G', 'T']
    elif nt == 'H':
        return ['A', 'C', 'T']
    elif nt == 'V':
        return ['A', 'C', 'G']
    elif nt == 'N':
        return ['A', 'T', 'G', 'C']
    else: 
        raise ValueError('Invalid nt code of "%s"' % nt)



def PurgeAmbiguousDNA(headers_sequences):
    """Removes all sequences with ambiguous positions from nucleotide sequences.

    This function takes a single calling argument `headers_sequences`, which
    is a list of tuples `(header, seq)` as would be returned by `Read`.
    These sequences should specify nucleotide sequences.  It returns
    a new list which is a copy of 'headers_sequences', except that all
    sequences that contain ambiguous nucleotide entries (i.e. characters
    that are not 'A', 'T', 'C', 'G', 'a', 't', 'c', or 'g') have been
    removed.
    """
    valid_bases = {'A':1, 'T':1, 'C':1, 'G':1, 'a':1, 't':1, 'c':1, 'g':1}
    assert isinstance(headers_sequences, list)
    purged_headers_sequences = []
    for (head, seq) in headers_sequences:
        for nt in seq:
            if nt not in valid_bases:
                break
        else:
            purged_headers_sequences.append((head, seq))
    return purged_headers_sequences


def GetEntries(namelist, fastafile, allow_substring=False):
    """Gets selected entries from a (potentially very large) FASTA file.

    This method is designed to extract sequences from a FASTA file.  It will work
    even if the FASTA file is very large, since it avoids reading the entire
    file into memory at once. 

    `namelist` : specifies names of the sequences that we want to extract from
    the FASTA file.  The name of a sequence is the string immediately following
    the ">" in the FASTA file header for a sequence, terminated by a space character
    (space, tab, or return).  So for example, the header::

        >E_coli_thioredoxin: the thioredoxin protein from E. coli

    would correspond to a name of "E_coli_thioredoxin".  `namelist` specifies
    a list of these names.

    `fastafile` : name of a FASTA file that contains the sequences we are searching
    for.  For this method to be guaranteed to work properly, each sequence in the FASTA
    file must contain a unique name, where a "name" is as defined above.  Note that this
    uniqueness of names is not rigorously checked for, so if there are not unique names,
    the function may raise an exception, or it may continue along and give no hint of
    the problem.

    `allow_substring` : Boolean switch that specifies that the name given in
    `namelist` need only be a substring of the first entry in the `fastafile` header.

    The function expects to find exactly one entry in `fastafile` for each name in
    `namelist`.  If it does not, it will raise an exception.  The returned variable
    is a list composed of 2-tuples.  Element `i` of this list corresponds to the name
    given by `namelist[i]`.  Each 2-tuple has the form `(header, sequence)` where
    `header` is the full FASTA header, but with the leading ">" character and any trailing
    linebreaks/spaces removed.  `sequence` is a string giving the sequence, again with the
    trailing linebreak removed.
    """
    namedict = {}
    for name in namelist:
        namedict[name] = None
    f = open(fastafile)
    line = f.readline()
    header = None
    seq = []
    while line:
        if line[0] == '>':
            if (not header) and (not seq):
                pass # first sequence in file
            elif header and not seq:
                raise ValueError, "Empty sequence for %s" % header
            elif seq and not header:
                raise ValueError, "File does not begin with header."
            else:
                name = header.split()[0]
                if name in namedict:
                    if namedict[name]:
                        raise ValueError, "Duplicate entries for name %s" % name
                    else:
                        namedict[name] = (header, ''.join(seq))
                elif allow_substring:
                    subnames = [iname for iname in namedict.iterkeys() if iname in name]
                    if len(subnames) > 1:
                        raise ValueError("Multiple subname matches for %s." % name)
                    elif subnames:
                        if namedict[subnames[0]]:
                            raise ValueError("Duplicate subname entries for name %s" % name)
                        else:
                            namedict[subnames[0]] = (header, ''.join(seq))
            header = line.strip()[1 : ]
            seq = []
        else:
            seq.append(line.strip())
        line = f.readline()
    f.close()
    if (not header) and (not seq):
        pass # first sequence in file
    elif header and not seq:
        raise ValueError, "Empty sequence for %s" % header
    elif seq and not header:
        raise ValueError, "File does not begin with header."
    else:
        name = header.split()[0]
        if name in namedict:
            if namedict[name]:
                raise ValueError, "Duplicate entries for name %s" % name
            else:
                namedict[name] = (header, ''.join(seq))
    header_seq_list = []
    for name in namelist:
        if not namedict[name]:
            raise ValueError, "No entry for %s" % name
        header_seq_list.append(namedict[name])
    return header_seq_list


def AAThreeToOne(aathree):
    """Converts a three letter amino acid code into a one letter code.

    The single input argument `aathree` is the three letter amino acid code.
    It can be upper or lower case.

    This function returns, in upper case, the one letter amino acid code.
    It raises an exception if `aathree` is not a valid one letter code.
    'Xaa' is converted to 'X'.

    EXAMPLES:

    >>> AAThreeToOne('Ala')
    'A'

    >>> AAThreeToOne('cys')
    'C'

    >>> AAThreeToOne('Xaa')
    'X'

    >>> AAThreeToOne('hi')
    Traceback (most recent call last):
       ...
    ValueError: Invalid amino acid code of hi.
    """
    aa_mapping = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY', 'H':'HIS',
                'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN',
                'R':'ARG', 'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR', 'X':'XAA'}
    aa_reverse_mapping = {}
    for (one, three) in aa_mapping.iteritems():
        aa_reverse_mapping[three] = one
    try:
        return aa_reverse_mapping[aathree.upper()]
    except KeyError:
        raise ValueError("Invalid amino acid code of %s." % aathree)


def ReverseComplement(seqs, use_csequtils=True):
    """Converts nucleotide sequences to their reverse complements.

    The input argument `seqs` specifies the sequence(s) to reverse
    complement. It can be:

    * a string giving a single sequence

    * a list of strings giving sequences

    * a list of 2-tuples of the format `(head, seq)` where `head` is the
      header and `seq` is the sequence (both strings).

    The sequences should be nucleotides A, T, C, or G
    (or their lowercase equivalents). N nucleotides are
    are reverse-complemented to N nucleotides. Other nucleotide
    ambiguous nucleotide codes are currently not accepted.

    This function also takes an optional argument `use_csequtils`, which
    is `True` by default. This specifies that we use the fast C extension
    available in `csequtils` if it is available. Set this to `False` if you
    want to use only the pure Python code.

    The returned variable is a copy of `seqs` in which any headers
    are unchanged but the sequences are converted to reverse complements.

    EXAMPLES:

    >>> ReverseComplement([('seq1', 'ATGCAA'), ('seq2', 'atgGCA')])
    [('seq1', 'TTGCAT'), ('seq2', 'TGCcat')]

    >>> ReverseComplement([('seq1', 'ATGCAA'), ('seq2', 'ntgGNA')])
    [('seq1', 'TTGCAT'), ('seq2', 'TNCcan')]

    >>> ReverseComplement([('seq1', 'ATGHAA')])
    Traceback (most recent call last):
        ...
    ValueError: Invalid nucleotide code.

    >>> ReverseComplement(['ATGCA', 'GCTga'])
    ['TGCAT', 'tcAGC']

    >>> ReverseComplement('ATGACAGTA')
    'TACTGTCAT'
    
    >>> ReverseComplement('ATGAcaGTA')
    'TACtgTCAT'
    """
    if isinstance(seqs, str):
        if _csequtilsimported and use_csequtils:
            return mapmuts.csequtils.ReverseComplement(seqs)
        mapping = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'n':'n', 'N':'N'}
        rc = [mapping[nt] for nt in seqs]
        rc.reverse()
        return ''.join(rc)
    elif seqs == []:
        return [] # no sequences
    elif isinstance(seqs, list) and isinstance(seqs[0], str):
        if _csequtilsimported and use_csequtils:
            return [mapmuts.csequtils.ReverseComplement(s) for s in seqs]
        mapping = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'n':'n', 'N':'N'}
        for seq in seqs:
            try:
                rc = [mapping[nt] for nt in seq]
            except KeyError:
                raise ValueError("Invalid nucleotide code.")
            rc.reverse()
            reversecomplements.append(''.join(rc))
    elif isinstance(seqs, list) and isinstance(seqs[0], tuple) and\
            len(seqs[0]) == 2:
        if _csequtilsimported and use_csequtils:
            return [(h, mapmuts.csequtils.ReverseComplement(s)) for (h, s) in seqs]
        mapping = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'n':'n', 'N':'N'}
        reversecomplements = []
        for (head, seq) in seqs:
            try:
                rc = [mapping[nt] for nt in seq]
            except KeyError:
                raise ValueError("Invalid nucleotide code.")
            rc.reverse()
            reversecomplements.append((head, ''.join(rc)))
    else:
        raise ValueError("Invalid value for seqs: %s"
                         % str(seqs))
    return reversecomplements


def FindMotifs(seq, motif):
    """Finds occurrences of a specific motif in a nucleotide sequence.

    `seq` is a string giving a nucleotide sequence.

    `motif` is a string giving the motif that we are looking for.  It should
    be a string of valid nucleotide characters:

        * A   Adenine

        * G   Guanine

        * C   Cytosine

        * T   Thymine

        * U   Uracil

        * R   Purine (A or G)

        * Y   Pyrimidine (C or T)
        
        * N   Any nucleotide

        * W   Weak (A or T)

        * S   Strong (G or C)

        * M   Amino (A or C)

        * K   Keto (G or T)

        * B   Not A (G or C or T)

        * H   Not G (A or C or T)

        * D   Not C (A or G or T)

        * V   Not T (A or G or C)

    The returned variable is a list `motif_indices` of the indices that
    each occurrence of `motif` in `seq` begins with.  For example,
    if there is a motif beginning at `seq[7]`, then 7 will be
    present in `motif_indices`.  So the number of occurrences of
    the motif will be equal to `len(motif_indices)`.

    This function is not case sensitive, so nucleotides can be either upper
    or lower case.  In addition, T (thymine) and U (uracil) nucleotides
    are treated identically, so the function can handle either DNA
    or RNA sequences.

    EXAMPLES:

    >>> FindMotifs('ATCGAA', 'WCGW')
    [1]
    """
    assert isinstance(seq, str)
    assert isinstance(motif, str)
    seq = seq.upper()
    motif = motif.upper()
    seq = seq.replace('U', 'T')
    motif = motif.replace('U', 'T')
    nt_mapping = {  # maps nucleotide codes to regular expression code
        'A' : 'A',
        'G' : 'G',
        'C' : 'C',
        'T' : 'T',
        'R' : '[A,G]',
        'Y' : '[C,T]',
        'N' : '[A,T,C,G]',
        'W' : '[A,T]',
        'S' : '[G,C]',
        'M' : '[A,C]',
        'K' : '[G,T]',
        'B' : '[G,C,T]',
        'H' : '[A,C,T]',
        'D' : '[A,G,T]',
        'V' : '[A,G,C]',
    }
    try:
        m = re.compile(''.join([nt_mapping[nt] for nt in motif]))
    except KeyError:
        raise ValueError("motif contains an invalid character: %s" % motif)
    motif_indices = [match.start() for match in m.finditer(seq)]
    return motif_indices


def CondenseSeqs(seqs, maxdiffs, exclude_positions):
    """Removes nearly identical protein sequences.

    `seqs` : list of sequences as `(head, seq)` 2-tuples. The sequences
    are assumed to be aligned such that they are all of the same length..

    `maxdiffs` : specifies the maximum number of differences that a sequence
    can have from another sequence in order to be removed.

    `exclude_positions` : list of integers specifying positions at which
    sequence can NOT differ and still be removed.  These integers
    are for numbering the sequences as 1, 2, ...

    The method proceeds as follows:

        1) For the first sequence iterates through the rest of the sequences,
           and removes any that differ at <= `maxdiffs` sites, and do not differ
           at the sites specified by `exclude_position`.

        2) Repeats this process for each of the remaining sequences. The
           returned variable is a copy of `seqs` with the removed sequences gone.

    EXAMPLES:

    >>> seqs = [('s1', 'ATGC'), ('s2', 'ATGA'), ('s3', 'ATC-'), ('s4', 'ATGC')]
    >>> CondenseSeqs(seqs, 0, [])
    [('s1', 'ATGC'), ('s2', 'ATGA'), ('s3', 'ATC-')]
    >>> CondenseSeqs(seqs, 1, [])
    [('s1', 'ATGC'), ('s3', 'ATC-')]
    >>> CondenseSeqs(seqs, 1, [4])
    [('s1', 'ATGC'), ('s2', 'ATGA'), ('s3', 'ATC-')]
    """
    i = 0
    while i < len(seqs):
        newseqs = list(seqs[ : i + 1])
        iseq = seqs[i]
        for jseq in seqs[i + 1 : ]:
            assert len(iseq[1]) == len(jseq[1])
            differ_at_excluded = False
            for position in exclude_positions:
                if iseq[1][position - 1] != jseq[1][position - 1]:
                    differ_at_excluded = True
                    break
            if differ_at_excluded:
                newseqs.append(jseq)
            else:
                ndiffs = len([k for k in range(len(iseq[1])) if iseq[1][k] != jseq[1][k]])
                if ndiffs > maxdiffs:
                    newseqs.append(jseq)
        seqs = newseqs
        i += 1
    return seqs



def KyteDoolittle(aa):
    """Returns the hydrophobicity of an amino acid based on the Kyte Doolittle scale.

    'aa' is the one letter code for an amino acid.

    The returned value is a number giving the hydrophobicity, as defined by Kyte
    and Doolittle in::

            J. Kyte & R. F. Doolittle: 
            "A simple method for displaying the hydropathic character of a protein." 
            J Mol Biol, 157, 105-132

    More positive values indicate higher hydrophobicity, while more negative values
    indicate lower hydrophobicity.

    >>> print round(KyteDoolittle('F'), 2)
    2.8

    >>> print round(KyteDoolittle('N'), 2)
    -3.5

    >>> print round(KyteDoolittle('X'), 2)
    Traceback (most recent call last):
        ...
    ValueError: Invalid amino acid code of X.
    """
    aa = aa.upper()
    d = {'A':1.8, 'C':2.5, 'D':-3.5, 'E':-3.5, 'F':2.8, 'G':-0.4,
         'H':-3.2, 'I':4.5, 'K':-3.9, 'L':3.8, 'M':1.9, 'N':-3.5,
         'P':-1.6, 'Q':-3.5, 'R':-4.5, 'S':-0.8, 'T':-0.7, 
         'V':4.2, 'W':-0.9, 'Y':-1.3}
    try:
        return d[aa]
    except KeyError:
        raise ValueError("Invalid amino acid code of %s." % aa)


# Test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()
