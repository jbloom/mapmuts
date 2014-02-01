"""Module for simulating data for ``mapmuts`` package.

This module contains function to simulate data that can be analyzed
by the ``mapmuts`` package. This simulated data can be useful for
testing the package's performance.

Use the `Seed` function if you want to ensure repeatable output.

Written by Jesse Bloom.


List of functions
--------------------
`Seed` : seeds random number generators used by module.

`MutateNT` : mutates a nucleotide.

`SimulateIlluminaPE` : simulates uniformly fragmented Illumina paired-end library.

`SimulateCounts` : simulates library counts for inference of enrichment ratio.

`SimulatePreferences` : simulates counts for inference of equilibrium preferences.


Documentation for functions
-----------------------------
Documentation for individual functions is provided in their definitions below.
"""


import random
import mapmuts.sequtils
import mapmuts.bayesian



def Seed(seed):
    """Seeds all random number generators used by this module.

    This function seeds all random number generators used by the
    functions in this module. This creates reproducible output if the
    seeding is done before calling other module functions. 

    `seed` : the integer seed used to seed the random number generators.

    Currently the generators that are seeded are:

        * random.seed
    """
    random.seed(seed)


def MutateNT(nt, mutrate):
    """Mutates a nucleotide with probability `mutrate`.

    `nt` : a DNA nucleotide character, either upper or lower case (i.e.
    'A', 'C', 'T', 'G', 'a', 'c', 't', 'g')

    `mutrate` : the probability with which the nucleotide is randomly
    mutated to some other nucleotide, so `0 <= mutrate <= 1`.

    Returns a character giving the original nucleotide (with probability
    `1 - mutrate`) or one of the other three nucleotide (with probability
    `mutrate`). Nucleotides are mutated to any of the other three
    nucleotides with equal probability. If the original nucleotide
    is upper-case, the mutant will be too. If the original
    nucleotide is lower-case, the mutant will be too.
    """
    assert 0 <= mutrate <= 1
    if random.random() < mutrate:
        return random.choice({'A':['T', 'C', 'G'], 'T':['A', 'C', 'G'],
                              'C':['T', 'A', 'G'], 'G':['A', 'T', 'C'],
                              'a':['t', 'c', 'g'], 't':['a', 'c', 'g'],
                              'c':['t', 'a', 'g'], 'g':['a', 't', 'c'],
                              }[nt])
    else:
        return nt


def SimulateIlluminaPE(gene, r1_trim3, r2_trim3, readlength,
                        insertlength, libsize, mutrate):
    """Simulates a uniformly fragmented Illumina paired-end library.

    This function simulates creating an Illumina paired-end read
    library. A target nucleotide sequence (given by `gene`) is randomly
    fragmented into small pieces as specified by `insertlength`. The
    inserts are then randomly reverse complemented with a probability of
    0.5. The adaptor sequence specified by `r1_trim3` is then appended to
    the 3' end of the fragment to generate the read 1 (R1) read. The
    read fragment is reverse complemented and the adaptor sequence
    specified by `r2_trim3` is then appended to the 3' of this fragment to
    generate the corresponding read 2 (R2 read). Each fragment-adaptor
    sequence is then subjected to random mutation as specified by
    `mutrate`. Finally, all fragment-adaptor sequences are truncated to
    the length specified by `readlength`. Overall, this creates data
    similar to that which would be expected from an Illumina PE read.
    
    If `readlength` exceeds the length of the gene fragment plus the
    adaptor sequence, random nucleotides are appended at the 3' end
    until `readlength` is reached.

    CALLING VARIABLES:

    `gene` : a string of DNA nucleotides giving the sequence of the gene that
    is randomly fragmented.

    `r1_trim3` : a string of DNA nucleotides giving the adaptor sequence that
    is appended to the 3' end of R1 fragments.

    `r2_trim3` : a string of DNA nucleotides giving the adaptor sequence that
    is appended to the 3' end of R2 fragments.

    `readlength` : the integer length of the Illumina reads.

    `insertlength` : specifies the length of the inserts (fragments of
    `gene`) that are included in each read between the adaptors. This
    argument should be a 3-tuple. It is of the form `('uniform',
    minlength, maxlength)` where `minlength` and `maxlength` are
    integers. The inserts lengths are then drawn from uniformly from
    the distribution of lengths >= `minlength` and <= `maxlength`.

    `libsize` : the number of read pairs to create.

    `mutrate` : the probabily that each base in the total read
    (insert and adaptor) is randomly mutated to another nucleotide.
    This is designed to simulate sequencing errors. Each nucleotide
    is randomly mutated to any other nucleotide. The mutations are
    applied separately to each read. So for example, if there are
    two 50 nt reads (`readlength == 50`) and `mutrate = 0.01`, then
    there will be an average of one mutation total between the two
    reads (100 nt times 0.01 mutation rate).

    RETURN VARIABLE:

    The returned variable is a 2-tuples of lists, `(r1_reads,
    r2_reads)`. `r1_reads` and `r2_reads` are lists of the same length (with
    ``len(r1_reads) == len(r2_reads) == libsize``). Element `i` of `r1_reads`
    is a string giving that R1 read, and element `i` of `r2_reads` is the
    corresponding string giving that R2 read. In the returned reads, the
    portion of the read that comes from the target gene is upper-case,
    while the portion that comes from the adaptor is lower-case,
    regardless of the case of the input variables for these sequences.
    """
    gene_rc = mapmuts.sequtils.ReverseComplement([gene])[0]
    n = len(gene)
    r1_trim3 = r1_trim3.lower()
    r2_trim3 = r2_trim3.lower()
    if isinstance(insertlength, tuple) and len(insertlength) == 3 and\
            insertlength[0] == 'uniform':
        (minlength, maxlength) = (insertlength[1], insertlength[2])
        assert 0 <= minlength <= maxlength
    else:
        raise ValueError("Invalid value of insertlength: %r" %
                         insertlength)
    assert 0 <= mutrate <= 1
    r1_reads = []
    r2_reads = []
    for i in range(int(libsize)):
        l = random.randint(minlength, maxlength) # length of insert
        fragstart = random.randint(0, n - l) # index start
        if random.random() < 0.5:
            r1 = gene[fragstart : fragstart + l]
        else:
            r1 = gene_rc[fragstart : fragstart + l]
        r1 = r1.upper()
        r2 = mapmuts.sequtils.ReverseComplement([r1])[0]
        r1 = ("%s%s" % (r1, r1_trim3))[ : readlength] 
        r2 = ("%s%s" % (r2, r2_trim3))[ : readlength] 
        r1 = ''.join([MutateNT(nt, mutrate) for nt in r1])
        r2 = ''.join([MutateNT(nt, mutrate) for nt in r2])
        r1_reads.append(r1)
        r2_reads.append(r2)
    assert len(r1_reads) == len(r2_reads)
    return (r1_reads, r2_reads)



def SimulatePreferences(pi, mu, epsilon, rho, wtcodon, nlibs, libsize, perturb, seed=1):
    """Simulates library counts to test inference of equilibrium preferences.

    This function simulates data for analysis by 
    *mapmuts.bayesian.InferPreferencesMCMC*. The variables and what is
    being simulated may be clearer after reading the documentation for
    that function.

    This function requires `` numpy`` and ``pymc``, so trying to run this
    function if *mapmuts.bayesian.PymcAvailable() == False* will raise
    an exception.

    CALLING VARIABLES:

    * *pi* : a dictionary keyed by all 21 amino acids in 
      *mapmuts.sequtils.AminoAcids(includestop=True)*, with the
      values giving the preferences for each amino acid (must sum to one).

    * *mu* : the mutagenesis rate in the *mutDNA* library.

    * *epsilon* : a dictionary keyed by 1, 2, and 3 with *epsilon[i]*
      giving the sequencing error rate to codons with *i* differences from *wtcodon*.

    * *rho* : a dictionary keyed by 1, 2, and 3 with *rho[i]* giving
      the reverse-transcription rate to codons with *i* differences from *wtcodon*.

    * *wtcodon* : string giving the wild-type codon, such as *GCA*.

    * *nlibs* : the number of libraries being analyzed.

    * *libsize* : the number of sequences in each library.

    * *perturb* : a number specifying how much we perturb the values of *libsize*
      for each sample, and how much we perturb *epsilon*, *rho*, and *mu* for
      each mutant codon, away from the set values. For each simulation of one of
      these, we set the simulation value to the calling value plus a number
      drawn uniformly from the calling value divided *perturb* to the calling
      value plus *perturb*. So *perturb* must be >= 1, and a reasonable value 
      might be 2.0. This tests if the inference still works even if the priors
      aren't quite accurate.

    * *seed* random number seed, 1 by default. Used to seed ``numpy`` seed,
      which is not seeded by *mapmuts.simulate.Seed*.

    OPERATION AND RESULT:

    This function simulates the indicated number and sizes of libraries using
    the specified calling values, but perturbing the appropriate ones as
    specified by *perturb*. The returned variable is *library_stats* which
    has all of the values that you would need to pass this variable
    to *mapmuts.bayesian.InferPreferencesMCMC* for inference. The priors
    for *mu*, *epsilon*, and *rho* are set to the calling values.

    The idea is that you can then apply the inference function to the results
    and test how well you actually infer the true enrichment ratio.
    """
    if not mapmuts.bayesian.PymcAvailable():
        raise ImportError("Cannot import numpy")
    import numpy
    numpy.random.seed(seed)
    codons = mapmuts.sequtils.Codons()
    assert numpy.allclose(1.0, sum(pi.values())), "pi does not sum to one"
    c_pi = []
    for codon in codons:
        aa = mapmuts.sequtils.Translate([('head', codon)])[0][1]
        if not aa:
            aa = '*'
        c_pi.append(pi[aa])
    c_pi = numpy.array(c_pi)
    perturb = float(perturb)
    assert perturb >= 1, "perturb must be >= 1"
    assert wtcodon in codons, "Not a valid wtcodon of %s" % wtcodon
    iwtcodon = codons.index(wtcodon)
    ndiffs = [len([i for i in range(len(codon)) if codon[i] != wtcodon[i]]) for codon in codons]
    rho[0] = 0.0
    epsilon[0] = 0.0
    library_stats = []
    for ilib in range(nlibs):
        libstats = {'mu_prior':mu, 'epsilon_prior':epsilon, 'rho_prior':rho, 'wtcodon':wtcodon}
        muvec = numpy.array([random.uniform(mu / perturb, mu * perturb) for codon in codons])
        muvec[iwtcodon] = 0
        muvec[iwtcodon] = 1.0 - numpy.sum(muvec)
        assert muvec[iwtcodon] > 0, "mu entry for wildtype not > 0, mu too large?"
        rhovec = numpy.array([random.uniform(rho[ndiff] / perturb, rho[ndiff] * perturb) for ndiff in ndiffs])
        rhovec[iwtcodon] = 1.0 - numpy.sum(rhovec)
        assert rhovec[iwtcodon] > 0, "rho entry for wildtype not > 0, rho too large?"
        epsilonvec = numpy.array([random.uniform(epsilon[ndiff] / perturb, epsilon[ndiff] * perturb) for ndiff in ndiffs])
        epsilonvec[iwtcodon] = 1.0 - numpy.sum(epsilonvec)
        assert epsilonvec[iwtcodon] > 0, "epsilon entry for wildtype not > 0, epsilon too large?"
        delta = numpy.zeros(len(codons))
        delta[iwtcodon] = 1.0
        libstats['nrdna_counts'] = dict(zip(codons, numpy.random.multinomial(random.uniform(libsize / perturb, libsize * perturb), epsilonvec)))
        libstats['nrrna_counts'] = dict(zip(codons, numpy.random.multinomial(random.uniform(libsize / perturb, libsize * perturb), epsilonvec + rhovec - delta)))
        libstats['nrmutdna_counts'] = dict(zip(codons, numpy.random.multinomial(random.uniform(libsize / perturb, libsize * perturb), epsilonvec + muvec - delta)))
        libstats['nrmutvirus_counts'] = dict(zip(codons, numpy.random.multinomial(random.uniform(libsize / perturb, libsize * perturb), epsilonvec + rhovec + c_pi * muvec / numpy.dot(c_pi, muvec) - 2.0 * delta)))
        library_stats.append(libstats)
    return library_stats



def SimulateCounts(phi, mu, epsilon, rho, wtcodon, codons, libsize, nlibs, perturb):
    """Simulates library counts to test inference of enrichment ratio.

    This function is designed to simulate data for analysis by
    *mapmuts.bayesian.InferEnrichmentMCMC*. The variables and what is
    being simulated may be clearer after reading the documentation for that
    function.

    This function requires ``numpy``, which is checked for by 
    *mapmuts.bayesian.PymcAvailable()*. Trying to run this function
    if that returns *False* will raise an exception, which is fine -- 
    because this function is designed to test the Bayesian inference,
    you wouldn't want to run it without being able to do the inference anyway!

    CALLING VARIABLES:

    * *phi* : the true enrichment ratio.

    * *mu* : the mutagenesis rate in the *mutDNA* library.

    * *epsilon* : a dictionary keyed by 1, 2, and 3 with *epsilon[i]*
      giving the sequencing error rate to codons with *i* differences from *wtcodon*.

    * *rho* : a dictionary keyed by 1, 2, and 3 with *rho[i]* giving
      the reverse-transcription rate to codons with *i* differences from *wtcodon*.

    * *wtcodon* : string giving the wild-type codon, such as *GCA*.

    * *codons* : a list of one or more mutant codons, such as
      *[GGC, GGG, GGA, GGT]*

    * *libsize* : the number of sequences in each library.

    * *nlibs* : the number of libraries being analyzed.

    * *perturb* : a number specifying how much we perturb the values of *libsize*
      for each sample, and how much we perturb *epsilon*, *rho*, and *mu* for
      each mutant codon, away from the set values. For each simulation of one of
      these, we set the simulation value to the calling value plus a number
      drawn uniformly from the calling value divided *perturb* to the calling
      value plus *perturb*. So *perturb* must be >= 1, and a reasonable value 
      might be 2.0. This tests if the inference still works even if the priors
      aren't quite accurate.

    OPERATION AND RESULT:

    This function simulates the indicated number and sizes of libraries using
    the specified calling values, but perturbing the appropriate ones as
    specified by *perturb*. The returned variable is *library_stats* which
    has all of the values that you would need to pass this variable
    to *mapmuts.bayesian.InferEnrichmentMCMC* for inference. The priors
    for *mu*, *epsilon*, and *rho* are set to the calling values.

    The idea is that you can then apply the inference function to the results
    and test how well you actually infer the true enrichment ratio.
    """
    if not mapmuts.bayesian.PymcAvailable():
        raise ImportError("Cannot import numpy")
    import numpy
    perturb = float(perturb)
    assert perturb >= 1, "perturb must be >= 1"
    libtypes = ['dna', 'rna', 'mutdna', 'mutvirus']
    library_stats = []
    for ilib in range(nlibs):
        libstats = {'mu_prior':mu, 'epsilon_prior':epsilon, 'rho_prior':rho}
        for libtype in libtypes:
            libstats["Nr%s" % libtype] = random.uniform(libsize / perturb, libsize * perturb)
            libstats['nr%s_list' % libtype] = []
        for codon in codons:
            ndiffs = len([i for i in range(len(codon)) if codon[i] != wtcodon[i]])
            assert 1 <= ndiffs <= 3
            imu = random.uniform(mu / perturb, mu * perturb)
            iepsilon = random.uniform(epsilon[ndiffs] / perturb, epsilon[ndiffs] * perturb)
            irho = random.uniform(rho[ndiffs] / perturb, rho[ndiffs] * perturb)
            for libtype in libtypes:
                ilibsize = libstats['Nr%s' % libtype]
                if libtype == 'dna':
                    nexpected = iepsilon * ilibsize
                elif libtype == 'rna':
                    nexpected = (iepsilon + irho) * ilibsize
                elif libtype == 'mutdna':
                    nexpected = (iepsilon + imu) * ilibsize
                elif libtype == 'mutvirus':
                    nexpected = (iepsilon + irho + phi * imu) * ilibsize
                else:
                    raise ValueError("Invalid libtype")
                n = numpy.random.poisson(nexpected)
                libstats['nr%s_list' % libtype].append(n)
        library_stats.append(libstats)
    return library_stats


if __name__ == '__main__':
    import doctest
    doctest.testmod()
