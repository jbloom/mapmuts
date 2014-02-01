"""Module for performing Bayesian inference for the ``mapmuts`` package.

This module uses the ``pymc`` package, and has been tested with version 2.2
and 2.3. It works with both, but will be faster with version 2.3.

The ``pymc`` package also requires ``numpy``. Many functions 
in this module will therefore fail if ``pymc`` and ``numpy`` are not available
for import. Before running any function in this module, you can use the
*PymcAvailable* function to determine if they are available.

In addition, some functions perform better if ``scipy`` is available, although
this is not strictly required.

Written by Jesse Bloom, 2013.


List of functions
--------------------
* *PymcAvailable* : returns *True* if and only if ``numpy`` / ``pymc`` are available.

* *ScipyAvailable* : returns *True* if and only if ``scipy`` is available.

* *Seed* : seeds the random number generator.

* *InferEnrichmentMCMC* : infers enrichment ratio from MCMC over posterior.

* *InferEnrichmentMCMC_2* : faster (usually) version of *InferEnrichmentMCMC_2*.

* *InferPreferencesMCMC* : infers equilibrium amino-acid preferences using MCMC.

* *EquilibriumFracs* : calculates equilibrium preferences for amino acids.

* *CredibleInterval* : computes median-centered credible interval.

* *SiteEntropy* : calculates site entropy from equilibrium fractions.

* *ShannonJensenDivergence* : calculates Shannon-Jensen Divergence.

* *PreferencesRemoveStop* : removes stop codon as possible amino acid from preferences.


Details of functions
----------------------
Documentation for functions is provided in their individual docstrings below.

"""


import os
import math
import sys
import copy
import warnings
import math
import random
import mapmuts.plot
import mapmuts.sequtils
# global variable _pmycavailable indicates if pymc and numpy are available
try:
    import numpy
    import pymc
    _pymcavailable = True
except ImportError:
    _pymcavailable = False

try:
    import scipy
    _scipyavailable = True
except ImportError:
    _scipyavailable = False


def PymcAvailable():
    """Tests if ``numpy`` and ``pymc`` are available.

    Returns *True* if these packages can be imported, and *False* otherwise."""
    return _pymcavailable


def ScipyAvailable():
    """Tests if ``scipy`` is available.

    Returns *True* if ``scipy`` can be imported, and *False* otherwise."""
    return _scipyavailable


def PreferencesRemoveStop(preferences):
    """Removes stop codon as a possible amino acid from site preferences.

    You should use this function if you have a set of amino acid
    *preferences* that include a stop codon (denoted by a * character) as a possibility,
    and you want to remove the possibility of that stop codon and
    renormalize all of the other preferences so that they sum to one.

    *preferences* is a variable that stores the preferences of sites
    for specific amino acids. Typically, these would be the values 
    written by ``mapmuts_inferpreferences.py`` or ``mapmuts_inferenrichment.py``
    as the ``*_equilibriumpreferences.txt`` or ``*_equilibriumfreqs.txt`` files,
    and then read by *mapmuts.io.ReadEntropyAndEquilFreqs*. So
    *preferences* is a dictionary keyed by site number *r*, and for each
    site *r*, *preferences[r]* is in turn a dictionary. This dictionary
    has the string keys 'WT_AA', 'SITE_ENTROPY', and then 'PI_A', 'PI_C',
    'PI_D', etc. The 'PI_A' entry gives the preference for amino acid
    *A*, etc. There may be 20 entries for the preferences for all 20 amino
    acids, or there may be 21 entries for all 20 amino acids plus for a 
    stop codon (this would be indicated by 'PI_*').

    This script goes through all sites *r* that key *preferences*. For a site,
    if there is **not** a stop codon key ('PI_*') in *preferences[r]*, then
    nothing is done. But if there is a stop codon entry, then the key
    'PI_*' is removed from the dictionary and all of the other preferences
    are renormalized so that they stay proportional to their old values
    but now sum to one.

    These changes are made to *preferences* in place, and there is no
    return value. Because *preferences* is a dictionary, it is mutable,
    so that is why the changes can be made in place.
    """
    aas = mapmuts.sequtils.AminoAcids(includestop=False)
    for r in preferences.keys():
        if 'PI_*' in preferences[r]:
            # has stop codon
            del preferences[r]['PI_*']
            rsum = float(sum([preferences[r]['PI_%s' % aa] for aa in aas]))
            if not (0 < rsum <= 1):
                raise ValueError("Non-stop codon PI sums not > 0 and <= 1")
            for aa in aas:
                preferences[r]['PI_%s' % aa] /= rsum
            pi_sum = sum([preferences[r]['PI_%s' % aa] for aa in aas]) 
            assert abs(pi_sum - 1) < 1e-7, 'PI values do not sum to close to one, sum is %g' % pi_sum


def ShannonJensenDivergence(p1, p2):
    """Returns the Shannon-Jensen divergence between *p1* and *p2*.

    Requires ``numpy``, which is available if *PymcAvailable()* evaluates
    to *True*. Otherwise you will get an error.

    *p1* and *p2* are two *numpy.ndarray* objects with elements that give
    the probability distributions for which we are computing the divergence.
    Specifically, *p1[i]* and *p2[i]* indicate the probability for state *i*
    in the two distributions. The constraints are:
    
        * *p1* and *p2* should have dimension of one
          (*1 == len(p1.shape) == len(p2.shape)*).
          
        * *p1* and *p2* should have the same
          nonzero length (*len(p1) == len(p2) > 0*)
          
        * *p1* and *p2* should both have entries that sum to one (since
          they are probability distributions), so
          *True == numpy.allclose(sum(p1), 1) == numpy.allclose(sum(p2), 1)*

        * *p1* and *p2* should have all entries >= 0, so 
          *True == numpy.all(p1 >= 0) == numpy.all(p2 >= 0)*.

    This function returns the Shannon-Jensen divergence between the 
    probability distributions specified by *p1* and *p2*. The logarithms
    are taken to the base 2, meaning that the returned divergence
    will always be between 0 and 1.
    """
    if not _pymcavailable:
        raise ValueError("numpy may not be available.")
    assert 1 == len(p1.shape) == len(p2.shape), "p1 and p2 not both numpy.ndarrays of dimension 1"
    assert len(p1) == len(p2) > 0, "p1 and p2 not both arrays of same nonzero length"
    assert numpy.allclose(sum(p1), 1), "p1 does not have entries summing to one"
    assert numpy.allclose(sum(p2), 1), "p2 does not have entries summing to one"
    assert numpy.all(p1 >= 0), "p1 does not have all entries >= 0"
    assert numpy.all(p2 >= 0), "p2 does not have all entries >= 0"
    m = (p1 + p2) / 2.0
    d_p1_m = d_p2_m = 0.0
    for i in range(len(p1)):
        if p1[i]:
            d_p1_m += p1[i] * math.log(p1[i] / m[i], 2)
        if p2[i]:
            d_p2_m += p2[i] * math.log(p2[i] / m[i], 2)
    jsd = (d_p1_m + d_p2_m) / 2.0
    assert 0 <= jsd <= 1, "Shannon-Jensen divergence should be between zero and one"
    return jsd


def Seed(seed):
    """Seeds the random number generator used for the MCMC.

    If you want to make the output of *InferEnrichmentMCMC* reproducible,
    you can seed the random number generator. Only call this
    method if you have confirmed that *PymcAvailable() == True*.
    The generator (*numpy.random.seed*) is set to *seed*,
    as is *random.seed*.
    """
    if not _pymcavailable:
        raise ImportError("pymc and numpy not both available")
    numpy.random.seed(seed)
    random.seed(seed)


def SiteEntropy(pi):
    """Calculates a site entropy (in bits) based on equilibrium fractions.

    This function takes as a calling argument a single dictionary *pi*
    which has the format of the return variable from *EquilibriumFracs*.
    It then calculates a site entropy *h* which is a measure of how
    variable the site is in terms of its tolerance for various amino acids.
    This entropy is defined as::

        h = - sum_aa pi[aa] * log2(pi[aa])

    where the sum is taken over all keys *aa* (presumably all amino acids)
    in *pi*, and *log2* indicates the logarithm to the base 2. The units
    on the returned site entropy are therefore bits.
    """
    h = 0.0
    for x in pi.itervalues():
        if x == 0:
            pass # 0 * log(0) is zero
        elif x < 0:
            raise ValueError("pi has an entry < 0, this should not happen!")
        else:
            h -= x * math.log(x, 2)
    return h


def EquilibriumFracs(wt, phi):
    """Calculates the equilibrium preference for an amino acid at a site.

    Most of the functions in this module are designed to calculate the
    enrichment ratio *phi* of some mutant amino acid relative to the wildtype
    amino acid at a site. Here we use the collection of all of these mutant
    amino acids enrichment ratios to calculate the equilibrium preference
    *pi* for each amino acid. 

    For each amino acid *aa* that is not equal to the wildtype amino acid
    *wt*, we define::

        pi[aa] = phi[aa] / (1.0 + sum_j phi[j])

    where the sum is taken over all non-wildtype amino acids *j*. For *wt*::

        pi[wt] = 1.0 / (1.0 + sum_j phi[j])

    In other words, *pi* is the expected probability that a site would have 
    some identity in a hypothetical equilibration to all of the enrichment
    ratios.

    CALLING VARIABLES:

    * *wt* : the one-letter amino acid code for the wildtype amino acid.

    * *phi* : a dictionary keyed be each of the non-wildtype amino acids,
      and with the values being the enrichment ratio for that amino acid.
      No strict checking is done to make sure that all of the non-wildtype
      amino acids are present, or even that the keys actually correspond
      to valid amino acid codes.

    RETURN VARIABLE:

    The return variable is the dictionary *pi*. It is keyed by all of the keys
    in *phi* (should be all non-wildtype amino acids) plus *wt* (should
    be the wildtype amino acid), and the values for each key is the 
    equilibrium preference for that amino acid.
    """
    if wt in phi:
        raise ValueError("wt is already in phi")
    denom = 1.0 + sum([x for x in phi.itervalues()])
    assert denom >= 1
    pi = {wt:1.0 / denom}
    for (aa, aaphi) in phi.iteritems():
        pi[aa] = aaphi / denom
    assert len(pi) == len(phi) + 1
    return pi


def DirectEnrichmentRatio(library_stats):
    """Directly calculates enrichment ratio from counts.

    The function *InferEnrichmentMCMC* infers the enrichment ratio
    *phi* from experimental counts in four different samples
    (*DNA*, *RNA*, *mutDNA*, and *mutvirus*) using a Bayesian approach.
    This function takes the much simpler approach of just calculating
    this ratio directly.

    It is simply done by calculating the frequencies of the mutation in
    each of the four samples, and then estimating *phi* as
    *(mutDNA - DNA) / (mutvirus - RNA)* for the different frequencies.
    If the denominator is zero and the numerator is zero, returns zero.
    Otherwise if the denominator is zero returns *float("inf")* if the
    numerator is > 0 and *-float("inf")* if the numerator is > 0.

    The calling variable is the list *library_stats* which should have
    the same properties as the calling variable of the same name
    to *InferEnrichmentMCMC*.

    The returned variable is the following 2-tuple: *(direct_ratio, mutdnacounts)*
    where *direct_ratio* is the ratio calculated as described above
    and *mutdnacounts* is the total number of counts for this mutations in the
    *mutDNA* library.

    """
    dnatotal = sum([libstats['Nrdna'] for libstats in library_stats])
    rnatotal = sum([libstats['Nrrna'] for libstats in library_stats])
    mutdnatotal = sum([libstats['Nrmutdna'] for libstats in library_stats])
    mutvirustotal = sum([libstats['Nrmutvirus'] for libstats in library_stats])
    dnacounts = sum([sum(libstats['nrdna_list']) for libstats in library_stats])
    rnacounts = sum([sum(libstats['nrrna_list']) for libstats in library_stats])
    mutdnacounts = sum([sum(libstats['nrmutdna_list']) for libstats in library_stats])
    mutviruscounts = sum([sum(libstats['nrmutvirus_list']) for libstats in library_stats])
    numerator = mutviruscounts / float(mutvirustotal) - rnacounts / float(rnatotal)
    denominator = mutdnacounts / float(mutdnatotal) - dnacounts / float(dnatotal)
    if denominator == numerator == 0:
        direct_ratio = 0
    elif denominator == 0:
        if numerator > 0:
            direct_ratio = float("inf")
        else:
            direct_ratio = -float("inf")
    else:
        direct_ratio = numerator / float(denominator)
    return (direct_ratio, mutdnacounts)


def InferEnrichmentMCMC(alpha, phi_prior, wtcodon, codons,\
        library_stats, nruns, nsteps, burn, thin, pickle=None,\
        loadfrompickle=False, convergence=None, plot_phi_traces=False,\
        plot_title='', progress_bar=False):
    """Infers enrichment ratio *phi* from MCMC over posterior.

    The function *InferEnrichmentMCMC_2* is preferred as it performs
    the same MCMC more quickly, particularly if ``scipy`` is available.

    Utilizes ``pymc`` to perform the MCMC.

    This function uses MCMC to infer the enrichment ratio at a specific site. 
    This enrichment ratio
    is the change in the frequency of a codon or amino acid in a selected library
    of mutants versus an unselected input library. Enrichment values > 1 indicate
    mutations that are favorable, and values < 1 indicate unfavorable. The values
    are calculated using four samples, named as follows:

        * *DNA* : this sample measures the sequencing error rate. Variants from
          wildtype in this library are due to sequencing errors.

        * *mutDNA* : this sample is the unselected library of mutants. Variants
          from wildtype in this library are due to the underlying rate of
          mutagenesis in the mutant library plus sequencing errors.

        * *RNA* : this sample measure any additional sources of error that affect
          the selected library. It is named based on what would be expected for
          selections performed with RNA viruses, where the genes must be reverse-
          transcribed which introduces additional errors. Variants from wildtype
          therefore arising from the sequencing errors plus these additional
          reverse-transcription errors.

        * *mutvirus* : this sample measures the frequency of variants in the
          selected library. Variants here arise due to selection on the initial
          mutants in the unselected library, and then there are the sequencing
          and reverse-transcription errors.

    For each non-wildtype codon at the site, the following rates are inferred:

        * *epsilon* : the sequencing error rate, as measured by variants in
          the *DNA* library.

        * *mu* : the mutagenesis rate, as measured by variants in the *mutDNA*
          library minus *epsilon*.

        * *rho* : the reverse-transcription error rate, as measured by variants
          in the *RNA* library minus *epsilon*.

        * *phi* : constrained to be the same for all codons in *codons*, this
          is the enrichment ratio in the *mutvirus* library relative to the
          *mutDNA* library after correcting for sequencing errors.

    The inference is done in a Bayesian fashion. A gamma distribution prior is
    used over each of the four unknown rates. The distributions all have the
    same shape parameter *alpha*, and the rate parameter is chosen so that
    the mean over the prior is equal to *phi_prior* for *phi*, and to the
    values specified in *library_stats* for that codon for *epsilon*, *mu*,
    and *rho*. The likelihood functions are Poisson likelihoods for observing
    the actual number of counts (specified in *library_stats* given the
    rates.

    The enrichment ratio *phi* for a given amino acid is calculated for
    multiple codons if *codons* has more than one entry. In that case,
    each mutant codon is constrained to have the same value of *phi*,
    but has different values of *epsilon*, *mu*, and *rho*. If there
    are multiple libraries in *library_stats*, then the enrichment is
    calculated for multiple libraries as well. Again, both libraries 
    share the same *phi*, but can have different values of the other
    three rate parameters for each codon.

    The MCMC over the posterior is performed using ``pymc``. The
    parameters *nruns*, *nsteps*, *burn*, and *thin* specify the MCMC
    parameters. You can use the *convergence* option to check for
    convergence.

    CALLING VARIABLES:

    * *alpha* is the gamma distribution shape parameter used for the
      priors. It must be a number > 0. The same shape parameter is
      used for all priors.

    * *phi_prior* is the prior estimate for the enrichment ratio *phi*.
      We use a gamma-distribution prior over *phi* with shape parameter
      *alpha* and mean *phi_prior*. It must be a number > 0.

    * *wtcodon* is a string giving the identity of the wildtype
      codon at the site for which we are computing the enrichment ratio.
      Should be upper case, for example 'GCA'.

    * *codons* is a list of the codons for the identity being inferred.
      For example, *['GGG', 'GGA', 'GGT', 'GGC']* for an inferring an
      enrichment ratio for glycine. There is not actually a strict
      enforcement that these represent the real genetic code,
      rather *codons* just must be a list with at least one entry that
      specifies the codon keys used for the statistics in *library_stats*.
      So if you want to infer the enrichment ratio for just one codon,
      just put that in *codons*. All strings should be upper case.

    * *library_stats* is a list of statistics for each independent library
      that is being used to calculate the enrichment *phi*. If there is only
      one library, then there will just be one entry in *library_stats*. Each
      entry itself in the list is a dictionary with the following keys:

      - 'mu_prior' : number > 0 giving prior estimate for the mutation rate
        in the mutDNA library. 

      - 'epsilon_prior' : a dictionary keyed by the integers 1, 2, and 3.
        *epsilon_prior[i]* gives the prior estimate for the sequencing
        error rate fo mutant codons with *i* nucleotide changes.

      - 'rho_prior' : a dictionary keyed by the integers 1, 2, and 3.
        *rho_prior[i]* gives the prior estimate for the reverse-transcription
        error rate for mutant codons with *i* nucleotide changes.

      - 'Nrdna' : integer giving the total number of counts for the DNA sample
        at the position in question.

      - 'Nrrna' : integer giving the total number of counts for the RNA sample
        at the position in question.

      - 'Nrmutdna' : integer giving the total number of counts for the mutDNA
        sample at the position in question.

      - 'Nrmutvirus' : integer giving the total number of counts for the 
        mutant virus sample at the position in question.

      - 'nrdna_list' : a list of the same length as *codons*, with entry *i*
        giving the number of counts for *codons[i]* in the DNA sample at the
        position in question.

      - 'nrrna_list' : a list of the same length as *codons*, with entry *i*
        giving the number of counts for *codons[i]* in the RNA sample at the
        position in question.

      - 'nrmutdna_list' : a list of the same length as *codons*, with entry *i*
        giving the number of counts for *codons[i]* in the mutDNA sample at the
        position in question.

      - 'nrmutvirus_list' : a list of the same length as *codons*, with entry *i*
        giving the number of counts for *codons[i]* in the mutvirus sample at the
        position in question.

    * *nruns* is the number of MCMC runs to perform, each starting from 
      different random values of the variables. Should be an integer >= 1.
      Each MCMC run is *nsteps*, with a burn-in of *burn* and thining of
      *thin*.

    * *nsteps* is an integer giving the number of MCMC steps to perform per run.
      Typically you might want values of 1e4 to 1e6 (depending on the chain
      is converging).

    * *burn* is a number giving the number of burn-in steps to perform per run
      before saving the MCMC values. Typically you might want a value
      equal to 10 to 20% of *nsteps*.

    * *thin* specifies the thinning of steps to perform (i.e. only sample
      every *thin* steps. Typicall you might want a value of 10 to 100.
      *thin* must be set so that *burn* and *nsteps* are both multiples
      of *thin*.

    * *pickle* specifies whether we save the entire ``pymc`` MCMC sampler using the
      Python *Pickle* database. It is *None* by default, which means
      that no such saving is done. If it is set to another value, then it
      should be the name of the file to which we want to save the sampler.
      This is the argument passed as *dbname* to the *pymc.MCMC* object
      created here. Typically you would want this file to have an extension
      of *.pickle*, although this is not required. If this file already exists,
      then it is overwritten if *loadfrompickle* is *False*. If *loadfrompickle*
      is not *False*, then no MCMC is run (see description of that option).

    * *loadfrompickle* is a Boolean switch that is meaningful only if *pickle*
      is set to some option other than *None*. In this case, if *loadfrompickle*
      is *True* and if the file specified by *pickle* already exists, then we
      do not re-run the MCMC but instead simply load the results already saved
      in the file *pickle*. Note that no checking is done to make sure that
      the calling arguments to this function match the ones used to create
      *pickle*, so only use this option if you know the origin of the data
      in *pickle*. However, this option can save you from re-running an MCMC
      whose results are already saved.

    * *convergence* is an option that lets you test for convergence of the
      MCMC. By default, it is *None*, meaning that no testing for convergence
      is done. However, if you set *nruns* to a value > 1, you can then test
      for MCMC convergence (this is recommended). Note however to setting
      *convergence* to a value other than *None* if *nruns* is 1 will
      raise an exception. To use this option for *nruns* > 1, set *convergence*
      to a number (typically you will want something between 1.0 and 1.1, a
      reasonable recommendation is 1.05). For the multiple runs, the Gelman-Rubin
      statistic for *phi* (the enrichment ratio)
      is computed (this statistic is theoretically 1.0 for perfect
      convergence). If this statistic is < *convergence*, then the chain is 
      considered to have converged. 

    * *plot_phi_traces* is an optional argument that can be used to generate
      PDF plots of the trace of *phi* as a function of the number of steps
      after thinning and burn-in (so there will be 900 points if *nsteps*
      is 1e4, *burn* is 1e3, and *thin* is 10) for each of the *nruns* runs.
      This PDF plot depends on ``pylab``, and so can only be created
      if this package is available (can be checked with 
      *mapmutsplot.PylabAvailable()*).
      By default *plot_phi_traces* is *None*, meaning no plot is created.
      If it is set to another value, it should be a string ending with the
      extension ``.pdf`` giving the name of the plot file to create. If you
      try to use this option when ``pylab`` is not available, an exception
      will be raised.

    * *plot_title* is an optional argument that is meaningful only if 
      *plot_phi_traces* is being used. In this case, *plot_title* is a string
      giving the name of the title placed above the plot in *plot_phi_traces*.

    * *progress_bar* is a Boolean switch specifying whether we display the
      MCMC progress bar.

    RETURN VARIABLES:

    This function returns the following 4-tuple: 
    *(phi_mean, phi_hpd95, phi_samples, converged)*
    In this tuple, the variables have the following meanings:

        - *phi_mean* is the mean enrichment ratio *phi* calculated
          over all of the *nruns* runs.

        - *phi_hpd95* is a tuple, list, or array with two entries:
          *phi_hpd95[0]* is the lower bound and *phi_hpd95[1]* is the
          upper bound of the highest posterior density 95% interval.

        - *phi_samples* is a list or array that gives all of the samples
          from the posterior distribution for the enrichment ratio *phi*
          after applying *burn* and *thin*.

        - *converged* specifies whether we have evidence that the
          chain converged. If *nruns* is 1 or *convergence* is *None*,
          this *converged* is *False*. Otherwise, *converged* is *True*
          if the criteria specified by *convergence* are met, and *False*
          if these criteria are not met.

    In addition, files might be created saving some of the MCMC results
    depending on the values of *pickle* and *plot_phi_traces*.

    """
    if not _pymcavailable:
        raise ImportError("Cannot import either numpy or pymc")
    # list of all MCMC variables
    variables = []
    # define phi (enrichment ratio) and its prior
    phi = pymc.Gamma(name='phi', alpha=alpha, beta=1.0 / phi_prior)
    variables.append(phi)
    # start adding variables for each library
    assert len(library_stats) >= 1
    ilibrary = 0
    for lstats in library_stats:
        ilibrary += 1
        icodon = 0
        # loop over mutant codons
        for codon in codons:
            ndiffs = len([i for i in range(len(codon)) if codon[i] != wtcodon[i]])
            assert 1 <= ndiffs <= 3, "%s and %s do not have 1 to 3 diffs" % (codon, wtcodon)
            # priors over rates 
            mu = pymc.Gamma(name='mu_lib%d_codon%s' % (ilibrary, codon), alpha=alpha, beta=1.0 / lstats['mu_prior'])
            rho = pymc.Gamma(name='rho_lib%d_codon%s' % (ilibrary, codon), alpha=alpha, beta=1.0 / lstats['rho_prior'][ndiffs])
            epsilon = pymc.Gamma(name='epsilon_lib%d_codon%s' % (ilibrary, codon), alpha=alpha, beta=1.0 / lstats['epsilon_prior'][ndiffs])
            variables += [mu, rho, epsilon]
            # likelihoods of counts for this codon
            ndna = pymc.Poisson(name='nrdna_codon%s_library%d' % (codon, ilibrary), mu=epsilon * lstats['Nrdna'], value=lstats['nrdna_list'][icodon], observed=True)
            nrna = pymc.Poisson(name='nrrna_codon%s_library%d' % (codon, ilibrary), mu=(epsilon + rho) * lstats['Nrrna'], value=lstats['nrrna_list'][icodon], observed=True)
            nmutdna = pymc.Poisson(name='nrmutdna_codon%s_library%d' % (codon, ilibrary), mu=(epsilon + mu) * lstats['Nrmutdna'], value=lstats['nrmutdna_list'][icodon], observed=True)
            nmutvirus = pymc.Poisson(name='nrmutvirus_codon%s_library%d' % (codon, ilibrary), mu=(epsilon + rho + phi * mu) * lstats['Nrmutvirus'], value=lstats['nrmutvirus_list'][icodon], observed=True)
            variables += [ndna, nrna, nmutdna, nmutvirus]
            icodon += 1
    assert len(variables) == len(library_stats) * 7 * len(codons) + 1
    assert thin < nsteps and burn < nsteps
    if not (nsteps % thin == 0 and burn % thin == 0):
        raise ValueError("nsteps = %d and burn = %d are not multiples of thin = %d" % (nsteps, burn, thin))
    n_per_run = (nsteps - burn) / thin # steps per run
    if pickle and loadfrompickle and os.path.isfile(pickle):
        mcmc = pymc.database.pickle.load(pickle)
    else:
        if pickle:
            mcmc = pymc.MCMC(variables, db='pickle', dbname=pickle)
        else:
            mcmc = pymc.MCMC(variables)
        for irun in range(nruns):
            mcmc.draw_from_prior()
            mcmc.sample(iter=nsteps, burn=burn, thin=thin, progress_bar=progress_bar)
        if pickle:
            mcmc.db.close()
    assert nruns * n_per_run == mcmc.stats(chain=None)['phi']['n']
    phi_means = [mcmc.stats(chain=irun)['phi']['mean'] for irun in range(nruns)]
    phi_samples = mcmc.trace('phi', chain=None)[:]
    phi_hpd95 = mcmc.stats(chain=None)['phi']['95% HPD interval']
    phi_mean = mcmc.stats(chain=None)['phi']['mean']
    converged = False
    if convergence != None:
        if nruns < 2:
            raise ValueError("You cannot use convergence unless nruns >= 2")
        if convergence < 1.0:
            raise ValueError("You cannot set convergence to less than 1.0")
        grstat = pymc.gelman_rubin(mcmc)['phi']
        if grstat <= convergence:
            converged = True
            for x in phi_means:
                if not (1.0 / convergence < x / phi_mean < convergence):
                    converged = False
                    break
    if plot_phi_traces:
        if not mapmuts.plot.PylabAvailable():
            raise ImportError("Cannot implement plot_phi_traces because pylab is not available")
        traces = [mcmc.trace('phi', chain=irun)[:] for irun in range(nruns)]
        mapmuts.plot.PlotTraces(traces, plot_phi_traces, xlabel='step number (after thinning and burn-in)', ylabel='enrichment ratio ($\phi$)', title=plot_title)
    return (phi_mean, phi_hpd95, phi_samples, converged)


def InferEnrichmentMCMC_2(alpha, phi_prior, wtcodon, codons, library_stats,\
        nruns, nsteps, burn, thin, pickle=None,\
        loadfrompickle=False, convergence=None, plot_phi_traces=False,\
        plot_title='', progress_bar=False):
    """Infers enrichment ratio *phi* from MCMC over posterior.

    This function is identical to *InferEnrichmentMCMC* in the calling
    and return arguments, so see that docstrings. However, this 
    function is faster because it places the variables in arrays
    which speeds the computations.

    It is especially superior if ``scipy`` is available, as
    ``scipy`` is then used to estimate maximum a posteriori
    values (in an empirical Bayes fashion) for key parameters
    prior to beginning the MCMC. Although it is not strictly
    required, it is HIGHLY ADVISED to install ``scipy``
    before running this function.
    """
    if not _pymcavailable:
        raise ImportError("Cannot import either numpy or pymc")
    nlibs = len(library_stats)
    assert nlibs >= 1
    ncodons = len(codons)
    # define phi (enrichment ratio) and its prior
    phi = pymc.Gamma(name='phi', alpha=alpha, beta=1.0 / phi_prior)
    # define epsilon, rho, and mu_rate: length = ncodons * nlibs
    # entries are icodon + ncodons * ilib
    # all are Gamma objects
    mu_betas = numpy.empty(ncodons * nlibs)
    epsilon_betas = numpy.empty(ncodons * nlibs)
    rho_betas = numpy.empty(ncodons * nlibs)
    for ilib in range(nlibs):
        for icodon in range(ncodons):
            codon = codons[icodon]
            ndiffs = len([i for i in range(len(codon)) if codon[i] != wtcodon[i]])
            assert 1 <= ndiffs <= 3, "%s and %s do not have 1 to 3 diffs" % (codon, wtcodon)
            mu_betas[icodon + ncodons * ilib] = 1.0 / library_stats[ilib]['mu_prior']
            epsilon_betas[icodon + ncodons * ilib] = 1.0 / library_stats[ilib]['epsilon_prior'][ndiffs]
            rho_betas[icodon + ncodons * ilib] = 1.0 / library_stats[ilib]['rho_prior'][ndiffs]
    mu = pymc.Gamma('mu', alpha, mu_betas)
    epsilon = pymc.Gamma('epsilon', alpha, epsilon_betas)
    rho = pymc.Gamma('rho', alpha, rho_betas)
    # define ntotal_dna, ntotal_rna, ntotal_mutdna, ntotal_mutvirus
    # lengths are ncodons * nlibs 
    # gives total number of reads at a site
    # entries are icodon + ncodons * ilib
    ntotal_dna = numpy.empty(ncodons * nlibs, dtype='int')
    ntotal_rna = numpy.empty(ncodons * nlibs, dtype='int')
    ntotal_mutdna = numpy.empty(ncodons * nlibs, dtype='int')
    ntotal_mutvirus = numpy.empty(ncodons * nlibs, dtype='int')
    for ilib in range(nlibs):
        ntotal_dna[ilib * ncodons : (ilib + 1) * ncodons] = library_stats[ilib]['Nrdna']
        ntotal_rna[ilib * ncodons : (ilib + 1) * ncodons] = library_stats[ilib]['Nrrna']
        ntotal_mutdna[ilib * ncodons : (ilib + 1) * ncodons] = library_stats[ilib]['Nrmutdna']
        ntotal_mutvirus[ilib * ncodons : (ilib + 1) * ncodons] = library_stats[ilib]['Nrmutvirus']
    # define mutrates_dna, mutrates_rna, mutrates_mutdna, mutrates_mutvirus
    # lengths are ncodons * nlibs 
    # entries are expected mutations for each sample
    # entries for ndna are icodon + ncodons * ilib
    @pymc.deterministic(plot=False)
    def mutrates_dna(epsilonx=epsilon):
        return epsilonx * ntotal_dna
    @pymc.deterministic(plot=False)
    def mutrates_rna(epsilonx=epsilon, rhox=rho):
        return (epsilonx + rhox) * ntotal_rna
    @pymc.deterministic(plot=False)
    def mutrates_mutdna(epsilonx=epsilon, mux=mu):
        return (epsilonx + mux) * ntotal_mutdna
    @pymc.deterministic(plot=False)
    def mutrates_mutvirus(epsilonx=epsilon, rhox=rho, mux=mu, phix=phi):
        return (epsilonx + rhox + phix * mux) * ntotal_mutvirus
    # define ncounts_dna, ncounts_rna, ncounts_mutdna, ncounts_mutvirus
    # lengths are ncodons * nlibs
    # entries are number of counts observed for each sample, Poisson likelihoods
    # entries are icodon + ncodons * ilib
    # nobserved gives the observed numbers
    nobserved_dna = numpy.empty(ncodons * nlibs, dtype='int')
    nobserved_rna = numpy.empty(ncodons * nlibs, dtype='int')
    nobserved_mutdna = numpy.empty(ncodons * nlibs, dtype='int')
    nobserved_mutvirus = numpy.empty(ncodons * nlibs, dtype='int')
    for ilib in range(nlibs):
        for icodon in range(ncodons):
            nobserved_dna[icodon + ncodons * ilib] = library_stats[ilib]['nrdna_list'][icodon]
            nobserved_rna[icodon + ncodons * ilib] = library_stats[ilib]['nrrna_list'][icodon]
            nobserved_mutdna[icodon + ncodons * ilib] = library_stats[ilib]['nrmutdna_list'][icodon]
            nobserved_mutvirus[icodon + ncodons * ilib] = library_stats[ilib]['nrmutvirus_list'][icodon]
    ncounts_dna = pymc.Poisson('ncounts_dna', mu=mutrates_dna,  value=nobserved_dna, observed=True)
    ncounts_rna = pymc.Poisson('ncounts_rna', mu=mutrates_rna,  value=nobserved_rna, observed=True)
    ncounts_mutdna = pymc.Poisson('ncounts_mutdna', mu=mutrates_mutdna,  value=nobserved_mutdna, observed=True)
    ncounts_mutvirus = pymc.Poisson('ncounts_mutvirus', mu=mutrates_mutvirus,  value=nobserved_mutvirus, observed=True)
    variables = [ncounts_dna, ncounts_rna, ncounts_mutdna, ncounts_mutvirus, epsilon, rho, mu, phi]
    assert thin < nsteps and burn < nsteps
    if not (nsteps % thin == 0 and burn % thin == 0):
        raise ValueError("nsteps = %d and burn = %d are not multiples of thin = %d" % (nsteps, burn, thin))
    n_per_run = (nsteps - burn) / thin # steps per run
    if pickle and loadfrompickle and os.path.isfile(pickle):
        mcmc = pymc.database.pickle.load(pickle)
    else:
        if pickle:
            mcmc = pymc.MCMC(variables, db='pickle', dbname=pickle)
        else:
            mcmc = pymc.MCMC(variables)
        for irun in range(nruns):
            if _scipyavailable:
                # initialize with maximum a posterior estimates 
                epsilon_map = pymc.MAP([epsilon, ncounts_dna])
                epsilon_map.fit()
                rho_map = pymc.MAP([rho, ncounts_rna])
                rho_map.fit()
                mu_map = pymc.MAP([mu, ncounts_mutdna])
                mu_map.fit()
                phi_map = pymc.MAP([phi, ncounts_mutvirus])
                phi_map.fit()
            else:
                mcmc.draw_from_prior()
            # now run MCMC
            mcmc.sample(iter=nsteps, burn=burn, thin=thin, progress_bar=progress_bar)
        if pickle:
            mcmc.db.close()
    assert nruns * n_per_run == mcmc.stats(chain=None)['phi']['n']
    phi_means = [mcmc.stats(chain=irun)['phi']['mean'] for irun in range(nruns)]
    phi_samples = mcmc.trace('phi', chain=None)[:]
    phi_hpd95 = mcmc.stats(chain=None)['phi']['95% HPD interval']
    phi_mean = mcmc.stats(chain=None)['phi']['mean']
    converged = False
    if convergence != None:
        if nruns < 2:
            raise ValueError("You cannot use convergence unless nruns >= 2")
        if convergence < 1.0:
            raise ValueError("You cannot set convergence to less than 1.0")
        grstat = pymc.gelman_rubin(mcmc)['phi']
        if grstat <= convergence:
            converged = True
    if plot_phi_traces:
        if not mapmuts.plot.PylabAvailable():
            raise ImportError("Cannot implement plot_phi_traces because pylab is not available")
        traces = [mcmc.trace('phi', chain=irun)[:] for irun in range(nruns)]
        mapmuts.plot.PlotTraces(traces, plot_phi_traces, xlabel='step number (after thinning and burn-in)', ylabel='enrichment ratio ($\phi$)', title=plot_title)
    return (phi_mean, phi_hpd95, phi_samples, converged)


def CredibleInterval(a, interval):
    """Computes the median-centered credible interval from a ``numpy`` array of samples.

    *a* is a numpy array of dimension 1 that lists samples from the posterior.

    *interval* is a number <= 1 and >= 0 that gives the credible interval. 
    For example, 0.95 corresponds to a 95% credible interval.

    This function returns a ``numpy`` array of length 2 giving the lower and
    upper bounds to the credible interval.
    """
    if len(a.shape) != 1:
        raise ValueError("a is not of dimension 1")
    assert 0 <= interval <= 1, "Invalid value of interval"
    acopy = a.copy()
    acopy.sort()
    return numpy.array([acopy[int((0.5 - interval / 2.) * len(acopy))], acopy[int((0.5 + interval / 2.) * len(acopy))]])


def InferPreferencesMCMC(library_stats, pi_concentration, epsilon_concentration, mu_concentration, rho_concentration, nruns, nsteps, burn, npreburns, thin, progress_bar=False, minvalue=1e-7, debugging=False):
    """Infers equilibrium amino-acid preferences :math:`\pi_{r,a}` by MCMC.

    This function utilizes ``pymc`` for the inference. 

    This function infers the equilibrium preference :math:`\pi_{r,a}` for
    all amino acids *a* for a given site *r*. There are 21 possible values
    for *a* since we allow all 20 amino acids plus a stop codon (denoted
    by the * character). The preferences sum to one:
    :math:`\sum_a \pi_{r,a} = 1`. Higher values of :math:`\pi_{r,a}`
    correspond to a larger preference for *a* at site *r*.

    The full inference algorithm is described in the documentation for
    the ``mapmuts`` function, and so it not recapitulated here.

    CALLING VARIABLES:

    * *library_stats* is a list of statistics for each replicate of
      the libraries. There must be at least one replicate, but there
      can be more. Each entry in the list is a dictionary giving the
      statistics for that replicate. This dictionary has the following
      strings for keys:

      * *wtcodon* is a string giving the wildtype codon at the site. Should
        be an upper-case string, such as *GGC*.

      * *mu_prior* : a number giving the prior estimate for the probability
        that the site is mutated from its wildtype codon to some specific other
        codon in the *mutDNA* library. Note that this is the probability to
        be mutated to a **specific** codon, and so is the overall
        codon mutation rate divided by 63 (since there are 63 other codons).

      * *epsilon_prior* : a dictionary keyed by the integers 1, 2, and 3
        giving the prior estimate for the probability that the site is
        erroneously read as some other **specific** codon in the *DNA*
        library. The value for the key of 1 is the probability to be
        mutated to some other specific codon that differs by 1 nucleotide,
        the value for the key of 2 is the probability to be mutated to 
        some other specific codon that differs by 2 nucleotides, etc.

      * *rho_prior* : like *epsilon_prior*, but the probabilities for
        reverse-transcription errors in the *RNA* library.

      * *nrdna_counts* : keyed by all 64 codons (upper-case strings), with values
        being the number of codons with that identity in the *DNA* library.
        The sum of entries should equal *Nrdna*.

      * *nrrna_counts* : like *nrdna_counts* but for the *RNA* library.

      * *nrmutdna_counts* : like *nrdna_counts* but for the *mutDNA* library.

      * *nrmutvirus_counts* : like *nrdna_counts* but for the *mutvirus* library.

    * *pi_concentration* is the concentration parameter (:math:`\sigma_{\pi}`)
      for the symmetric Dirichlet-distribution prior over the :math:`\pi_{r,a}` 
      values.

    * *epsilon_concentration* is the concentration parameter
      (:math: `\sigma_{\epsilon}`) for the Dirichlet-distribution prior
      over error rate in the *DNA* library.

    * *mu_concentration* is the concentration parameter (:math:`\sigma_{\mu}`)
      for the Dirichlet-distribution prior over the mutagenesis rate in
      the *mutDNA* library.

    * *rho_concentration* is the concentration parameter (:math:`\sigma_{\\rho}`)
      for the Dirichlet-distribution prior over the reverse-transcription error
      rata in the *RNA* library.

    * *nruns* is the number of MCMC runs to perform, each starting from 
      different random values of the variables. Should be an integer >= 1.
      Each MCMC run is *nsteps*, with a burn-in of *burn* and thining of
      *thin*.

    * *nsteps* is an integer giving the number of MCMC steps to perform per run.
      Typically you might want values of 1e4 to 1e6 (depending on the chain
      is converging).

    * *burn* is a number giving the number of burn-in steps to perform per run
      before saving the MCMC values. Typically you might want a value
      equal to 10 to 20% of *nsteps*.

    * *npreburns* is the number of "pre-burn" runs that are performed. 
      Setting this variable to at least one, and preferably two, is 
      **strongly** recommended. The MCMC becomes very inefficient as it
      approaches its maximum because the sampler currently does not choose
      the step sizes in the Dirichlet-distributed variables in a particularly
      intelligent way. This option specifies that *npreburns* MCMC runs of
      *burn* steps are performed before the main MCMC, with the results
      of these pre-burn runs used to tune the step sizes for the main
      MCMC. The other heuristic (which is always applied) is to place
      the wildtype codon as the last entry in the Dirichlet-distributed
      variable since it will typically be the largest entry. However, 
      it is still not clear if this is the most efficient way to do
      MCMC updates of Dirichlet variables, and it is possible that
      using a delta exchange type operator (as implemented in ``BEAST``)
      might be substantially more efficient). In addition,
      the setting of the wildtype codon to the last entry is only done
      based on the entries in the first library in *library_stats* --
      so if there are multiple libraries, this efficiency is lost.

    * *thin* specifies the thinning of steps to perform (i.e. only sample
      every *thin* steps. Typically you might want a value of 10 to 100.
      *thin* must be set so that *burn* and *nsteps* are both multiples
      of *thin*.

    * *progress_bar* is a Boolean switch specifying whether we display the 
      MCMC progress bar.

    * *minvalue* is the minimum value assigned to any of the priors (over *mu*,
      *rho*, *epsilon*) over the rates. This prevents the priors from being set to
      zero. A reasonable value is *1e-7* or *1e-8*. Set to *1e-7* by default.

    * *debugging* is a Boolean switch specifying whether we display
      fairly extensive debugging information. Is *False* by
      default. If *True*, the following is done:

        - MCMC is run with *verbose=2* rather than the default
          of *verbose=0*.

        - Information about the pre-burn runs and the step methods
          are printed to standard output.

    RETURN VALUE:

    This function returns the following tuple:
    *(pi_mean, pi_cred95, pi_traces, run_diff)*

    The tuple-elements are as follows:

    * *pi_mean* is a ``numpy`` array of length 21. Element *pi_mean[i]* holds
      the mean :math:`\pi_{r,a}` value over all non-burnin steps and all
      *nruns* runs for amino-acid *a*, where 
      *i = mapmuts.sequtils.AminoAcids(includestop=True).index(a)*.

    * *pi_cred95* gives the median-centered 95% credible interval.
      *pi_cred95[i]* is the interval corresponding to the estimate with 
      the mean of *pi_mean[i]*. The interval is represented such
      that *pi_cred95[i][0]* is the lower bound, and *pi_cred95[i][1]*
      is the upper bound. Overall, *pi_cred95* is a ``numpy`` array
      of dimension *(21, 2)*.

    * *pi_traces* is a ``numpy`` array of dimension 
      *(nruns, (nsteps - burn) / thin, 21)*. This give values of the
      :math:`\pi_{r,a}` values during the thinned non-burnin MCMC steps.
      Specifically, *pi_traces[irun][istep][iaa]* gives the value
      of :math:`\pi_{r,a}` for run *irun* (where *0 <= irun < nruns*)
      for step *istep* (where *0 <= istep < (nsteps - burn) / thin*)
      for amino acid *a* where 
      *iaa = mapmuts.sequtils.AminoAcids(includestop=True).index(a)*.

    * *run_diff* is a measure of the difference between the mean 
      :math:`\pi_{r,a}` values for all non-burnin steps between the
      *nruns* different runs. If *nruns* is 1, then this cannot be
      calculated, so *run_diff = None*. Otherwise, we calculate
      the Shannon-Jensen divergence between the :math:`\pi_{r,a}`
      values for each pair of runs, taking the logarithm to the base
      2. If the two runs are identical, this divergence is zero,
      while its maximum value is one. We then calculate the average
      of the pairwise Shannon-Jensen divergence for all pairs of
      runs, and return that as *run_diff*. If *run_diff* is close
      to zero (say less than 0.01), that indicates that the runs are converging
      to similar values. If *run_diff* is large (say > 0.05), that indicates
      they are converging to different values -- in that case,
      you might want to increase *nsteps* (and perhaps *burn*).
    """
    if not _pymcavailable:
        raise ImportError("Cannot import either numpy or pymc")
    if debugging:
        verbose = 2
    else:
        verbose = 0
    assert isinstance(npreburns, int) and npreburns >= 0, "npreburns must be integer >= 0"
    nlibs = len(library_stats)
    assert nlibs >= 1, "No libraries specified"
    assert minvalue > 0, 'minvalue must be > 0'
    wtcodon = library_stats[0]['wtcodon'] # wildtype codon in the first library.
    wtaa = mapmuts.sequtils.Translate([('wt', wtcodon)])[0][1]
    if not wtaa:
        wtaa = "*"
    # We re-order codons and aas so that the wtcodon is last, and therefore
    # is the incompleted entry in the Dirichlet distributions. This should
    # improve MCMC performance. Note that this improvement will be lost
    # if not all libraries have the same wildtype codon -- further re-writing
    # of this code would be necessary to optimize for that possibility.
    codons = [codon for codon in mapmuts.sequtils.Codons() if codon != wtcodon] + [wtcodon] # list of codons with wtcodon last
    ncodons = len(codons)
    assert ncodons == len(mapmuts.sequtils.Codons())
    aas_original = mapmuts.sequtils.AminoAcids(includestop=True)
    aas = [aa for aa in aas_original if aa != wtaa] + [wtaa]
    naas = len(aas)
    assert len(aas_original) == naas
    codon_to_aa_indices = numpy.ndarray(ncodons, dtype='int')
    icodon = 0
    for codon in codons:
        aa = mapmuts.sequtils.Translate([('wt', codon)])[0][1]
        if not aa:
            aa = '*'
        assert aa in aas, "Failed to find aa %s" % aa
        codon_to_aa_indices[icodon] = aas.index(aa)
        icodon += 1
    # define pi, the equilibrium preferences, with Dirichlet prior
    pi_incomplete = pymc.Dirichlet('pi_incomplete', pi_concentration * numpy.ones(naas), value=numpy.array([1.0 / naas] * (naas - 1))) # amino acid version of pi, not completed. Important to set value to uniform or convergence will be poor for entries initially seeded with tiny values.
    pi = pymc.CompletedDirichlet('pi', pi_incomplete) # amino acid version of pi
    @pymc.deterministic(plot=False)
    def c_pi(pix=pi): 
        """Codon version of pi."""
        return pix.take(codon_to_aa_indices)
    # define epsilon, rho, and mu with Dirichlet priors for each library
    alpha_mu = {}
    alpha_rho = {}
    alpha_epsilon = {}
    mu_incomplete = {}
    rho_incomplete = {}
    epsilon_incomplete = {}
    mu = {}
    rho = {}
    epsilon = {}
    for ilib in range(nlibs):
        if codons.count(library_stats[ilib]['wtcodon']) != 1: 
            raise ValueError("Failed to find wtcodon of %s in codons" % library_stats[ilib]['wtcodon'])
        iwtcodon = codons.index(library_stats[ilib]['wtcodon']) # index of wildtype codon
        alpha_mu[ilib] = numpy.zeros(ncodons)
        alpha_rho[ilib] = numpy.zeros(ncodons)
        alpha_epsilon[ilib] = numpy.zeros(ncodons)
        for icodon in range(ncodons):
            codon = codons[icodon]
            ndiffs = len([i for i in range(len(codon)) if codon[i] != library_stats[ilib]['wtcodon'][i]])
            if icodon != iwtcodon:
                assert ndiffs > 0, "Can't be wildtype if not different codon"
                alpha_mu[ilib][icodon] = max(minvalue, library_stats[ilib]['mu_prior'])
                alpha_rho[ilib][icodon] = max(minvalue, library_stats[ilib]['rho_prior'][ndiffs])
                alpha_epsilon[ilib][icodon] = max(minvalue, library_stats[ilib]['epsilon_prior'][ndiffs])
        alpha_mu[ilib][iwtcodon] = 1 - numpy.sum(alpha_mu[ilib])
        alpha_rho[ilib][iwtcodon] = 1 - numpy.sum(alpha_rho[ilib])
        alpha_epsilon[ilib][iwtcodon] = 1 - numpy.sum(alpha_epsilon[ilib])
        assert alpha_mu[ilib][iwtcodon] > 0, "Non-wildtype mu_prior sums to >= 1"
        assert alpha_rho[ilib][iwtcodon] > 0, "Non-wildtype rho_prior sums to >= 1"
        assert alpha_epsilon[ilib][iwtcodon] > 0, "Non-wildtype epsilon_prior sums to >= 1"
        mu_incomplete[ilib] = pymc.Dirichlet('mu_incomplete_lib%d' % ilib, mu_concentration * ncodons * alpha_mu[ilib], value=alpha_mu[ilib][ : -1])
        mu[ilib] = pymc.CompletedDirichlet('mu_lib%d' % ilib, mu_incomplete[ilib])
        rho_incomplete[ilib] = pymc.Dirichlet('rho_incomplete_lib%d' % ilib, rho_concentration * ncodons * alpha_rho[ilib], value=alpha_rho[ilib][ : -1])
        rho[ilib] = pymc.CompletedDirichlet('rho_lib%d' % ilib, rho_incomplete[ilib])
        epsilon_incomplete[ilib] = pymc.Dirichlet('epsilon_incomplete_lib%d' % ilib, epsilon_concentration * ncodons * alpha_epsilon[ilib], value=alpha_epsilon[ilib][ : -1])
        epsilon[ilib] = pymc.CompletedDirichlet('epsilon_lib%d' % ilib, epsilon_incomplete[ilib])
    # define the likelihoods for each library
    pr_nrdna = {}
    pr_nrrna = {}
    pr_nrmutdna = {}
    pr_nrmutvirus = {}
    nrdna = {}
    nrrna = {}
    nrmutdna = {}
    nrmutvirus = {}
    normalization = {}
    def normalization_eval(c_pi_x, mu_x):
        """Evaluates normalization for nrmutvirus."""
        return numpy.dot(c_pi_x, mu_x.ravel())
    for ilib in range(nlibs):
        delta = numpy.zeros(ncodons)
        iwtcodon = codons.index(library_stats[ilib]['wtcodon']) # index of wildtype codon
        delta[iwtcodon] = 1.0
        for (d, label) in [(nrdna, 'nrdna_counts'), (nrrna, 'nrrna_counts'), (nrmutdna, 'nrmutdna_counts'), (nrmutvirus, 'nrmutvirus_counts')]:
            d[ilib] = numpy.empty(ncodons, dtype='int')
            for icodon in range(ncodons):
                d[ilib][icodon] = library_stats[ilib][label][codons[icodon]]
        pr_nrdna[ilib] = pymc.Multinomial('pr_nrdna_lib%d' % ilib, n=numpy.sum(nrdna[ilib]), p=epsilon[ilib], value=nrdna[ilib], observed=True)
        pr_nrrna[ilib] = pymc.Multinomial('pr_nrrna_lib%d' % ilib, n=numpy.sum(nrrna[ilib]), p=epsilon[ilib] + rho[ilib] - delta, value=nrrna[ilib], observed=True)
        pr_nrmutdna[ilib] = pymc.Multinomial('pr_nrmutdna_lib%d' % ilib, n=numpy.sum(nrmutdna[ilib]), p=epsilon[ilib] + mu[ilib] - delta, value=nrmutdna[ilib], observed=True)
        normalization[ilib] = pymc.Deterministic(eval=normalization_eval, name='normalization_lib%d' % ilib, parents={'c_pi_x':c_pi, 'mu_x':mu[ilib]}, doc='Evaluates numpy.dot(c_pi, mu[ilib])')
        pr_nrmutvirus[ilib] = pymc.Multinomial('pr_nrmutvirus_lib%d' % ilib, n=numpy.sum(nrmutvirus[ilib]), p=epsilon[ilib] + rho[ilib] + c_pi * mu[ilib] / normalization[ilib] - 2 * delta, value=nrmutvirus[ilib], observed=True)
    variables = [pi_incomplete, pi, c_pi]
    for ilib in range(nlibs):
        for d in [epsilon_incomplete, epsilon, mu_incomplete, mu, rho_incomplete, rho, pr_nrdna, pr_nrrna, pr_nrmutdna, normalization, pr_nrmutvirus]:
            variables.append(d[ilib])
    assert thin < nsteps and burn < nsteps, "nsteps must be greater than both thin and burn"
    if not (nsteps % thin == 0 and burn % thin == 0):
        raise ValueError("nsteps = %d and burn = %d are not multiples of thin = %d" % (nsteps, burn, thin))
    pi_traces = [] # save traces of pi values for runs
    initial_stochastic_values = {} # initial values of stochastics before any burn-in. Reset to these prior to each run.
    for irun in range(nruns):
        # To avoid carry-over from previous runs, we set each stochastic
        # variable back to its very initial values, which are just the initial 
        # estimates before MCMC.
        mcmc = pymc.MCMC(variables, verbose=0)
        for stochasticvariable in mcmc.step_method_dict.iterkeys():
            if stochasticvariable in initial_stochastic_values:
                stochasticvariable.value = copy.deepcopy(initial_stochastic_values[stochasticvariable])
                if debugging:
                    print "Setting initial value of stochastic variable %s to %s." % (stochasticvariable, str(stochasticvariable.value))
            else:
                initial_stochastic_values[stochasticvariable] = copy.deepcopy(stochasticvariable.value)
        # The proposal_sd after each pre-burn run should be the values at the end
        # of the previous pre-burn run. This will make the steps for the Dirichlet
        # variable much more efficient, since small-valued variables will have
        # small steps and large-valued variables with have large steps.
        for ipreburn in range(npreburns):
            if debugging:
                print "\nPerforming pre-burn run %d" % ipreburn
            mcmc = pymc.MCMC(variables, verbose=verbose)
            for stochasticvariable in mcmc.step_method_dict.iterkeys():
                mcmc.use_step_method(pymc.Metropolis, stochasticvariable, verbose=verbose)
                if debugging:
                    print "For stochastic variable %s:\n\tinitial values = %s\n\tinitial proposal_sd = %s" % (stochasticvariable, str(stochasticvariable.value), str(mcmc.step_method_dict[stochasticvariable][0].proposal_sd))
            mcmc.sample(iter=burn, burn=0, thin=thin, progress_bar=progress_bar)
            if debugging:
                print "At completion of pre-burn run %d:" % ipreburn
                for stochasticvariable in mcmc.step_method_dict.iterkeys():
                    print "For stochastic variable %s:\n\tfinal values = %s\n\tfinal proposal_sd = %s" % (stochasticvariable, str(stochasticvariable.value), str(mcmc.step_method_dict[stochasticvariable][0].proposal_sd))
        # now set up the actual sampling run
        mcmc = pymc.MCMC(variables, verbose=verbose)
        if debugging:
            print "\nNow beginning the actual sampling MCMC."
        for stochasticvariable in mcmc.step_method_dict.iterkeys():
            mcmc.use_step_method(pymc.Metropolis, stochasticvariable, verbose=verbose)
            if debugging:
                print "For stochastic variable %s:\n\tinitial values = %s\n\tinitial proposal_sd = %s" % (stochasticvariable, str(stochasticvariable.value), str(mcmc.step_method_dict[stochasticvariable][0].proposal_sd))
        mcmc.sample(iter=nsteps, burn=burn, thin=thin, progress_bar=progress_bar)
        i_trace = mcmc.trace('pi', chain=None)[:][:,0,:]
        if debugging:
            print "At completion of the actual sampling MCMC."
            for stochasticvariable in mcmc.step_method_dict.iterkeys():
                print "For stochastic variable %s:\n\tfinal values = %s\n\tfinal proposal_sd = %s" % (stochasticvariable, str(stochasticvariable.value), str(mcmc.step_method_dict[stochasticvariable][0].proposal_sd))
        pi_traces.append(i_trace)
    merged_traces = numpy.concatenate(pi_traces)
    pi_mean = numpy.mean(merged_traces, axis=0) 
    pi_cred95 = numpy.array([CredibleInterval(merged_traces.transpose()[iaa], 0.95) for iaa in range(naas)])
    pi_traces = numpy.array(pi_traces)
    # Now reorder the values in pi_mean, pi_cred95, pi_traces so that the amino-acid
    # indexing is changed in that from aas to that in aas_original
    reindex = numpy.array([aas.index(aa) for aa in aas_original])
    pi_mean = pi_mean.take(reindex)
    assert pi_cred95.shape == (naas, 2), "pi_cred95 has invalid shape of %s" % str(pi_cred95.shape)
    pi_cred95 = pi_cred95.take(reindex, axis=0)
    assert pi_cred95.shape == (naas, 2), "Re-indexed pi_cred95 has invalid shape of %s" % str(pi_cred95.shape)
    assert pi_traces.shape == (nruns, (nsteps - burn) / thin, naas), "pi_traces has invalid shape of %s" % pi_traces.shape
    pi_traces = pi_traces.take(reindex, axis=2)
    assert pi_traces.shape == (nruns, (nsteps - burn) / thin, naas), "Re-indexed pi_traces has invalid shape of %s" % pi_traces.shape
    # Calculate the run_diff
    if nruns == 1:
        run_diff = None
    else:
        run_diff = []
        for irun in range(nruns):
            for jrun in range(irun + 1, nruns):
                i_pi_mean = numpy.mean(pi_traces[irun], axis=0)
                j_pi_mean = numpy.mean(pi_traces[jrun], axis=0)
                run_diff.append(ShannonJensenDivergence(i_pi_mean, j_pi_mean))
        run_diff = numpy.mean(run_diff)
    # Return the results with the reordered variables
    return (pi_mean, pi_cred95, pi_traces, run_diff)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
