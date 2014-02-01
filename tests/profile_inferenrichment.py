"""Profiles the enrichment ratio inference.

Figures out what is taking the CPU time.

Written by Jesse Bloom, 2013.
"""

import time
import sys
import cProfile
import pstats
import mapmuts.simulate
import mapmuts.plot
import mapmuts.bayesian


def DoInferences():
    """Do the inference on simulated data."""
    mapmuts.simulate.Seed(1)
    mapmuts.bayesian.Seed(1)
    phis = [0.1, 0.6, 2.7]
    alpha = 4.0
    phi_prior = 0.5
    mu = 9.7e-5
    rho = {1:2.1e-5, 2:5.6e-7, 3:1.0e-7}
    epsilon = {1:6.4e-5, 2:3.2e-7, 3:1.5e-7}
    wtcodon = 'TCA'
    codons = ['TTA', 'TTG', 'CTT', 'CTC', 'CTG', 'CTA']
    libsize = 6e5
    nlibs = 2
    perturb = 1.5
    convergence = 1.1
    nruns = 3
    nsteps = [2e4]
    thin = 10
    burnfrac = 0.2
    print "\nTesting MCMC enrichment value inference (this may take a while)."
    for phi in phis:
        print "\nTesting enrichment ratio inference for phi = %.2f..." % phi
        library_stats = mapmuts.simulate.SimulateCounts(phi, mu, epsilon, rho, wtcodon, codons, libsize, nlibs, perturb)
        for nstep in nsteps:
            print "For %d steps..." % (nstep),
            (phi_mean, phi_hpd95, phi_samples, converged) = mapmuts.bayesian.InferEnrichmentMCMC_2(alpha, phi_prior, wtcodon, codons, library_stats, nruns, nstep, int(nstep * burnfrac), thin, convergence=convergence)
            print "inferred phi = %.3f (%.3f to %.3f) compared to true value of %.3f, converged=%s." % (phi_mean, phi_hpd95[0], phi_hpd95[1], phi, converged)


# main body of script
if not mapmuts.bayesian.PymcAvailable():
    print "FAILED: Cannot do profiling if pymc / numpy are not available."
else:
    profilefile = 'pstats-inferenrichment'
    print "Profiling inference of enrichment ratios to %s..." % profilefile
    start = time.clock()
    sys.stdout.flush()
    cProfile.run('DoInferences()', profilefile)
    t = time.clock() - start
    print "Completed inferences and profiling in %.3f seconds." % t
    s = pstats.Stats(profilefile)
    s.strip_dirs()
    s.sort_stats('time')
    print "Here are the profiling stats:"
    s.print_stats()
