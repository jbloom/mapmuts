"""Tests Bayesian inference of enrichment ratios.

Written by Jesse Bloom, 2013.
"""

import time
import sys
import unittest
import mapmuts.simulate
import mapmuts.plot
import mapmuts.bayesian


class TestInferEnrichment(unittest.TestCase):
    """Tests Bayesian inference of enrichment ratios."""

    def test_InferEnrichment(self):
        """Tests inference of enrichment ratios.
        
        Uses a large library size such that inferred enrichment should
        be very close to the true value.
        """
        mapmuts.simulate.Seed(1)
        mapmuts.bayesian.Seed(1)
        phis = [0.3, 3.0]
        alpha = 4.0
        phi_prior = 0.5
        mu = 9.7e-5
        rho = {1:2.1e-5, 2:5.6e-7, 3:1.0e-7}
        epsilon = {1:6.4e-5, 2:3.2e-7, 3:1.5e-7}
        wtcodon = 'TCA'
        codons = ['TTA', 'TTG', 'CTT', 'CTC', 'CTG', 'CTA']
        libsize = 1e8
        nlibs = 2
        perturb = 1.5
        convergence = 1.1
        nruns = 3
        nsteps = [4e4]
        thin = 10
        burnfrac = 0.2
        match_criterion = 1.1 # require actual and inferred to be at least this close
        sys.stderr.write("\nTesting MCMC enrichment value inference (this may take a while).")
        for phi in phis:
            sys.stderr.write("\nTesting enrichment ratio inference for phi = %.2f..." % phi)
            library_stats = mapmuts.simulate.SimulateCounts(phi, mu, epsilon, rho, wtcodon, codons, libsize, nlibs, perturb)
            for nstep in nsteps:
                sys.stderr.write("\nFor %d steps..." % (nstep))
                (phi_mean, phi_hpd95, phi_samples, converged) = mapmuts.bayesian.InferEnrichmentMCMC_2(alpha, phi_prior, wtcodon, codons, library_stats, nruns, nstep, int(nstep * burnfrac), thin, convergence=convergence)
                sys.stderr.write("inferred phi = %.3f (%.3f to %.3f) compared to true value of %.3f, converged=%s." % (phi_mean, phi_hpd95[0], phi_hpd95[1], phi, converged)) 
                msg = 'FAILED: the enrichment ratio inference did not converge for phi = %.4f' % phi
                if not converged:
                    sys.stderr.write("\n%s" % msg)
                self.assertTrue(converged, msg)
                msg = 'FAILED: inferred enrichment ratio of %.4f was not close to actual of %.4f' % (phi_mean, phi)
                if not (1.0 / match_criterion < phi_mean / phi < match_criterion):
                    sys.stderr.write("\n%s" % msg)
                self.assertTrue(1.0 / match_criterion < phi_mean / phi < match_criterion, msg)
        sys.stderr.write("\nEnrichment ratio inference testing complete.\n")


if __name__ == '__main__':
    if not mapmuts.bayesian.PymcAvailable():
        sys.stderr.write("\nFAILED: Cannot test enrichment inference without numpy / pymc.\n")
    else:
        runner = unittest.TextTestRunner()
        unittest.main(testRunner=runner)
