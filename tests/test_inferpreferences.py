"""Tests Bayesian inference of equilibrium preferences.

Written by Jesse Bloom, 2013.
"""

import time
import sys
import os
import unittest
import copy
import mapmuts.simulate
import mapmuts.plot
import mapmuts.bayesian


class TestInferPreferences(unittest.TestCase):
    """Tests Bayesian inference of equilibrium preferences."""

    def test_InferPreferences(self):
        """Tests inference of equilibrium preferences.
        
        Uses a large library size such that inferred preferences should
        be very close to the true value.
        """
        import numpy
        numpy.random.seed(0)
        aas = mapmuts.sequtils.AminoAcids(includestop=True)
        wtcodon = 'GTA'
        wtaa = mapmuts.sequtils.Translate([('head', wtcodon)])[0][1]
        alphas = numpy.random.uniform(1e-7, 1.0, len(aas))
        alphas[aas.index(wtaa)] = 2.0
        pi = dict(zip(aas, numpy.random.dirichlet(alphas)))
        pivec = numpy.array([pi[aa] for aa in aas])
        mu = 9.7e-5
        rho = {1:2.1e-5, 2:5.6e-7, 3:1.0e-7}
        epsilon = {1:6.4e-5, 2:3.2e-7, 3:1.5e-7}
        libsize = 1e9 # large value to reduce sampling error
        nlibs = 2
        perturb = 1.5
        nruns = 2
        npreburns = 2
        nsteps = [3e4]
        burnfrac = 0.15
        match_criterion = 0.001 # require actual and inferred to be at least this close 
        sys.stderr.write("\nTesting MCMC equilibrium preference inference (this may take a while).")
        mapmuts.simulate.Seed(1)
        library_stats = mapmuts.simulate.SimulatePreferences(pi, mu, epsilon, rho, wtcodon, nlibs, libsize, perturb)
        for nstep in nsteps:
            sys.stderr.write("\n\nPerforming two replicates of preference inference of %d steps... " % (nstep))
            sys.stdout.flush()
            pi_by_seed = {}
            for iseed in [0, 1]:
                mapmuts.bayesian.Seed(iseed)
                thin = max(1, int(5e-4 * nstep))
                (pi_mean, pi_cred95, pi_traces, run_diff) = mapmuts.bayesian.InferPreferencesMCMC(copy.deepcopy(library_stats), 1.0, 1.0, 1.0, 1.0, nruns, nstep, int(burnfrac * nstep), npreburns, thin, debugging=False)
                pi_by_seed[iseed] = pi_mean
                match_diff = mapmuts.bayesian.ShannonJensenDivergence(pi_mean, pivec)
                sys.stderr.write("for seed %d, overall difference between actual and inferred pi vectors is %g (compared to ideal value of zero) and run_diff is %g (compared to ideal value of zero)... " % (iseed, match_diff, run_diff))
                sys.stderr.flush()
            seed_diff = mapmuts.bayesian.ShannonJensenDivergence(pi_by_seed[0], pi_by_seed[1])
            sys.stderr.write("the overall difference in the inferred values between the two seeds is %g (compared to ideal value of zero)." % seed_diff)
        if match_diff > match_criterion or seed_diff > match_criterion:
            msg = 'FAILED: the equilibrium preferences inference did not converge even after %d steps. The difference (summed elements of absolute value of difference vector) from actual values is %g, and the difference between replicates is %g.' % (nstep, match_diff, seed_diff)
            self.assertTrue(False, msg)
            sys.stderr.write("\n%s" % msg)
        else:
            sys.stderr.write("\nEquilibrium preference inference converged to specified tolerance of %g." % match_criterion)
        sys.stderr.write("\nEquilibrium preference inference testing complete.\n")



if __name__ == '__main__':
    if not mapmuts.bayesian.PymcAvailable():
        sys.stderr.write("\nFAILED: Cannot test equilibrium preference inference without numpy / pymc.\n")
    else:
        runner = unittest.TextTestRunner()
        unittest.main(testRunner=runner)
