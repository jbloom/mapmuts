#!python

"""Infers equilibrium preferences for each amino acid at each site.

Written by Jesse Bloom, 2013.
"""


import re
import sys
import os
import time
import tempfile
import math
import copy
import traceback
import multiprocessing
import warnings
import cPickle
import mapmuts
import mapmuts.io
import mapmuts.sequtils
import mapmuts.bayesian
import mapmuts.plot



def Entropy(pi_mean):
    """Computes site entropy in bits from array of probabilities."""
    h = 0.0 # calculate entropy
    for pi in pi_mean:
        if pi == 0:
            pass
        elif pi < 0:
            raise ValueError("Negative pi value of %g" % pi)
        else:
            h -= pi * math.log(pi, 2)
    return h


def RemoveOutlierRunMCMC(ires, library_stats, pi_concentration, epsilon_concentration, mu_concentration, rho_concentration, nruns, nsteps, burn, npreburns, thin, convergence, stepincrease, pickleresults, minvalue, seed):
    """Like *RunMCMC* but removes outlier library.

    This function can only be called if there is data for at least three
    libraries (if *len(library_stats) >= 3*). Otherwise you will get
    an exception.

    The calling and eventual result (creation of output pickle file) is
    exactly the same for this function as for *RunMCMC*. 

    However, this function analyzes the case where there is data for
    at least three libraries, identifies the outlier library for this
    site (in terms of the inferred preferences), and then does the inference
    on the *len(library_stats) - 1* libraries that do not include this
    outlier, and finally returns the inferred preferences for those
    non-outlier libraries.

    To identify the outlier library, the following algorithm is applied:

        1) For each of the *len(library_stats)* libraries, the preferences
           are inferred for that library individually.

        2) For each library, we compute the average pairwise Shannon-Jensen
           divergence between the preferences for that library and all of the other
           libraries.

        3) We remove as the outlier the library with the largest
           average pairwise Shannon-Jensen divergence.
    """
    mapmuts.bayesian.Seed(seed)
    nlibs = len(library_stats)
    if nlibs < 3:
        raise ValueError("Cannot remove outlier for less than 3 libraries.")
    pi_means = {}
    logstring = ['\nBeginning inferences for site %d.\nFirst trying to identify the library by performing the inference individually for each of the %d libraries.' % (ires, nlibs)]
    for ilib in range(nlibs):
        try:
            (fd, ilib_pickleresults) = tempfile.mkstemp()
            os.close(fd)
            RunMCMC(ires, [library_stats[ilib]], pi_concentration, epsilon_concentration, mu_concentration, rho_concentration, nruns, nsteps, burn, npreburns, thin, convergence, stepincrease, ilib_pickleresults, minvalue)
            (pi_mean, pi_cred95, pi_samples, run_diff, ilib_logstring) = cPickle.load(open(ilib_pickleresults))
            pi_means[ilib] = pi_mean
            logstring.append('\nPerforming inference for individual library %d of %d.\n%s' % (ilib + 1, nlibs, ilib_logstring.strip()))
        finally:
            if os.path.isfile(ilib_pickleresults):
                os.remove(ilib_pickleresults)
    logstring.append('\nComparing the average pairwise difference of each library to all other libraries.')
    maxdiff = maxdifflib = None
    for ilib in range(nlibs):
        diffs = []
        for jlib in range(nlibs):
            if jlib == ilib:
                continue
            diff = mapmuts.bayesian.ShannonJensenDivergence(pi_means[ilib], pi_means[jlib])
            diffs.append(diff)
        assert len(diffs) == nlibs - 1
        diff = sum(diffs) / float(nlibs - 1)    
        if maxdiff == None or diff > maxdiff:
            maxdiff = diff
            maxdifflib = ilib
        logstring.append('.. the average difference of library %d to the other libraries is %g.' % (ilib + 1, diff))
    logstring.append('\nThe library with the maximum difference is %d, so classifying this as the outlier and performing the inference using the other %d libraries.' % (maxdifflib + 1, nlibs - 1))
    library_stats = [library_stats[ilib] for ilib in range(nlibs) if ilib != maxdifflib]
    assert len(library_stats) == nlibs - 1
    try:
        (fd, combined_pickleresults) = tempfile.mkstemp()
        os.close(fd)
        RunMCMC(ires, library_stats, pi_concentration, epsilon_concentration, mu_concentration, rho_concentration, nruns, nsteps, burn, npreburns, thin, convergence, stepincrease, combined_pickleresults, minvalue)
        (pi_mean, pi_cred95, pi_samples, run_diff, combined_logstring) = cPickle.load(open(combined_pickleresults))
        logstring.append(combined_logstring)
    finally:
        if os.path.isfile(combined_pickleresults):
            os.remove(combined_pickleresults)
    logstring = ''.join(logstring)
    cPickle.dump((pi_mean, pi_cred95, pi_samples, run_diff, logstring), open(pickleresults, 'w'))
    time.sleep(1)



def RunMCMC(ires, library_stats, pi_concentration, epsilon_concentration, mu_concentration, rho_concentration, nruns, nsteps, burn, npreburns, thin, convergence, stepincrease, pickleresults, minvalue, seed):
    """Runs MCMC to infer equilibrium preferences.

    Calling variables have the same meaning as in the *main* calling function.

    This function immediately creates a copy of *library_stats*
    for subsequent use, so it is not affected if that dictionary is
    later changed.

    The result tuple is: *(pi_mean, pi_cred95, pi_samples, run_diff, logstring)*
    where the first four have the same meaning as the return variables
    from *mapmuts.bayesian.InferPreferencesMCMC* and the last is a string
    to be written to the log file. However, these results are **not** returned
    by this method. Instead, they are written using *cPickle* to the file
    specified by *pickleresults*.
    """
    mapmuts.bayesian.Seed(seed)
    library_stats = copy.deepcopy(library_stats)
    logstring = ['\nPerforming inference for site %d...' % ires]
    start_t = time.clock()
    (pi_mean, pi_cred95, pi_samples, run_diff) = \
            mapmuts.bayesian.InferPreferencesMCMC(library_stats,\
            pi_concentration, epsilon_concentration, mu_concentration,\
            rho_concentration, nruns, nsteps, burn, npreburns, thin,
            minvalue=minvalue)
    t = time.clock() - start_t
    logstring.append(" completed MCMC of %d steps in %.1f seconds" % (nsteps, t))
    if nruns > 1 and run_diff > convergence:
        logstring.append('; inference FAILED to converge (run difference of %.3g).\n' % run_diff)
        if stepincrease > 1:
            start_t = time.clock()
            logstring.append('Trying again with %d-fold more steps...' % stepincrease)
            (pi_mean, pi_cred95, pi_samples, run_diff) = \
                    mapmuts.bayesian.InferPreferencesMCMC(library_stats,\
                    pi_concentration, epsilon_concentration, mu_concentration,\
                    rho_concentration, nruns, stepincrease * nsteps,\
                    stepincrease * burn, npreburns, thin,\
                    minvalue=minvalue)
            t = time.clock() - start_t
            if run_diff <= convergence:
                logstring.append(' this time MCMC converged in %.1f seconds (run difference of %.3g).\n' % (t, run_diff))
            else:
                logstring.append(' MCMC still FAILED to converge in %.1f seconds (run difference of %.2g).\n' % (t, run_diff))
    elif nruns > 1: 
        logstring.append("; inference converged (run difference of %.3g).\n" % run_diff)
    else:
        logstring.append('.\n')
    logstring = ''.join(logstring)
    cPickle.dump((pi_mean, pi_cred95, pi_samples, run_diff, logstring), open(pickleresults, 'w'))
    time.sleep(1)


def main():
    """Main body of script."""
    # hard-coded variables
    includestop = True # include stop codons as a possible amino acid
    burnfrac = 0.1 # set burn-in to this times nsteps
    npreburns = 2 # perform this many pre-burn runs
    # check on module availability
    if not mapmuts.bayesian.PymcAvailable():
        raise ImportError("Cannot run this script as pymc or numpy are not available.")
    # read input variables
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument"\
                + ' specifying the name of the input file.')
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = mapmuts.io.ParseInfile(open(infilename))
    outfileprefix = mapmuts.io.ParseStringValue(d, 'outfileprefix')
    logfile = "%s_inferpreferences_log.txt" % outfileprefix
    log = open(logfile, 'w')
    equilibriumpreferencesfile = open('%s_equilibriumpreferences.txt' % outfileprefix, 'w')
    credibleintervalsfile = open('%s_equilibriumpreferences_credibleintervals_95.txt' % outfileprefix, 'w')
    equilibriumpreferencesfile.write('#SITE\tWT_AA\tSITE_ENTROPY')
    credibleintervalsfile.write('#SITE')
    aas = mapmuts.sequtils.AminoAcids(includestop=includestop)
    for aa in aas:
        equilibriumpreferencesfile.write('\tPI_%s' % aa)
        credibleintervalsfile.write('\tPI_%s_95cred' % aa)
    equilibriumpreferencesfile.write('\n')
    credibleintervalsfile.write('\n')
    try:
        log.write("Beginning execution of mapmuts_inferpreferences.py"\
                " in directory %s" % (os.getcwd()))
        mapmuts.io.PrintVersions(log)
        log.write("Input data being read from infile %s\n\n" % infilename)
        log.write("Progress being logged to this file, %s\n\n" % logfile)
        log.write("Read the following key/value pairs from infile %s:"\
                % (infilename))
        for (key, value) in d.iteritems():
            log.write("\n%s %s" % (key, value))
        dnafiles = mapmuts.io.ParseFileList(d, 'DNA_files')
        rnafiles = mapmuts.io.ParseFileList(d, 'RNA_files')
        mutdnafiles = mapmuts.io.ParseFileList(d, 'mutDNA_files')
        mutvirusfiles = mapmuts.io.ParseFileList(d, 'mutvirus_files')
        if not (len(dnafiles) == len(rnafiles) == len(mutdnafiles) == len(mutvirusfiles) >= 1):
            raise IOError("Failed to find four file lists (DNA, RNA, mutDNA, mutvirus) all of the same length with at least one file each.")
        epsilon_concentration = mapmuts.io.ParseFloatValue(d, 'epsilon_concentration')
        assert epsilon_concentration > 0, "epsilon_concentration must be > 0"
        mu_concentration = mapmuts.io.ParseFloatValue(d, 'mu_concentration')
        assert mu_concentration > 0, "mu_concentration must be > 0"
        rho_concentration = mapmuts.io.ParseFloatValue(d, 'rho_concentration')
        assert rho_concentration > 0, "rho_concentration must be > 0"
        pi_concentration = mapmuts.io.ParseFloatValue(d, 'pi_concentration')
        assert pi_concentration > 0, "pi_concentration must be > 0"
        minvalue = mapmuts.io.ParseFloatValue(d, 'minvalue')
        assert minvalue > 0, "minvalue must be > 0"
        seed = mapmuts.io.ParseIntValue(d, 'seed')
        nruns = mapmuts.io.ParseIntValue(d, 'nruns')
        if 'sites' in d:
            sites = mapmuts.io.ParseStringValue(d, 'sites')
            if sites.upper() in ['ALL', 'NONE']:
                sites = 'all'
            else:
                try:
                    (start, end) = sites.split()
                    (start, end) = (int(start), int(end))
                    sites = [r for r in range(start, end + 1)]
                except:
                    raise ValueError("problem parsing line for sites")
        else:
            sites = 'all'
        assert nruns >= 1, "nruns must be >= 1"
        if nruns < 2:
            warnings.warn('Will not be able to check for convergence since nruns < 2. You are suggested to use nruns >= 2.')
        nsteps = mapmuts.io.ParseIntValue(d, 'nsteps')
        burn = int(burnfrac * nsteps)
        assert nsteps >= 1 and burn >= 1, "nsteps must be set to a larger value than %d" % nsteps
        thin = mapmuts.io.ParseIntValue(d, 'thin')
        assert thin >= 1, "thin must be >= 1"
        convergence = mapmuts.io.ParseFloatValue(d, 'convergence')
        assert convergence > 0, "convergence must be > 0"
        stepincrease = mapmuts.io.ParseIntValue(d, 'stepincrease')
        assert stepincrease >= 1, "stepincrease must be >= 1"
        MCMC_traces = mapmuts.io.ParseStringValue(d, 'MCMC_traces')
        if MCMC_traces in ['None', 'False']:
            MCMC_traces = None
        elif not mapmuts.plot.PylabAvailable():
            log.write("\nWARNING: cannot create posterior plots as pylab / matplotlib are not available.\n")
            MCMC_traces = None
        elif not os.path.isdir(MCMC_traces):
            raise IOError("MCMC_traces directory of %s does not already exist. You must create this directory before running this script." % MCMC_traces)
        preference_plots = mapmuts.io.ParseStringValue(d, 'preference_plots')
        if preference_plots in ['None', 'False']:
            preference_plots = None
        elif not mapmuts.plot.PylabAvailable():
            log.write("\nWARNING: cannot create preference plots as pylab / matplotlib are not available.\n")
            preference_plots = None
        elif not os.path.isdir(preference_plots):
            raise IOError("preference_plots directory of %s does not already exist. You must create it before running this script." % preference_plots)
        log.write('\n\n')
        ncpus = mapmuts.io.ParseIntValue(d, 'ncpus')
        assert ncpus >= 1, "ncpus must be >= 1"
        removeoutlier = mapmuts.io.ParseBoolValue(d, 'removeoutlier')
        if removeoutlier and len(dnafiles) < 3:
            raise ValueError("You cannot set removeoutlier to True unless you have specified at least three libraries for analysis.")

        # read codon counts and set up priors
        log.write("Reading in the codon counts...")
        log.flush()
        nlibs = len(dnafiles)
        dna_libs = [mapmuts.io.ReadCodonCounts(open(f)) for f in dnafiles]
        rna_libs = [mapmuts.io.ReadCodonCounts(open(f)) for f in rnafiles]
        mutdna_libs = [mapmuts.io.ReadCodonCounts(open(f)) for f in mutdnafiles]
        mutvirus_libs = [mapmuts.io.ReadCodonCounts(open(f)) for f in mutvirusfiles]
        for x in dna_libs + rna_libs + mutdna_libs + mutvirus_libs:
            mapmuts.sequtils.ClassifyCodonCounts(x)
        log.write(" completed reading the codon counts.\n")
        assert nlibs == len(dna_libs) == len(rna_libs) == len(mutdna_libs) == len(mutvirus_libs)
        mu_priors = []
        rho_priors = []
        epsilon_priors = []
        for ilib in range(nlibs):
            mu_priors.append((mutdna_libs[ilib]['TOTAL_MUT'] / float(mutdna_libs[ilib]['TOTAL_COUNTS']) - dna_libs[ilib]['TOTAL_MUT'] / float(dna_libs[ilib]['TOTAL_COUNTS'])) / 63)
            rho = {}
            epsilon = {}
            for (ndiffs, denom) in [(1, 9.0), (2, 27.0), (3, 27.0)]:
                epsilon[ndiffs] = dna_libs[ilib]['TOTAL_N_%dMUT' % ndiffs] / float(dna_libs[ilib]['TOTAL_COUNTS']) / denom
                rho[ndiffs] = rna_libs[ilib]['TOTAL_N_%dMUT' % ndiffs] / float(rna_libs[ilib]['TOTAL_COUNTS']) / denom - epsilon[ndiffs]
            epsilon_priors.append(epsilon)
            rho_priors.append(rho)
        # get the protein length
        integer_keys = [key for key in dna_libs[0].keys() if isinstance(key, int)]
        protlength = max(integer_keys)
        assert protlength == len(integer_keys), "Problem with residue numbering?"
        if sites == 'all':
            sites = [r for r in range(1, protlength + 1)]
        # get list of all codons
        codons = mapmuts.sequtils.Codons()
        # begin inference
        log.write('\nNow beginning inference of the equilibrium preferences...\n')
        log.flush()
        # set up processes for all residues
        processes = {} # keyed by residue number, value is multiprocessing.Process
        wtaa_d = {} # keyed by residue number, value is wtaa
        pickleresults = {} # keyed by residue number, value is pickle file name
        for ires in sites:
            library_stats = []
            wtaas = []
            for ilib in range(nlibs):
                wtcodon = dna_libs[ilib][ires]['WT']
                assert wtcodon == rna_libs[ilib][ires]['WT'] == mutdna_libs[ilib][ires]['WT'] == mutvirus_libs[ilib][ires]['WT'], "Mismatches of wildtype codon for residue %d" % ires
                wtaa = mapmuts.sequtils.Translate([('wt', wtcodon)])[0][1]
                if not wtaa:
                    wtaa = '*'
                wtaas.append(wtaa)
                ilib_stats = {'wtcodon':wtcodon,
                              'mu_prior':mu_priors[ilib],
                              'epsilon_prior':epsilon_priors[ilib],
                              'rho_prior':rho_priors[ilib],
                              'nrdna_counts':dict([(codon, dna_libs[ilib][ires][codon]) for codon in codons]),
                              'nrrna_counts':dict([(codon, rna_libs[ilib][ires][codon]) for codon in codons]),
                              'nrmutdna_counts':dict([(codon, mutdna_libs[ilib][ires][codon]) for codon in codons]),
                              'nrmutvirus_counts':dict([(codon, mutvirus_libs[ilib][ires][codon]) for codon in codons]),
                             }
                assert dna_libs[ilib][ires]['PAIRED_COUNTS'] == sum(ilib_stats['nrdna_counts'].values()), "Didn't get the right number of DNA counts."
                assert rna_libs[ilib][ires]['PAIRED_COUNTS'] == sum(ilib_stats['nrrna_counts'].values()), "Didn't get the right number of RNA counts."
                assert mutdna_libs[ilib][ires]['PAIRED_COUNTS'] == sum(ilib_stats['nrmutdna_counts'].values()), "Didn't get the right number of mutDNA counts."
                assert mutvirus_libs[ilib][ires]['PAIRED_COUNTS'] == sum(ilib_stats['nrmutvirus_counts'].values()), "Didn't get the right number of mutvirus counts."
                library_stats.append(ilib_stats)
            # determine the wildtype aa - either all the same, or list of ones for each lib
            assert nlibs == len(wtaas) > 0
            if wtaas.count(wtaas[0]) == nlibs:
                wtaa = wtaas[0]
            else:
                wtaa = ','.join(wtaas)
            wtaa_d[ires] = wtaa
            (fd, pickleresults[ires]) = tempfile.mkstemp()
            os.close(fd)
            if removeoutlier:
                processes[ires] = multiprocessing.Process(target=RemoveOutlierRunMCMC, args=(ires, copy.deepcopy(library_stats), pi_concentration, epsilon_concentration, mu_concentration, rho_concentration, nruns, nsteps, burn, npreburns, thin, convergence, stepincrease, pickleresults[ires], minvalue, seed))
            else:
                processes[ires] = multiprocessing.Process(target=RunMCMC, args=(ires, copy.deepcopy(library_stats), pi_concentration, epsilon_concentration, mu_concentration, rho_concentration, nruns, nsteps, burn, npreburns, thin, convergence, stepincrease, pickleresults[ires], minvalue, seed))
        # now start running the processes in order of residues. Don't start the next
        # until the first residue still running is done.
        processes_running = dict([(ires, False) for ires in processes.iterkeys()])
        processes_finished = dict([(ires, False) for ires in processes.iterkeys()])
        for ires in sites:
            i = 0
            while (ires + i  <= sites[-1]) and (processes_running.values().count(True) < ncpus):
                if (not processes_finished[ires + i]) and (not processes_running[ires + i]):
                    processes[ires + i].start()
                    processes_running[ires + i] = True
                    time.sleep(1)
                i += 1
            if not processes_running[ires]:
                raise ValueError("Process for ires %d should be running" % ires)

            while processes[ires].is_alive():
                time.sleep(1)
            if processes[ires].exitcode:
                raise ValueError("Error running MCMC for residue %d" % ires)
            processes_finished[ires] = True
            processes_running[ires] = False
            (pi_mean, pi_cred95, pi_samples, run_diff, logstring) = cPickle.load(open(pickleresults[ires]))
            os.remove(pickleresults[ires])
            log.write(logstring)
            log.flush()
            # make plots
            h = Entropy(pi_mean)
            if MCMC_traces:
                traces = []
                trace_labels = []
                for irun in range(nruns):
                    for iaa in range(len(aas)):
                        traces.append(pi_samples[irun].transpose()[iaa])
                        trace_labels.append("%s" % aas[iaa])
                mapmuts.plot.PlotTraces(traces, "%s/%s_residue_%d.pdf" % (MCMC_traces, outfileprefix, ires), xlabel='MCMC step', ylabel='equilibrium preference', title='residue %d MCMC traces' % ires, trace_labels=trace_labels)
            if preference_plots:
                pis = dict(zip(aas, pi_mean))
                pi_errs = dict(zip(aas, pi_cred95))
                mapmuts.plot.PlotEquilibriumFreqs(pis, "%s/%s_residue_%d.pdf" % (preference_plots, outfileprefix, ires), title='residue %d, site entropy of %.2f bits' % (ires, h), pi_errs=pi_errs)
            # now write pi values to file
            assert abs(sum(pi_mean) - 1.0) < 1e-7, "Sum of pi values of %g not close to one." % sum(pi_mean)
            equilibriumpreferencesfile.write('%d\t%s\t%g' % (ires, wtaa_d[ires], h))
            credibleintervalsfile.write('%d' % ires)
            assert len(aas) == len(pi_mean), "Didn't get the right number of pi values"
            for i in range(len(aas)):
                equilibriumpreferencesfile.write('\t%g' % pi_mean[i])
                credibleintervalsfile.write('\t%g,%g' % tuple(pi_cred95[i]))
            equilibriumpreferencesfile.write('\n')
            equilibriumpreferencesfile.flush()
            credibleintervalsfile.write('\n')
            credibleintervalsfile.flush()
    except:
        (exc_type, exc_value, exc_traceback) = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_traceback, file=log)
        log.write("\n\nPrematurely closing log due to execution error.")
        raise
    finally:
        log.write("\n\nExecution completed at %s." % time.ctime())
        log.close()
        equilibriumpreferencesfile.close()
        credibleintervalsfile.close()
        try:
            for (i, is_running) in processes_running.iteritems():
                if is_running and processes[i].is_alive():
                    processes[i].terminate()
        except NameError:
            pass
        try:
            for f in pickleresults.itervalues():
                if os.path.isfile(f):
                    os.remove(f)
        except NameError:
            pass



if __name__ == '__main__':
    main() # run the script
