#!python

"""Infers differential preferences for each amino acid at each site.

Written by Jesse Bloom, 2014.
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


def RMS(dpi_mean):
    """Computes root mean square value of entries."""
    rms = 0.0
    for x in dpi_mean:
        rms += x**2
    rms = math.sqrt(rms)
    return rms


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


def RunMCMC(ires, error_control_counts, starting_sample_counts, control_selection_counts, selection_counts, wtcodon, f_prior, epsilon_prior, pi_concentration, epsilon_concentration, f_concentration, deltapi_concentration, nruns, nsteps, burn, thin, minvalue, convergence, stepincrease, pickleresults, seed):
    """Runs MCMC to infer differential preferences.

    Calling variables have same meaning as in the *main* function.

    The results are written using *cPickle* to the file
    specified by *pickleresults*.
    """
    mapmuts.bayesian.Seed(seed)
    logstring = ['\nPerforming inference for site %d...' % ires]
    start_t = time.clock()
    returnvalue = \
            mapmuts.bayesian.InferDifferentialPreferencesMCMC(\
            error_control_counts, starting_sample_counts, control_selection_counts, selection_counts, wtcodon, f_prior, epsilon_prior, pi_concentration, epsilon_concentration, f_concentration, deltapi_concentration, nruns, nsteps, burn, thin, minvalue=minvalue)
    t = time.clock() - start_t
    run_diffs = [(selection, returnvalue[selection][3]) for selection in returnvalue.iterkeys()]
    logstring.append(" completed MCMC of %d steps in %.1f seconds" % (nsteps, t))
    if nruns > 1 and max([tup[1] for tup in run_diffs]) > convergence:
        logstring.append('; inference FAILED to converge (run differences of: %s).\n' % ', '.join(['%g for %s' % (tup[1], tup[0]) for tup in run_diffs]))
        if stepincrease > 1:
            start_t = time.clock()
            logstring.append('Trying again with %d-fold more steps...' % stepincrease)
            returnvalue = \
                mapmuts.bayesian.InferDifferentialPreferencesMCMC(\
                error_control_counts, starting_sample_counts, control_selection_counts, selection_counts, wtcodon, f_prior, epsilon_prior, pi_concentration, epsilon_concentration, f_concentration, deltapi_concentration, nruns, int(stepincrease * nsteps), burn, thin, minvalue=minvalue)
            assert len(returnvalue) >= 2, "Should be at least the control preferences and one differential preference"
            t = time.clock() - start_t
            run_diffs = [(selection, returnvalue[selection][3]) for selection in returnvalue.iterkeys()]
            if max([tup[1] for tup in run_diffs]) <= convergence:
                logstring.append(' this time MCMC converged in %.1f seconds (run differences of: %s).\n' % (t, ', '.join(['%g for %s' % (tup[1], tup[0]) for tup in run_diffs])))
            else:
                logstring.append(' MCMC still FAILED to converge in %.1f seconds (run differences of: %s).\n' % (t, ', '.join(['%g for %s' % (tup[1], tup[0]) for tup in run_diffs])))
    elif nruns > 1: 
        logstring.append("; inference converged (run differences of: %s).\n" % ', '.join(['%g for %s' % (tup[1], tup[0]) for tup in run_diffs]))
    else:
        logstring.append('.\n')
    logstring = ''.join(logstring)
    cPickle.dump((logstring, returnvalue), open(pickleresults, 'w'))
    time.sleep(1)


def main():
    """Main body of script."""
    # hard-coded variables
    includestop = True # include stop codons as a possible amino acid
    burnfrac = 0.2 # set burn-in to this times nsteps
    # check on module availability
    if not mapmuts.bayesian.PymcAvailable():
        raise ImportError("Cannot run this script as pymc or numpy are not available.")
    aas = mapmuts.sequtils.AminoAcids(includestop=includestop)
    codons = mapmuts.sequtils.Codons()

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
    if outfileprefix.upper() == 'NONE':
        outfileprefix = ''
    logfile = "%sinferdifferentialpreferences_log.txt" % outfileprefix
    log = open(logfile, 'w')
    log.write("Beginning execution of mapmuts_inferpreferences.py"\
            " in directory %s" % (os.getcwd()))
    mapmuts.io.PrintVersions(log)
    log.write("Input data being read from infile %s\n\n" % infilename)
    log.write("Progress being logged to this file, %s\n\n" % logfile)
    log.write("Read the following key/value pairs from infile %s:"\
            % (infilename))
    for (key, value) in d.iteritems():
        log.write("\n%s %s" % (key, value))
    codoncounts_data = {} # dictionary keyed by sample type
    for sample in ['error_control', 'starting_sample', 'control_selection']:
        fname = mapmuts.io.ParseStringValue(d, sample)
        if not os.path.isfile(fname):
            raise IOError("Failed to find file %s specified by %s" % (fname, sample))
        codoncounts_data[sample] = mapmuts.io.ReadCodonCounts(open(fname))
    selections = []
    for (key, value) in d.iteritems():
        m = re.search('^selection_(?P<sample>\S+)$', key)
        if m:
            sample = m.group('sample')
            if sample in codoncounts_data:
                raise ValueError("Duplicate selection sample of name %s" % sample)
            if sample in ['error_control', 'starting_selection', 'control_selection']:
                raise ValueError("Selection sample cannot have name %s" % sample)
            if not os.path.isfile(value):
                raise IOError("Failed to find file %s specified for selection sample %s" % (value, sample))
            codoncounts_data[sample] = mapmuts.io.ReadCodonCounts(open(value))
            selections.append(sample)
    if not selections:
        raise ValueError("Failed to find any selected samples with keys of form selection_???")
    for value in codoncounts_data.itervalues():
        mapmuts.sequtils.ClassifyCodonCounts(value)
    log.write("\n\nRead codon counts for the following samples: %s" % ', '.join(codoncounts_data.keys()))
    epsilon_concentration = mapmuts.io.ParseFloatValue(d, 'epsilon_concentration')
    assert epsilon_concentration > 0, "epsilon_concentration must be > 0"
    f_concentration = mapmuts.io.ParseFloatValue(d, 'f_concentration')
    assert f_concentration > 0, "f_concentration must be > 0"
    deltapi_concentration = mapmuts.io.ParseFloatValue(d, 'deltapi_concentration')
    assert deltapi_concentration > 0, "deltapi_concentration must be > 0"
    pi_concentration = mapmuts.io.ParseFloatValue(d, 'pi_concentration')
    assert pi_concentration > 0, "pi_concentration must be > 0"
    minvalue = mapmuts.io.ParseFloatValue(d, 'minvalue')
    assert minvalue > 0, "minvalue must be > 0"
    seed = mapmuts.io.ParseIntValue(d, 'seed')
    nruns = mapmuts.io.ParseIntValue(d, 'nruns')
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
    ncpus = mapmuts.io.ParseIntValue(d, 'ncpus')
    assert ncpus >= 1, "ncpus must be >= 1"
    MCMC_traces = mapmuts.io.ParseStringValue(d, 'MCMC_traces')
    if MCMC_traces in ['None', 'False']:
        MCMC_traces = None
    elif not mapmuts.plot.PylabAvailable():
        log.write("\nWARNING: cannot create posterior plots as pylab / matplotlib are not available.\n")
        MCMC_traces = None
    elif not os.path.isdir(MCMC_traces):
        os.mkdir(MCMC_traces)
    preference_plots = mapmuts.io.ParseStringValue(d, 'preference_plots')
    if preference_plots in ['None', 'False']:
        preference_plots = None
    elif not mapmuts.plot.PylabAvailable():
        log.write("\nWARNING: cannot create preference plots as pylab / matplotlib are not available.\n")
        preference_plots = None
    elif not os.path.isdir(preference_plots):
        os.mkdir(preference_plots)

    # Now set up to run the MCMC
    # first, compute the parameters needed for the priors
    starting_mutrate = codoncounts_data['starting_sample']['TOTAL_MUT'] / float(codoncounts_data['starting_sample']['TOTAL_COUNTS'])
    error_rate = codoncounts_data['error_control']['TOTAL_MUT'] / float(codoncounts_data['error_control']['TOTAL_COUNTS'])
    f_prior = (starting_mutrate - error_rate) / float(len(codons) - 1)
    log.write('\nThe prior estimate for the frequency of any specific mutation in starting_sample is %g (overall mutation rate of %g in starting_sample minus overall error rate of %g for error_control, divided by number of codons).' % (f_prior, starting_mutrate, error_rate))
    epsilon_prior = {}
    for (ndiffs, denom) in [(1, 9.0), (2, 27.0), (3, 27.0)]:
        epsilon_prior[ndiffs] = codoncounts_data['error_control']['TOTAL_N_%dMUT' % ndiffs] / float(codoncounts_data['error_control']['TOTAL_COUNTS']) / denom
        log.write('\nThe prior estimate for the error rate in error_control for a mutation with %d nucleotide changes is %g.' % (ndiffs, epsilon_prior[ndiffs]))

    # now get a list of all (sites, wildtype_codon)
    sites = [(key, codoncounts_data['error_control'][key]['WT']) for key in codoncounts_data['error_control'].keys() if isinstance(key, int)]
    sites.sort()
    for (sample, sampledata) in codoncounts_data.iteritems():
        samplesites = [(key, sampledata[key]['WT']) for key in sampledata.keys() if isinstance(key, int)]
        samplesites.sort()
        if sites != samplesites:
            raise ValueError("Not all samples specify the same set of sites / wildtype codons.")
    log.write("\nData is specified for %d sites.\n" % len(sites))
    preferencesfile = '%spreferences_control_selection.txt' % outfileprefix
    preferencescred95file = '%spreferences_control_selection_credibleintervals_95.txt' % outfileprefix
    log.write('\nPreferences for control selection will be written to %s and %s.\n' % (preferencesfile, preferencescred95file))
    preferencesfile = open(preferencesfile, 'w')
    preferencesfile.write('#SITE\tWT_AA\tSITE_ENTROPY\t%s\n' % '\t'.join(['PI_%s' % aa for aa in aas]))
    preferencescred95file = open(preferencescred95file, 'w')
    preferencescred95file.write('SITE\t%s\n' % '\t'.join(['PI_%s_95cred' % aa for aa in aas]))
    meanfiles = {}
    cred95files = {}
    for selection in selections:
        meanfiles[selection] = '%sdifferentialpreferences_selection_%s.txt' % (outfileprefix, selection)
        cred95files[selection] = '%sdifferentialpreferences_selection_%s_credibleintervals_95.txt' % (outfileprefix, selection)
        log.write('\nDifferential preferences for selection %s will be written to %s and %s.\n' % (selection, meanfiles[selection], cred95files[selection]))
        meanfiles[selection] = open(meanfiles[selection], 'w')
        cred95files[selection] = open(cred95files[selection], 'w')
        meanfiles[selection].write('#SITE\tWT_AA\tRMS_dPI\t%s\n' % '\t'.join(['dPI_%s' % aa for aa in aas]))
        cred95files[selection].write('#SITE\t%s\n' % '\t'.join(['dPI_%s_95cred' % aa for aa in aas]))
    log.write('\nNow beginning inferences...\n')
    log.flush()
    processes = {} # keyed by residue number, value is multiprocessing.Process
    wtaa_d = {} # keyed by residue number, value is wtaa
    pickleresults = {} # keyed by residue number, value is pickle file name
    try:
        # set up the processes
        for (ires, wtcodon) in sites:
            wtaa = mapmuts.sequtils.Translate([('wt', wtcodon)])[0][1]
            if not wtaa:
                wtaa = '*'
            wtaa_d[ires] = wtaa
            (fd, pickleresults[ires]) = tempfile.mkstemp()
            os.close(fd)
            error_control_counts = dict([(codon, codoncounts_data['error_control'][ires][codon]) for codon in codons])
            starting_sample_counts = dict([(codon, codoncounts_data['starting_sample'][ires][codon]) for codon in codons])
            control_selection_counts = dict([(codon, codoncounts_data['control_selection'][ires][codon]) for codon in codons])
            selection_counts = {}
            for selection in selections:
                selection_counts[selection] = dict([(codon, codoncounts_data[selection][ires][codon]) for codon in codons])
            processes[ires] = multiprocessing.Process(target=RunMCMC, args=(ires, error_control_counts, starting_sample_counts, control_selection_counts, selection_counts, wtcodon, f_prior, epsilon_prior, pi_concentration, epsilon_concentration, f_concentration, deltapi_concentration, nruns, nsteps, burn, thin, minvalue, convergence, stepincrease, pickleresults[ires], seed))
        # start running processes. Don't start the next
        # until the first residue still running is done.
        processes_running = dict([(ires, False) for ires in processes.iterkeys()])
        processes_finished = dict([(ires, False) for ires in processes.iterkeys()])
        for (ires, wtcodon) in sites:
            i = 0
            while (ires + i  <= sites[-1][0]) and (processes_running.values().count(True) < ncpus):
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
            (logstring, returnvalue) = cPickle.load(open(pickleresults[ires]))
            os.remove(pickleresults[ires])
            log.write(logstring)
            log.flush()
            (mean, cred95, traces, run_diff) = returnvalue['control_selection']
            assert len(aas) == len(mean) == len(cred95), "Not right number of entries"
            assert abs(sum(mean) - 1.0) < 1e-7, "Sum of control preferences of %g not close to one." % sum(mean)
            preferencesfile.write('%d\t%s\t%g\t%s\n' % (ires, wtaa_d[ires], Entropy(mean), '\t'.join(['%g' % pi for pi in mean])))
            preferencescred95file.write('%d\t%s\n' % (ires, '\t'.join(['%g,%g' % (x[0], x[1]) for x in cred95])))
            preferencesfile.flush()
            preferencescred95file.flush()
            for selection in selections:
                (mean, cred95, traces, run_diff) = returnvalue[selection]
                assert len(aas) == len(mean) == len(cred95), "Not right number of entries"
                assert abs(sum(mean)) < 1e-7, "Sum of differential preferences of %g not close to one." % sum(mean)
                meanfiles[selection].write('%d\t%s\t%g\t%s\n' % (ires, wtaa_d[ires], RMS(mean), '\t'.join(['%g' % dpi for dpi in mean])))
                cred95files[selection].write('%d\t%s\n' % (ires, '\t'.join(['%g,%g' % (x[0], x[1]) for x in cred95])))
                meanfiles[selection].flush()
                cred95files[selection].flush()
            if MCMC_traces:
                for selection in ['control_selection'] + selections:
                    plottraces = []
                    trace_labels = []
                    for irun in range(nruns):
                        plottraces += [returnvalue[selection][2][irun].transpose()[iaa] for iaa in range(len(aas))]
                        trace_labels += [aas[iaa] for iaa in range(len(aas))]
                    if selection == 'control_selection':
                        plotname = '%s/%spreferences_control_selection_%d.pdf' % (MCMC_traces, outfileprefix, ires)
                        ylabel = 'preference'
                        title = 'control selection preferences, residue %d' % ires
                    else:
                        plotname = '%s/%sdifferentialpreferences_selection_%s_%d.pdf' % (MCMC_traces, outfileprefix, selection, ires)
                        ylabel = 'differential preference'
                        title = 'selection %s differential preferences, residue %d' % (selection, ires)
                    mapmuts.plot.PlotTraces(plottraces, plotname, xlabel='MCMC step', ylabel=ylabel, title=title, trace_labels=trace_labels)
                    log.write('Wrote MCMC traces to %s\n' % plotname)
                    log.flush()
            if preference_plots:
                for selection in ['control_selection'] + selections:
                    if selection == 'control_selection':
                        plotname = '%s/%spreferences_control_selection_%d.pdf' % (preference_plots, outfileprefix, ires)
                        differentialpreferences = False
                    else:
                        plotname = '%s/%sdifferentialpreferences_selection_%s_%d.pdf' % (preference_plots, outfileprefix, selection, ires)
                        differentialpreferences = True
                    (mean, cred95) = (dict(zip(aas, returnvalue[selection][0])), dict(zip(aas, returnvalue[selection][1])))
                    mapmuts.plot.PlotEquilibriumFreqs(mean, plotname, 'residue %d' % ires, pi_errs=cred95, differentialpreferences=differentialpreferences)
    except:
        (exc_type, exc_value, exc_traceback) = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_traceback, file=log)
        log.write("\n\nPrematurely closing log due to execution error.")
        raise
    finally:
        log.write("\n\nExecution completed at %s." % time.ctime())
        log.close()
        preferencesfile.close()
        preferencescred95file.close()
        for selection in selections:
            meanfiles[selection].close()
            cred95files[selection].close()
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
