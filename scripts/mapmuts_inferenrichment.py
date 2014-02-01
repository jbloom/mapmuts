#!python

"""Infers enrichment ratios.

Written by Jesse Bloom, 2013.
"""


import re
import sys
import os
import time
import warnings
import mapmuts
import mapmuts.io
import mapmuts.sequtils
import mapmuts.bayesian
import mapmuts.plot


def main():
    """Main body of script."""
    # check on module availability
    if not mapmuts.bayesian.PymcAvailable():
        raise ImportError("Cannot run this script as pymc or numpy are not available.")
    if not mapmuts.bayesian.ScipyAvailable():
        warnings.warn("Cannot import scipy. The MCMC in this script will be less efficient. Installation of scipy is strongly recommended to improve performance.\n")
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
    logfile = "%s_inferenrichment_log.txt" % outfileprefix
    log = open(logfile, 'w')
    equilibriumfreqsfile = open('%s_equilibriumfreqs.txt' % outfileprefix, 'w')
    equilibriumfreqsfile.write('#SITE\tWT_AA\tSITE_ENTROPY')
    for aa in mapmuts.sequtils.AminoAcids():
        equilibriumfreqsfile.write('\tPI_%s' % aa)
    equilibriumfreqsfile.write('\n')
    enrichmentratiosfile = open('%s_enrichmentratios.txt' % outfileprefix, 'w')
    enrichmentratiosfile.write('#MUTATION\tPHI\tPHI_HPD95_LOW\tPHI_HPD95_HIGH\tDIRECT_RATIO\tMUTDNA_COUNTS\n')
    try:
        log.write("Beginning execution of mapmuts_inferenrichment.py"\
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
        # parse any excludesite_NNN entries
        excludekeymatch = re.compile('excludesite\_(?P<site>\d+)')
        excludekeys = [key for key in d.iterkeys() if excludekeymatch.search(key)]
        exclude_libs = {} # keyed by site, values dictionary keyed by excluded num (0, 1, 2, ..)
        for key in excludekeys:
            site = int(excludekeymatch.search(key).group('site'))
            exclude_libs[site] = {}
            entries = d[key].split()
            if len(entries) != len(dnafiles):
                raise ValueError("Invalid number of entries for %s" % key)
            for i in range(len(entries)):
                if entries[i].strip() == 'use':
                    continue
                elif entries[i].strip() == 'exclude':
                    exclude_libs[site][i] = True                    
                else:
                    raise ValueError("Entries for %s must be 'use' or 'exclude', not %s" % (key, entries[i]))
        #
        alpha = mapmuts.io.ParseFloatValue(d, 'alpha')
        phi_prior = mapmuts.io.ParseFloatValue(d, 'phi_prior')
        if not phi_prior > 0:
            raise ValueError("phi_prior must be greater than zero")
        log.write("\nUsing a prior for phi of phi_prior = %.4f\n" % phi_prior)
        if not (alpha > 0):
            raise ValueError("alpha must be greater than zero")
        minbeta = mapmuts.io.ParseFloatValue(d, 'minbeta')
        if not (minbeta > 0):
            raise ValueError("minbeta must be greater than zero")
        seed = mapmuts.io.ParseIntValue(d, 'seed')
        mapmuts.bayesian.Seed(seed)
        nruns = mapmuts.io.ParseIntValue(d, 'nruns')
        if nruns < 2:
            warnings.warn('Will not be able to check for convergence since nruns < 2')
        nsteps = mapmuts.io.ParseIntValue(d, 'nsteps')
        burn = mapmuts.io.ParseIntValue(d, 'burn')
        thin = mapmuts.io.ParseIntValue(d, 'thin')
        convergence = mapmuts.io.ParseFloatValue(d, 'convergence')
        if not (convergence > 1):
            raise ValueError("convergence must be greater than one")
        stepincrease = mapmuts.io.ParseIntValue(d, 'stepincrease')
        convergencewarning = mapmuts.io.ParseBoolValue(d, 'convergencewarning')
        MCMC_traces = mapmuts.io.ParseStringValue(d, 'MCMC_traces')
        if MCMC_traces in ['None', 'False']:
            MCMC_traces = None
        elif not mapmuts.plot.PylabAvailable():
            log.write("\nWARNING: cannot create posterior plots as pylab / matplotlib are not available.\n")
            MCMC_traces = None
        elif not os.path.isdir(MCMC_traces):
            raise IOError("MCMC_traces directory of %s does not already exist. You must create it before running this script." % MCMC_traces)
        enrichmentratio_plots = mapmuts.io.ParseStringValue(d, 'enrichmentratio_plots')
        if enrichmentratio_plots in ['None', 'False']:
            enrichmentratio_plots = None
        elif not mapmuts.plot.PylabAvailable():
            log.write("\nWARNING: cannot create enrichment ratio plots as pylab / matplotlib are not available.\n")
            enrichmentratio_plots = None
        elif not os.path.isdir(enrichmentratio_plots):
            raise IOError("enrichmentratio_plots directory of %s does not already exist. You must create it before running this script." % enrichmentratio_plots)
        equilibriumfreqs_plots = mapmuts.io.ParseStringValue(d, 'equilibriumfreqs_plots')
        if equilibriumfreqs_plots in ['None', 'False']:
            equilibriumfreqs_plots = None
        elif not mapmuts.plot.PylabAvailable():
            log.write("\nWARNING: cannot create equilibrium frequency plots as pylab / matplotlib are not available.\n")
            equilibriumfreqs_plots = None
        elif not os.path.isdir(equilibriumfreqs_plots):
            raise IOError("equilibriumfreqs_plots directory of %s does not already exist. You must create it before running this script." % equilibriumfreqs_plots)
        log.write('\n\n')
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
            mu_priors.append(max(minbeta, (mutdna_libs[ilib]['TOTAL_MUT'] / float(mutdna_libs[ilib]['TOTAL_COUNTS']) - dna_libs[ilib]['TOTAL_MUT'] / float(dna_libs[ilib]['TOTAL_COUNTS'])) / 63))
            rho = {}
            epsilon = {}
            for (ndiffs, denom) in [(1, 9.0), (2, 27.0), (3, 27.0)]:
                epsilon[ndiffs] = max(minbeta, dna_libs[ilib]['TOTAL_N_%dMUT' % ndiffs] / float(dna_libs[ilib]['TOTAL_COUNTS']) / denom)
                rho[ndiffs] = max(minbeta, rna_libs[ilib]['TOTAL_N_%dMUT' % ndiffs] / float(rna_libs[ilib]['TOTAL_COUNTS']) / denom - epsilon[ndiffs])
            epsilon_priors.append(epsilon)
            rho_priors.append(rho)
        # get the protein length
        integer_keys = [key for key in dna_libs[0].keys() if isinstance(key, int)]
        protlength = max(integer_keys)
        assert protlength == len(integer_keys), "Problem with residue numbering?"
        # loop over all residues
        for ires in range(1, protlength + 1):
            phis = {}
            wtcodon = dna_libs[0][ires]['WT']
            wtaa = mapmuts.sequtils.Translate([('wt', wtcodon)])[0][1]
            if not wtaa:
                wtaa = '*'
            decorated_list = [] 
            for (mutaa, mutcodons) in mapmuts.sequtils.MutAAsCodons(wtcodon):
                if MCMC_traces:
                    plot_phi_traces = "%s/%s_%s%d%s.pdf" % (MCMC_traces, outfileprefix, wtaa, ires, mutaa)
                else:
                    plot_phi_traces = False
                start_t = time.clock()
                log.write('\nPerforming inference for %s%d%s...' % (wtaa, ires, mutaa))
                log.flush()
                library_stats = []
                for ilib in range(nlibs):
                    if (ires in exclude_libs) and (ilib in exclude_libs[site]):
                        continue
                    assert dna_libs[ilib][ires]['WT'] == rna_libs[ilib][ires]['WT'] == mutdna_libs[ilib][ires]['WT'] == mutvirus_libs[ilib][ires]['WT'] == wtcodon
                    libstats = {'mu_prior':mu_priors[ilib],
                                'epsilon_prior':epsilon_priors[ilib],
                                'rho_prior':rho_priors[ilib],
                                'Nrdna':dna_libs[ilib][ires]['PAIRED_COUNTS'],
                                'Nrrna':rna_libs[ilib][ires]['PAIRED_COUNTS'],
                                'Nrmutdna':mutdna_libs[ilib][ires]['PAIRED_COUNTS'],
                                'Nrmutvirus':mutvirus_libs[ilib][ires]['PAIRED_COUNTS'],
                                'nrdna_list':[dna_libs[ilib][ires][codon] for codon in mutcodons],
                                'nrrna_list':[rna_libs[ilib][ires][codon] for codon in mutcodons],
                                'nrmutdna_list':[mutdna_libs[ilib][ires][codon] for codon in mutcodons],
                                'nrmutvirus_list':[mutvirus_libs[ilib][ires][codon] for codon in mutcodons],
                               }
                    library_stats.append(libstats)
                if not library_stats:
                    raise ValueError("No included libraries for %s%d" % (mutaa, ires))
                (phi_mean, phi_hpd95, phi_samples, converged) = mapmuts.bayesian.InferEnrichmentMCMC_2(alpha, phi_prior, wtcodon, mutcodons, library_stats, nruns, nsteps, burn, thin, convergence=convergence, plot_phi_traces=plot_phi_traces)
                t = time.clock() - start_t
                log.write(" completed MCMC of %d steps in %.1f seconds; inferred phi of %.4f (%.4f to %.4f)" % (nsteps, t, phi_mean, phi_hpd95[0], phi_hpd95[1]))
                if nruns > 1 and not converged: 
                    log.write('; inference FAILED to converge.\n')
                    log.flush()
                    if stepincrease > 1:
                        start_t = time.clock()
                        log.write('Trying again with %d-fold more steps...' % stepincrease)
                        log.flush()
                        (phi_mean, phi_hpd95, phi_samples, converged) = mapmuts.bayesian.InferEnrichmentMCMC_2(alpha, phi_prior, wtcodon, mutcodons, library_stats, nruns, stepincrease * nsteps, burn * stepincrease, thin, convergence=convergence, plot_phi_traces=plot_phi_traces)
                        t = time.clock() - start_t
                        if converged:
                            log.write(' this time MCMC converged in %.1f seconds; inferred phi of %.4f (%.4f to %.4f).\n' % (t, phi_mean, phi_hpd95[0], phi_hpd95[1]))
                        else:
                            log.write(' MCMC still FAILED to converge in %.1f seconds; inferred phi of %.4f (%.4f to %.4f).\n' % (t, phi_mean, phi_hpd95[0], phi_hpd95[1]))
                    if not converged:
                        warnings.warn('Inference failed to converge for %s%d%s. Using non-converged estimate.' % (wtaa, ires, mutaa), RuntimeWarning)
                elif converged:
                    log.write("; inference converged.\n")
                else:
                    log.write('.\n')
                (direct_ratio, mutdnacounts) = mapmuts.bayesian.DirectEnrichmentRatio(library_stats)
                log.write("For comparison, the direct ratio is %.4f with %d mutDNA counts.\n" % (direct_ratio, mutdnacounts))
                log.flush()
                enrichmentratiosfile.write("%s%d%s\t%f\t%f\t%f\t%f\t%d\n" % (wtaa, ires, mutaa, phi_mean, phi_hpd95[0], phi_hpd95[1], direct_ratio, mutdnacounts))
                enrichmentratiosfile.flush()
                decorated_list.append(("%s%d%s" % (wtaa, ires, mutaa), phi_mean, phi_hpd95[0], phi_hpd95[1]))
                if mutaa != '*':
                    phis[mutaa] = phi_mean
            decorated_list.sort()
            if enrichmentratio_plots:
                mutations = [x[0] for x in decorated_list]
                ratios = [x[1] for x in decorated_list]
                ratio_low = [x[2] for x in decorated_list]
                ratio_high = [x[3] for x in decorated_list]
                plotfile = "%s/%s_%s%d.pdf" % (enrichmentratio_plots, outfileprefix, wtaa, ires)
                mapmuts.plot.PlotEnrichmentRatios(mutations, ratios, [ratio_low, ratio_high], plotfile)
            pis = mapmuts.bayesian.EquilibriumFracs(wtaa, phis) # equilibrium fracs
            assert len(pis) == 20 and abs(sum(pis.values()) - 1.0) < 1e-6
            h = mapmuts.bayesian.SiteEntropy(pis) # site entropy
            equilibriumfreqsfile.write('%d\t%s\t%f' % (ires, wtaa, h))
            for aa in mapmuts.sequtils.AminoAcids():
                equilibriumfreqsfile.write('\t%f' % pis[aa])
            equilibriumfreqsfile.write('\n')
            equilibriumfreqsfile.flush()
            if equilibriumfreqs_plots:
                plotfile = '%s/%s_%s%d.pdf' % (equilibriumfreqs_plots, outfileprefix, wtaa, ires)
                title = 'residue %s%d, site entropy of %.2f bits' % (wtaa, ires, h)
                mapmuts.plot.PlotEquilibriumFreqs(pis, plotfile, title)
    except:
        for x in sys.exc_info():
            log.write("\n\n%s" % str(x))
        log.write("\n\nPrematurely closing log due to execution error.")
        raise
    finally:
        log.write("\n\nExecution completed at %s." % time.ctime())
        log.close()
        enrichmentratiosfile.close()
        equilibriumfreqsfile.close()



if __name__ == '__main__':
    main() # run the script
