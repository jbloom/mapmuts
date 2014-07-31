#!python

"""Plots correlations between amino-acid preferences or differential preferences.

Written by Jesse Bloom.
"""


import sys
import os
import time
import warnings
import copy
import mapmuts.io
import mapmuts.sequtils
import mapmuts.plot
import mapmuts.bayesian


def main():
    """Main body of script."""

    # check on module availability
    if not mapmuts.plot.PylabAvailable():
        raise ImportError("Cannot run this script as pylab is not available.")
    if not mapmuts.bayesian.ScipyAvailable():
        warnings.warn("Cannot import scipy. Pearson correlation coefficients cannot be computed without this package.\n")
        corr = None
    else:
        corr = 'Pearson'

    # read input variables
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument"\
                + ' specifying the name of the input file.')
    infilename = sys.argv[1]
    sys.stdout.write("\nBeginning execution of mapmuts_preferencescorrelate.py, reading input from %s.\n" % infilename)
    mapmuts.io.PrintVersions(sys.stdout)
    sys.stdout.flush()
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = mapmuts.io.ParseInfile(open(infilename))
    if 'preferencesfiles' in d and 'differentialpreferencesfiles' in d:
        raise ValueError("You can specify only one of preferencesfile or differentialpreferencesfiles, not both")
    elif 'preferencesfiles' in d:
        preferencesfiles = mapmuts.io.ParseFileList(d, 'preferencesfiles')
        preferences = [mapmuts.io.ReadEntropyAndEquilFreqs(f) for f in preferencesfiles]
        use_diffprefs = False
    elif 'differentialpreferencesfiles' in d:
        preferencesfiles = mapmuts.io.ParseFileList(d, 'differentialpreferencesfiles')
        preferences = [mapmuts.io.ReadDifferentialPreferences(f) for f in preferencesfiles]
        use_diffprefs = True
    else: 
        raise ValueError("Failed to specify either preferencesfiles or differentialpreferencesfiles")
    if len(preferencesfiles) < 2:
        raise ValueError("Failed to specify at least two preferencesfiles or differentialpreferencesfiles")
    samplenames = mapmuts.io.ParseStringValue(d, 'samplenames')
    if 'alpha' in d:
        alpha = mapmuts.io.ParseFloatValue(d, 'alpha')
        if not 0 < alpha <= 1:
            raise ValueError("alpha must be > 0 and <= 1")
    else:
        alpha = 1.0
    if 'logscale' in d:
        logscale = mapmuts.io.ParseBoolValue(d, 'logscale')
    else:
        logscale = False
    if 'plot_simpsondiversity' in d:
        plot_simpsondiversity = mapmuts.io.ParseBoolValue(d, 'plot_simpsondiversity')
    else:
        plot_simpsondiversity = False
    if 'plot_RMSdiffpref' in d:
        plot_RMSdiffpref = mapmuts.io.ParseBoolValue(d, 'plot_RMSdiffpref')
    else:
        plot_RMSdiffpref = False
    samplenames = [x.strip() for x in samplenames.split()]
    if len(preferencesfiles) != len(samplenames):
        raise ValueError("samplenames and preferencesfiles/differentialpreferencesfiles do not specify the same number of entries")
    plotdir = mapmuts.io.ParseStringValue(d, 'plotdir')
    if not os.path.isdir(plotdir):
        raise IOError("plotdir directory of %s does not already exist. You must create it before running this script." % plotdir)

    # make the sample comparison plots
    for isample in range(len(samplenames)):
        (d1, sample1) = (preferences[isample], samplenames[isample])
        for (d2, sample2) in zip(preferences, samplenames)[isample + 1 : ]:
            plotfile = "%s/%s_vs_%s.pdf" % (plotdir, sample1, sample2)
            sys.stdout.write("Making plot %s...\n" % plotfile)
            sys.stdout.flush()
            sharedresidues = [r for r in d1.iterkeys() if r in d2]
            if not sharedresidues:
                raise ValueError("No shared residues")
            if (('PI_*' in d1[sharedresidues[0]]) and ('PI_*' in d2[sharedresidues[0]])) or (('dPI_*' in d1[sharedresidues[0]]) and ('dPI_*' in d2[sharedresidues[0]])):
                includestop = True
            else:
                includestop = False
                if 'PI_*' in d1[sharedresidues[0]] or 'dPI_*' in d1[sharedresidues[0]]:
                    d1 = copy.deepcopy(d1)
                    if use_diffprefs:
                        mapmuts.bayesian.DifferentialPreferencesRemoveStop(d1)
                    else:
                        mapmuts.bayesian.PreferencesRemoveStop(d1)
                if 'PI_*' in d2[sharedresidues[0]] or 'dPI_*' in d2[sharedresidues[0]]:
                    d2 = copy.deepcopy(d2)
                    if use_diffprefs:
                        mapmuts.bayesian.DifferentialPreferencesRemoveStop(d2)
                    else:
                        mapmuts.bayesian.PreferencesRemoveStop(d2)
            sample1data = []
            sample2data = []
            aas = mapmuts.sequtils.AminoAcids(includestop=includestop)
            for r in sharedresidues:
                if plot_simpsondiversity and not use_diffprefs:
                    pi1 = dict([(aa, d1[r]['PI_%s' % aa]) for aa in aas])
                    pi2 = dict([(aa, d2[r]['PI_%s' % aa]) for aa in aas])
                    sample1data.append(mapmuts.bayesian.GiniSimpson(pi1))
                    sample2data.append(mapmuts.bayesian.GiniSimpson(pi2))
                elif plot_RMSdiffpref and use_diffprefs:
                    sample1data.append(d1[r]['RMS_dPI'])
                    sample2data.append(d2[r]['RMS_dPI'])
                else:
                    for aa in aas:
                        if use_diffprefs:
                            sample1data.append(d1[r]['dPI_%s' % aa])
                            sample2data.append(d2[r]['dPI_%s' % aa])
                        else:
                            sample1data.append(d1[r]['PI_%s' % aa])
                            sample2data.append(d2[r]['PI_%s' % aa])
            if logscale:
                if min(sample1data + sample2data) <= 0:
                    raise ValueError("Cannot use logscale as there is data <= 0")
                mapmuts.plot.PlotCorrelation(sample1data, sample2data, plotfile, xlabel=sample1.replace('_', ' '), ylabel=sample2.replace('_', ' '), logx=True, logy=True, corr=corr, alpha=alpha, symmetrize=False, fixaxes=False)
            else:
                fixaxes = not use_diffprefs
                mapmuts.plot.PlotCorrelation(sample1data, sample2data, plotfile, xlabel=sample1.replace('_', ' '), ylabel=sample2.replace('_', ' '), corr=corr, alpha=alpha, symmetrize=True, fixaxes=fixaxes)

    sys.stdout.write("\nScript execution completed.\n")



if __name__ == '__main__':
    main() # run the script
