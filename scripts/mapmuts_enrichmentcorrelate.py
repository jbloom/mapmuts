#!python

"""Determines correlations between enrichment ratios.

Written by Jesse Bloom, 2013.
"""


import sys
import os
import time
import warnings
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
    sys.stdout.write("\nBeginning execution of mapmuts_enrichmentcorrelate.py, reading input from %s.\n" % infilename)
    mapmuts.io.PrintVersions(sys.stdout)
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = mapmuts.io.ParseInfile(open(infilename))
    enrichmentratios = mapmuts.io.ParseFileList(d, 'enrichmentratios')
    limitadjustment = mapmuts.io.ParseFloatValue(d, 'limitadjustment')
    if limitadjustment < 1:
        raise ValueError("limitadjustment should be >= 1")
    samplenames = mapmuts.io.ParseStringValue(d, 'samplenames')
    samplenames = [x.strip() for x in samplenames.split()]
    if len(enrichmentratios) != len(samplenames):
        raise ValueError("samplenames and enrichmentratios do not specify the same number of entries")
    plotdir = mapmuts.io.ParseStringValue(d, 'plotdir')
    if not os.path.isdir(plotdir):
        raise IOError("plotdir directory of %s does not already exist. You must create it before running this script." % plotdir)
    mutdnacounts_cutoff = mapmuts.io.ParseStringValue(d, 'mutdnacounts_cutoff')
    cutoffs = []
    for x in mutdnacounts_cutoff.split():
        try:
            x = int(x)
        except ValueError:
            raise ValueError("mutdnacounts_cutoff specifies non-integer values, such as %s" % x)
        if x < 0:
            raise ValueError("mutdnacounts_cutoff specifies integers less than zero; they must be >= 0")
        cutoffs.append(x)
    enrichmentratios = [mapmuts.io.ReadEnrichmentRatios(f) for f in enrichmentratios]
    # apply limitadjustment
    directratios = []
    for d in enrichmentratios:
        for mut_d in d.itervalues():
            directratio = mut_d['DIRECT_RATIO']
            if directratio > 0 and directratio != float('inf'):
                directratios.append(directratio)
    (mindirectratio, maxdirectratio) = (min(directratios), max(directratios))
    for d in enrichmentratios:
        for mut_d in d.itervalues():
            directratio = mut_d['DIRECT_RATIO']
            if directratio <= 0:
                mut_d['DIRECT_RATIO'] = mindirectratio / float(limitadjustment)
            if directratio == float('inf'):
                mut_d['DIRECT_RATIO'] = maxdirectratio * limitadjustment
    # make the inferred_vs_direct plots
    for (d, sample) in zip(enrichmentratios, samplenames):
        phis = []
        direct = []
        counts = []
        for mut_d in d.itervalues():
            phis.append(mut_d['PHI'])
            direct.append(mut_d['DIRECT_RATIO'])
            counts.append(mut_d['MUTDNA_COUNTS'])
        n = len(phis)
        for cutoff in cutoffs:
            xs = [phis[i] for i in range(n) if counts[i] >= cutoff]
            ys = [direct[i] for i in range(n) if counts[i] >= cutoff]
            plotfile = '%s/%s_inferred_vs_direct_cutoff%d.pdf' % (plotdir, sample, cutoff)
            if cutoff > 0:
                title = '%s (count cutoff of %d)' % (sample.replace('_', ' '), cutoff)
            else:
                title = '%s (no count cutoff)' % sample.replace('_', ' ')
            sys.stdout.write("Making plot %s...\n" % plotfile)
            mapmuts.plot.PlotCorrelation(xs, ys, plotfile, xlabel='enrichment ratio (inferred)', ylabel='direct ratio', logx=True, logy=True, corr=corr, title=title)
    # make the replicate comparison plots
    for isample in range(len(samplenames)):
        (d1, sample1) = (enrichmentratios[isample], samplenames[isample])
        for (d2, sample2) in zip(enrichmentratios, samplenames)[isample + 1 : ]:
            sharedmuts = [mut for mut in d1.iterkeys() if mut in d2]
            for (plot_type, key, titlelabel) in [('inferred', 'PHI', 'inferred enrichment ratio'), ('direct', 'DIRECT_RATIO', 'direct ratio')]:
                for cutoff in cutoffs:
                    sample1data = []
                    sample2data = []
                    for mut in sharedmuts:
                        if d1[mut]['MUTDNA_COUNTS'] >= cutoff and d2[mut]['MUTDNA_COUNTS'] >= cutoff:
                            sample1data.append(d1[mut][key])
                            sample2data.append(d2[mut][key])
                    plotfile = '%s/%s_vs_%s_%s_cutoff%d.pdf' % (plotdir, sample1, sample2, plot_type, cutoff)
                    if cutoff > 0:
                        title = '%s,\ncount cutoff of %d' % (titlelabel, cutoff)
                    else:
                        title = '%s,\nno count cutoff' % titlelabel
                    sys.stdout.write("Making plot %s...\n" % plotfile)
                    mapmuts.plot.PlotCorrelation(sample1data, sample2data, plotfile, xlabel=sample1.replace('_', ' '), ylabel=sample2.replace('_', ' '), logx=True, logy=True, corr=corr, title=title)
    sys.stdout.write("\nScript execution completed.\n")



if __name__ == '__main__':
    main() # run the script
