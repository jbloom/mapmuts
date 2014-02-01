#!python

"""Calculates geometric means of enrichment ratios.

Written by Jesse Bloom, 2013.
"""


import re
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
    print "Beginning execution of mapmuts_enrichmentgeomeans.py..."
    mapmuts.io.PrintVersions(sys.stdout)
    if not mapmuts.plot.PylabAvailable():
        raise ImportError("Cannot run this script as matplotlib / pylab are not available.")
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument"\
                + ' specifying the name of the input file.')
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = mapmuts.io.ParseInfile(open(infilename))
    outfileprefix = mapmuts.io.ParseStringValue(d, 'outfileprefix')
    equilibriumfreqsfile = open('%s_equilibriumfreqs.txt' % outfileprefix, 'w')
    equilibriumfreqsfile.write('#SITE\tWT_AA\tSITE_ENTROPY')
    for aa in mapmuts.sequtils.AminoAcids():
        equilibriumfreqsfile.write('\tPI_%s' % aa)
    equilibriumfreqsfile.write('\n')
    enrichmentratiosfile = open('%s_enrichmentratios.txt' % outfileprefix, 'w')
    enrichmentratiosfile.write('#MUTATION\tGEOMEAN_PHI\n')
    enrichmentratios = mapmuts.io.ParseFileList(d, 'enrichmentratios')
    enrichmentratios = [mapmuts.io.ReadEnrichmentRatios(f) for f in enrichmentratios]
    # parse any excludesite_NNN entries
    excludekeymatch = re.compile('excludesite\_(?P<site>\d+)')
    excludekeys = [key for key in d.iterkeys() if excludekeymatch.search(key)]
    exclude_libs = {} # keyed by site, values dictionary keyed by excluded num (0, 1, 2, ..)
    for key in excludekeys:
        site = int(excludekeymatch.search(key).group('site'))
        exclude_libs[site] = {}
        entries = d[key].split()
        if len(entries) != len(enrichmentratios):
            raise ValueError("Invalid number of entries for %s" % key)
        for i in range(len(entries)):
            if entries[i].strip() == 'use':
                continue
            elif entries[i].strip() == 'exclude':
                exclude_libs[site][i] = True                    
            else:
                raise ValueError("Entries for %s must be 'use' or 'exclude', not %s" % (key, entries[i]))
    enrichmentratio_plots = mapmuts.io.ParseStringValue(d, 'enrichmentratio_plots')
    if enrichmentratio_plots in ['None', 'False']:
        enrichmentratio_plots = None
    elif not os.path.isdir(enrichmentratio_plots):
        raise IOError("enrichmentratio_plots directory of %s does not already exist. You must create it before running this script." % enrichmentratio_plots)
    equilibriumfreqs_plots = mapmuts.io.ParseStringValue(d, 'equilibriumfreqs_plots')
    if equilibriumfreqs_plots in ['None', 'False']:
        equilibriumfreqs_plots = None
    elif not os.path.isdir(equilibriumfreqs_plots):
        raise IOError("equilibriumfreqs_plots directory of %s does not already exist. You must create it before running this script." % equilibriumfreqs_plots)
    # get the list of all sites and the wildtype identities
    mutmatch = re.compile('^(?P<wt>[A-Z\*])(?P<site>\d+)(?P<mut>[A-Z\*])$')
    wt_by_site = {}
    ilib = 0
    for d in enrichmentratios:
        for key in d.iterkeys():
            m = mutmatch.search(key)
            if not m:
                raise ValueError("Failed to match mutation %s" % key)
            wt = m.group('wt')
            site = int(m.group('site'))
            if site in exclude_libs and ilib in exclude_libs[site]:
                continue
            if site in wt_by_site:
                if wt_by_site[site] != wt:
                    raise ValueError("Discrepancy on wildtype identity at site %d" % site)
            else:
                wt_by_site[site] = wt
        ilib += 1
    (minsite, maxsite) = (min(wt_by_site.keys()), max(wt_by_site.keys()))
    sites = []
    for site in range(minsite, maxsite + 1):
        sites.append(site)
        if site not in wt_by_site:
            raise ValueError("No data from any libraries for site %d" % site)
    # calculate the enrichment ratio geometric means
    for site in sites:
        sys.stdout.write("Performing analysis for site %d...\n" % site)
        sys.stdout.flush()
        wt = wt_by_site[site]
        site_phi = {}
        mutlist = []
        philist = []
        for mut in [aa for aa in mapmuts.sequtils.AminoAcids() + ['*'] if aa != wt]:
            mut_phis = []
            key = "%s%d%s" % (wt, site, mut)
            for ilib in range(len(enrichmentratios)):
                if site in exclude_libs and ilib in exclude_libs[site]:
                    continue
                else:
                    mut_phis.append(enrichmentratios[ilib][key]['PHI'])
            if not mut_phis:
                raise ValueError("No enrichment ratios for %s" % key)
            geomean = 1.0
            for phi in mut_phis:
                geomean *= phi
            geomean = geomean**(1.0 / len(mut_phis))
            mutlist.append(key)
            philist.append(geomean)
            enrichmentratiosfile.write('%s\t%f\n' % (key, geomean))
            if mut != '*':
                site_phi[mut] = geomean
        pi = mapmuts.bayesian.EquilibriumFracs(wt, site_phi)
        h = mapmuts.bayesian.SiteEntropy(pi)
        assert len(pi) == 20 and abs(sum(pi.values()) - 1.0) < 1e-6
        equilibriumfreqsfile.write('%d\t%s\t%f' % (site, wt, h))
        for aa in mapmuts.sequtils.AminoAcids():
            equilibriumfreqsfile.write('\t%f' % pi[aa])
        equilibriumfreqsfile.write('\n')
        if enrichmentratio_plots:
            plotfile = '%s/%s_%s%d.pdf' % (enrichmentratio_plots, outfileprefix, wt, site)
            mapmuts.plot.PlotEnrichmentRatios(mutlist, philist, None, plotfile)
        if equilibriumfreqs_plots:
            plotfile = '%s/%s_%s%d.pdf' % (equilibriumfreqs_plots, outfileprefix, wt, site)
            title = 'residue %s%d, site entropy of %.2f bits' % (wt, site, h)
            mapmuts.plot.PlotEquilibriumFreqs(pi, plotfile, title)
    enrichmentratiosfile.close()
    equilibriumfreqsfile.close()
    print "Script complete."



if __name__ == '__main__':
    main() # run the script
