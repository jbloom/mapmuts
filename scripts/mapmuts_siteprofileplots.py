#!python

"""Plots profiles of entropies and preferences along primary sequence.

Written by Jesse Bloom, 2013.
"""


import sys
import os
import time
import warnings
import mapmuts.io
import mapmuts.sequtils
import mapmuts.plot
import mapmuts.dssp
import mapmuts.bayesian
import mapmuts.weblogo


def main():
    """Main body of script."""

    # check on module availability
    if not mapmuts.weblogo.WebLogoAvailable():
        raise ImportError("Cannot run this script as weblogo is not available.")
    if not mapmuts.weblogo.PyPdfAvailable():
        raise ImportError("Cannot run this script as PyPdf is not available.")
    if not mapmuts.plot.PylabAvailable():
        raise ImportError("Cannot run this script at pylab is not available.")

    # read input variables
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument"\
                + ' specifying the name of the input file.')
    infilename = sys.argv[1]
    sys.stdout.write("\nBeginning execution of mapmuts_siteprofileplots.py, reading input from %s.\n" % infilename)
    mapmuts.io.PrintVersions(sys.stdout)
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = mapmuts.io.ParseInfile(open(infilename))
    sitepreferences = mapmuts.io.ParseFileList(d, 'sitepreferences')
    if len(sitepreferences) != 1:
        raise ValueError("sitepreferences did not specify exactly one file")
    sitepreferences = sitepreferences[0]
    outfileprefix = mapmuts.io.ParseStringValue(d, 'outfileprefix')
    siterange = mapmuts.io.ParseStringValue(d, 'siterange')
    if siterange.upper() == 'ALL':
        siterange = None
    else:
        tup = siterange.split()
        if len(tup) != 2:
            raise ValueError("Invalid siterange")
        siterange = (int(tup[0]), int(tup[1]))
        if siterange[1] <= siterange[0]:
            raise ValueError("Empty siterange")
    dsspfile = mapmuts.io.ParseStringValue(d, 'dsspfile')
    if dsspfile.upper() in ['NONE', 'FALSE']:
        dsspfile = False
    else:
        dsspfile = dsspfile.strip()
        if not os.path.isfile(dsspfile):
            raise IOError("Cannot find dsspfile of %s" % dsspfile)
        dsspchain = mapmuts.io.ParseStringValue(d, 'dsspchain')
        if dsspchain.upper() in ['NONE', 'FALSE']:
            dsspchain = None
    add_rsa = mapmuts.io.ParseBoolValue(d, 'add_rsa')
    if add_rsa and not dsspfile:
        raise ValueError("Cannot use add_rsa without dsspfile")
    add_ss = mapmuts.io.ParseBoolValue(d, 'add_ss')
    if add_ss and not dsspfile:
        raise ValueError("Cannot use add_ss without dsspfile")
    nperline = mapmuts.io.ParseIntValue(d, 'nperline')
    if nperline < 1:
        raise ValueError("nperline must be an integer >= 1")
    includestop = mapmuts.io.ParseBoolValue(d, 'includestop')
    sitenumbermapping = None
    if 'sitenumbermapping' in d:
        sitenumbermapping = mapmuts.io.ParseStringValue(d, 'sitenumbermapping')
        if sitenumbermapping.upper() == 'NONE':
            sitenumbermapping = None
        else:
            if not os.path.isfile(sitenumbermapping):
                raise IOError("Failed to find sitenumbermapping file %s" % sitenumbermapping)
            sitenumbermapping = dict([(int(line.split(',')[0]), line.split(',')[1]) for line in open(sitenumbermapping).readlines() if line[0] != '#' and not line.isspace()])

#    add_entropy = mapmuts.io.ParseBoolValue(d, 'add_entropy')

    # start running script
    d = mapmuts.io.ReadEntropyAndEquilFreqs(sitepreferences)
    if not siterange:
        siterange = (min(d.keys()), max(d.keys()))
    sites = [site for site in range(siterange[0], siterange[1] + 1)]
    if not includestop:
        mapmuts.bayesian.PreferencesRemoveStop(d)
    if sitenumbermapping:
        for site in sites:
            if site not in sitenumbermapping:
                raise ValueError("sitenumbermapping file fails to specify a value for site %d" % site)

    # make plots of entropies along primary sequence
    try:
        entropies = [d[site]['SITE_ENTROPY'] for site in sites]
    except KeyError:
        raise KeyError("The sitepreferences file may not contain all of the consecutive required site numbers")
    plotfile = "%ssite_entropy_plot.pdf" % outfileprefix
    sys.stdout.write('Creating plot %s...\n' % plotfile)
    mapmuts.plot.PlotLinearDensity([('site entropy', zip(sites, entropies))], plotfile, xlabel='residue number', ylabel='site entropy (bits)')

    # plot entropy / RSA correlation
    if dsspfile:
        dsspplotfile = '%sentropy_rsa_correlationplot.pdf' % outfileprefix
        sys.stdout.write('Creating plot %s...\n' % dsspplotfile)
        dssp = mapmuts.dssp.ReadDSSP(dsspfile, 'Tien2013', chain=dsspchain)
        if mapmuts.bayesian.ScipyAvailable():
            corr = 'Pearson'
        else:
            warnings.warn("Will not be able to display correlation on %s as scipy is unavailable.\n" % dsspplotfile)
            corr = None
        entropies = [d[site]['SITE_ENTROPY'] for site in sites if site in dssp]
        rsas = [dssp[site]['RSA'] for site in sites if site in dssp]
        mapmuts.plot.PlotCorrelation(entropies, rsas, dsspplotfile, xlabel='site entropy (bits)', ylabel='RSA', corr=corr)

    # make dictionaries defining other properties
    otherprops = []
    if add_rsa:
        rsa_d = dict([(site, dssp[site]['RSA']) for site in sites if site in dssp])
        otherprops.append(('RSA', rsa_d))
    if add_ss:
        ss_d = dict([(site, dssp[site]['SS_CLASS']) for site in sites if site in dssp])
        ss_categories = ['helix', 'strand', 'loop']
        otherprops.append(('SS', ss_d, ss_categories))
#    if add_entropy:
#        otherprops.append(('entropy', dict([(site, d[site]['SITE_ENTROPY']) for site in sites])))
    # make heatmap CURRENTLY THIS IS COMMENTED OUT AND NOT DONE
#    heatmapfile = '%ssite_preferences_heatmap.pdf' % outfileprefix
#    sys.stdout.write('Creating the heat map plot %s...\n' % heatmapfile)
#    mapmuts.plot.EquilibriumFreqsHeatMap(sites, d, 'hydrophobicity', heatmapfile, otherprops=otherprops)

    # make sequence logo plot
    logoplotfile = '%ssite_preferences_logoplot.pdf' % outfileprefix
    sys.stdout.write("Creating plot %s...\n" % logoplotfile)
    if add_rsa and add_ss:
        overlay = [rsa_d, ss_d]
    else:
        overlay = None
    mapmuts.weblogo.EquilibriumFreqsLogo(sites, d, logoplotfile, nperline=nperline, overlay=overlay, sitenumbermapping=sitenumbermapping)

    sys.stdout.write("Script execution completed.\n")



if __name__ == '__main__':
    main() # run the script
