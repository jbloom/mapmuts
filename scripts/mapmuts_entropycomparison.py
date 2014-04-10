#!python

"""Compares site entropies for two sets of sites controlling for solvent accessibility.

Jesse Bloom, 2014."""


import sys
import os
import random
import rpy2.robjects
import mapmuts.plot
import mapmuts.bayesian
import mapmuts.io
import mapmuts.dssp
import mapmuts.sequtils


def main():
    """Main body of script."""
    if not mapmuts.plot.PylabAvailable():
        raise ImportError("Cannot run this script as pylab / matplotlib are not available")

    # read input variables
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the input file")
    infilename = sys.argv[1]
    sys.stdout.write("\nBeginning execution of mapmuts_entropycomparison.py, reading input from %s" % infilename)
    mapmuts.io.PrintVersions(sys.stdout)
    if not os.path.isfile(infilename):
        raise IOError("Cannot find input file %s" % infilename)
    d = mapmuts.io.ParseInfile(open(infilename))
    aapreferences = mapmuts.io.ParseFileList(d, 'aapreferences')
    if len(aapreferences) != 1:
        raise ValueError("aapreferences failed to specify exactly one file")
    aapreferences = mapmuts.io.ReadEntropyAndEquilFreqs(aapreferences[0])
    dsspfile = mapmuts.io.ParseFileList(d, 'dsspfile')
    if len(dsspfile) != 1:
        raise ValueError("dsspfile failed to specify exactly one file")
    dsspchain = mapmuts.io.ParseStringValue(d, 'dsspchain')
    if dsspchain.upper() in ['NONE', 'FALSE']:
        dsspchain = None
    dssp = mapmuts.dssp.ReadDSSP(dsspfile[0], 'Tien2013', chain=dsspchain)
    siterange = mapmuts.io.ParseStringValue(d, 'siterange')
    if siterange.upper() == 'ALL':
        sites = aapreferences.keys()
    else:
        tup = siterange.split()
        if len(tup) != 2:
            raise ValueError("Invalid siterange")
        siterange = (int(tup[0]), int(tup[1]))
        if siterange[1] <= siterange[0]:
            raise ValueError("Empty siterange")
        sites = [r for r in range(siterange[0], siterange[1] + 1) if r in aapreferences]
    if 'includestop' in d:
        includestop = mapmuts.io.ParseBoolValue(d, 'includestop')
    else:
        includestop = False
    if includestop:
        mapmuts.bayesian.PreferencesRemoveStop(aapreferences)
    if 'PI_*' in aapreferences[sites[0]]:
        aminoacids = mapmuts.sequtils.AminoAcids(includestop=True)
    else:
        aminoacids = mapmuts.sequtils.AminoAcids(includestop=False)
    plotfile = mapmuts.io.ParseStringValue(d, 'plotfile')
    linearmodelfile = mapmuts.io.ParseStringValue(d, 'linearmodelfile')
    if 'alpha' in d:
        alpha = mapmuts.io.ParseFloatValue(d, 'alpha')
        if not (0 < alpha <= 1.0):
            raise ValueError("alpha must be > 0 and <= 1")
    else:
        alpha = 1.0
    selectedsites = mapmuts.io.ParseStringValue(d, 'selectedsites')
    if not os.path.isfile(selectedsites):
            raise ValueError("selectedsites does not specify a valid file, cannot find %s" % selectedsites)
    lines = [line for line in open(selectedsites).readlines() if (not line.isspace()) and line[0] != '#']
    selectedsites = [int(line.split()[0]) for line in lines]

    # make plot
    entropies = {}
    for r in sites:
        pir = dict([(aa, aapreferences[r]['PI_%s' % aa]) for aa in aminoacids])
        entropies[r] = mapmuts.bayesian.SiteEntropy(pir)
    rsas  = dict([(r, dssp[r]['RSA']) for r in sites if r in dssp])
    sharedsites = [r for r in sites if (r in entropies) and (r in dssp)]
    selectedsites = [r for r in selectedsites if r in sharedsites]
    if len(sharedsites) < 2:
        raise ValueError("Failed to find at least two sites with both RSA and amino-acid preferences specified")
    entropy_list = [entropies[r] for r in sharedsites if r not in selectedsites]
    rsa_list = [rsas[r] for r in sharedsites if r not in selectedsites]
    if selectedsites:
        selected_entropy_list = [entropies[r] for r in selectedsites]
        selected_rsa_list = [rsas[r] for r in selectedsites]
        additionalxy = [(selected_rsa_list, selected_entropy_list)]
    else:
        additionalxy = []
    sys.stdout.write("\nPlotting the correlation in %s\n" % plotfile)
    mapmuts.plot.PlotCorrelation(rsa_list, entropy_list, plotfile, xlabel='RSA', ylabel='site entropy', corr=None, alpha=alpha, additionalxy=additionalxy)

    # do multiple linear regression
    r_entropy = rpy2.robjects.FloatVector([entropies[r] for r in sharedsites])
    r_rsa = rpy2.robjects.FloatVector([rsas[r] for r in sharedsites])
    r_antigenic = rpy2.robjects.FloatVector([int(r in selectedsites) for r in sharedsites])
    rpy2.robjects.globalenv['r_entropy'] = r_entropy
    rpy2.robjects.globalenv['r_rsa'] = r_rsa
    rpy2.robjects.globalenv['r_antigenic'] = r_antigenic
    r_linmodel = rpy2.robjects.r['lm']
    results = r_linmodel(formula='r_entropy ~ r_rsa + r_antigenic')
    base = rpy2.robjects.packages.importr('base')
    sys.stdout.write("\nWriting the results of the linear model analysis to %s\n" % linearmodelfile)
    old_out = sys.stdout
    sys.stdout = open(linearmodelfile, 'w')
    print(base.summary(results))
    sys.stdout.close()
    sys.stdout = old_out

    sys.stdout.write('\nScript execution complete.\n')



if __name__ == '__main__':
    main() # run the main program
