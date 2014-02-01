#!python

"""Computes mean of amino-acid preferences from several files.

Written by Jesse Bloom, 2013.
"""


import sys
import os
import math
import time
import mapmuts
import mapmuts.io
import mapmuts.sequtils
import mapmuts.bayesian



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



def main():
    """Main body of script."""
    # hard-coded variables
    log = sys.stdout
    # read input variables
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument"\
                + ' specifying the name of the input file.')
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = mapmuts.io.ParseInfile(open(infilename))
    log.write("Beginning execution of mapmuts_preferencemeans.py"\
            " in directory %s" % (os.getcwd()))
    mapmuts.io.PrintVersions(log)
    log.write("Input data being read from infile %s\n\n" % infilename)
    log.write("Read the following key/value pairs from infile %s:" % infilename)
    for (key, value) in d.iteritems():
        log.write("\n%s %s" % (key, value))
    preferencefiles = mapmuts.io.ParseFileList(d, 'preferencefiles')
    nfiles = len(preferencefiles)
    if nfiles < 2:
        raise ValueError('preferencefiles must specify at least 2 files')
    outfile = mapmuts.io.ParseStringValue(d, 'outfile')
    includestop = mapmuts.io.ParseBoolValue(d, 'includestop')

    # read the preferences
    log.write("\nReading in the preferences...")
    log.flush()
    preferences = [mapmuts.io.ReadEntropyAndEquilFreqs(f) for f in preferencefiles]
    residues = preferences[0].keys()
    residues.sort()
    if not residues:
        raise ValueError("No residues specified in the first preference file")
    if includestop and ('PI_*' in preferences[0][residues[0]]):
        includestop = True
    else:
        includestop = False
        for x in preferences:
            mapmuts.bayesian.PreferencesRemoveStop(x)
    aas = mapmuts.sequtils.AminoAcids(includestop=includestop)
    log.write("\nNow writing mean preferences to %s..." % outfile)
    log.flush()
    out = open(outfile, 'w')
    out.write('#SITE\tWT_AA\tSITE_ENTROPY\t%s\n' % '\t'.join(['PI_%s' % aa for aa in aas]))
    for r in residues:
        # get wildtype(s) for all files
        wts = [preferences[ifile][r]['WT_AA'] for ifile in range(nfiles)]
        if len(dict([(wt, True) for wt in wts])) == 1:
            wts = wts[0]
        else:
            wts = ','.join(wts)
        # check that all files are consistent with stop codon presence/absence
        for ifile in range(nfiles):
            if includestop != ('PI_*' in preferences[ifile][r]):
                raise ValueError("Not all files and residues are consistent with regard to the presence or absence of a key for stop codons (PI_*). All files and residues must either have or lack this key.")
        # start writing preference sums
        meanpis = {}
        for aa in aas:
            meanpis[aa] = sum([preferences[ifile][r]['PI_%s' % aa] for ifile in range(nfiles)]) / float(nfiles)
        if abs(1.0 - sum(meanpis.values())) > 1.0e-5:
            raise ValueError("The mean preferences do not sum to nearly one for residue %d. Check that the preferences in the individual files correctly sum to one." % r)
        h = Entropy(meanpis.values())
        out.write("%d\t%s\t%g\t%s\n" % (r, wts, h, '\t'.join(['%g' % meanpis[aa] for aa in aas])))
    out.close()
    log.write("\n\nExecution completed at %s." % time.ctime())



if __name__ == '__main__':
    main() # run the script
