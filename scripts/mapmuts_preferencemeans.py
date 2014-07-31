#!python

"""Computes mean of amino-acid preferences or differential preferences.

Written by Jesse Bloom.
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
    includestop = mapmuts.io.ParseBoolValue(d, 'includestop')
    if 'preferencefiles' in d and 'differentialpreferencefiles' in d:
        raise ValueError("Input file can only specify one of preferencefiles or differentialpreferencefiles, not both")
    elif 'preferencefiles' in d:
        use_diffprefs = False
        preferencefiles = mapmuts.io.ParseFileList(d, 'preferencefiles')
        log.write("\nReading in the preferences...")
        log.flush()
        preferences = [mapmuts.io.ReadEntropyAndEquilFreqs(f) for f in preferencefiles]
        residues = preferences[0].keys()
        residues.sort()
        if includestop and ('PI_*' in preferences[0][residues[0]]):
            includestop = True
        else:
            includestop = False
            for x in preferences:
                mapmuts.bayesian.PreferencesRemoveStop(x)
    elif 'differentialpreferencefiles' in d:
        use_diffprefs = True
        preferencefiles = mapmuts.io.ParseFileList(d, 'differentialpreferencefiles')
        log.write("\nReading in the differential preferences...")
        log.flush()
        preferences = [mapmuts.io.ReadDifferentialPreferences(f) for f in preferencefiles]
        residues = preferences[0].keys()
        residues.sort()
        if includestop and ('dPI_*' in preferences[0][residues[0]]):
            includestop = True
        else:
            includestop = False
            for x in preferences:
                mapmuts.bayesian.DifferentialPreferencesRemoveStop(x)
    else:
        raise ValueError("Input file failed to specify either preferencefiles or differentialpreferencefiles")
    nfiles = len(preferencefiles)
    if nfiles < 2:
        raise ValueError('preferencefiles / differentialpreferencefiles must specify at least 2 files')
    outfile = mapmuts.io.ParseStringValue(d, 'outfile')

    # read the preferences
    residues = preferences[0].keys()
    residues.sort()
    if not residues:
        raise ValueError("No residues specified in the first preference file")
    aas = mapmuts.sequtils.AminoAcids(includestop=includestop)
    log.write("\nNow writing means to %s..." % outfile)
    log.flush()
    out = open(outfile, 'w')
    if use_diffprefs:
        out.write('#SITE\tWT_AA\tRMS_dPI\t%s\n' % '\t'.join(['dPI_%s' % aa for aa in aas]))
    else:
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
            if includestop != (('PI_*' in preferences[ifile][r]) or ('dPI_*' in preferences[ifile][r])):
                raise ValueError("Not all files and residues are consistent with regard to the presence or absence of a key for stop codons (PI_* or dPI_*). All files and residues must either have or lack this key.")
        # start writing preference sums
        means = {}
        for aa in aas:
            if use_diffprefs:
                means[aa] = sum([preferences[ifile][r]['dPI_%s' % aa] for ifile in range(nfiles)]) / float(nfiles)
            else:
                means[aa] = sum([preferences[ifile][r]['PI_%s' % aa] for ifile in range(nfiles)]) / float(nfiles)
        if use_diffprefs and abs(sum(means.values())) > 1.0e-5:
            raise ValueError("The mean differential preferences do not sum to nearly zero for residue %d. Check that the preferences in the individual files correctly sum to zero." % r)
        elif not use_diffprefs and abs(1.0 - sum(means.values())) > 1.0e-5:
            raise ValueError("The mean differential preferences do not sum to nearly one for residue %d. Check that the preferences in the individual files correctly sum to one." % r)
        if use_diffprefs:
            h = math.sqrt(sum([dpi**2 for dpi in means.values()]))
        else:
            h = Entropy(means.values())
        out.write("%d\t%s\t%g\t%s\n" % (r, wts, h, '\t'.join(['%g' % means[aa] for aa in aas])))
    out.close()
    log.write("\n\nExecution completed at %s." % time.ctime())



if __name__ == '__main__':
    main() # run the script
