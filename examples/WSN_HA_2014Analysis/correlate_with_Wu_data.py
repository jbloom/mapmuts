"""Correlates amino-acid preferences with data from Wu et al 2014.

Written by Jesse Bloom."""


import mapmuts.io
import mapmuts.plot


def main():
    """Main body of script."""
    # input / output files and settings
    wudatafile = 'Wu_et_al_Table_S1.csv'
    aapreferencesfile = 'average_equilibriumpreferences.txt'
    plotfile = 'correlation_with_Wu_et_al.pdf'

    # Read in the data from Wu et al
    muts_by_res = {}
    nwt = nmut = nstop = 0
    for line in open(wudatafile).readlines()[1 : ]:
        (r, mut, rf) = line.split(',')
        r = int(r)
        wt = mut[0]
        if r not in muts_by_res:
            muts_by_res[r] = {'WT':wt}
        elif wt != muts_by_res[r]['WT']:
            raise ValueError("wildtype mismatch at residue %d" % r)
        m = mut[-1]
        if m == '_':
            nstop += 1
            continue
        elif m == wt:
            nwt += 1
        else:
            nmut += 1
        assert int(mut[1 : -1]) == r, "Mismatch at residue %d" % r
        rf = float(rf)
        assert (wt, m) not in muts_by_res[r], "Duplicate for %s" % line
        muts_by_res[r][(wt, m)] = rf
    print "Read from %s RF values for a total of %d identites for %d sites. Of these, %d were for wildtype identities, %d were for amino-acid mutations, and %d were for stop-codon mutations and were discarded." % (wudatafile, nwt + nmut + nstop, len(muts_by_res), nwt, nmut, nstop)
    wu_wtnormalizedrf = {}
    for (r, rd) in muts_by_res.iteritems():
        wt = muts_by_res[r]['WT']
        if (wt, wt) not in muts_by_res[r]:
            continue # no RF for wildtype
        wtrf = muts_by_res[r][(wt, wt)]
        for (key, rf) in muts_by_res[r].iteritems():
            if key not in ['WT', (wt, wt)]:
                m = key[1]
                wu_wtnormalizedrf['%s%d%s' % (wt, r, m)] = rf / float(wtrf)
    print "Computed a total of %d wildtype-normalized RF values." % len(wu_wtnormalizedrf)

    # Get the equivalent wildtype-normalized preference ratios
    preferences = mapmuts.io.ReadEntropyAndEquilFreqs(aapreferencesfile)
    preferences_normalizedrf = {}
    for (r, rd) in preferences.iteritems():
        wt = rd['WT_AA']
        for aa in mapmuts.sequtils.AminoAcids():
            if aa != wt:
                preferences_normalizedrf['%s%d%s' % (wt, r, aa)] = rd['PI_%s' % aa] / float(rd['PI_%s' % wt])
    print "\nComputed a total of %d wildtype normalized RF values from preferences in %s" % (len(preferences_normalizedrf), aapreferencesfile)

    # Now correlate
    wulist = []
    preflist = []
    for (mut, wurf) in wu_wtnormalizedrf.iteritems():
        if mut in preferences_normalizedrf:
            wulist.append(wurf)
            preflist.append(preferences_normalizedrf[mut])
    print "\nThere are paired values to correlate for %d mutations" % (len(wulist))
    nmuts = len(wulist)
    nonzeros = [rf for rf in wulist if rf]
    minvalue = min(nonzeros)
    print "The minimum non-zero value from Wu et al is %g; setting all zero values to this value" % minvalue
    wulist = [max(rf, minvalue) for rf in wulist]
    mapmuts.plot.PlotCorrelation(preflist, wulist, plotfile, 'our study', 'Wu \\textit{et al.}', logx=True, logy=True, corr='Pearson', alpha=0.1, symmetrize=True, bigmargin=0.27, xsize=2.2, title='Relative preference\n(%d of %d mutations)' % (nmuts, len(preferences_normalizedrf)))
    print "Plotted data is in %s" % plotfile


if __name__ == '__main__':
    main()
