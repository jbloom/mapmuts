#!python

"""Counts number of mutations that occur >= a set number of times.

Designed to analyze how extensively a library samples all mutations, and
all synonymous mutations (at the codon level).

To run, type::

    mapmuts_countparsedmuts.py infile.txt

after creating the appropriate input file ``infile.txt``.
"""


import sys
import os
import time
import mapmuts.plot
import mapmuts.io


def main():
    """Main body of script."""

    print "Beginning execution of mapmuts_countparsedmuts.py..."

    if not mapmuts.plot.PylabAvailable():
        raise ImportError("This script requires matplotlib / pylab, which are not available.")

    mapmuts.io.PrintVersions(sys.stdout)
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument"\
                + ' specifying the name of the input file.')
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    lines = [line for line in open(infilename).readlines() if \
            line[0] != '#' and not line.isspace()]
    plotfileprefix = maxn = None
    legendloc = 'bottom' # default
    writecounts = True
    samples = [] # tuples are (name, all_counts, syn_counts)
    for line in lines:
        entries = line.split()
        if entries[0].strip() == 'plotfileprefix':
            if plotfileprefix != None:
                raise ValueError("Duplicate plotfileprefix keys")
            plotfileprefix = entries[1].strip()
        elif entries[0].strip() == 'maxn':
            if maxn != None:
                raise ValueError("Duplicate maxn key")
            if entries[1].strip().upper() == 'NONE':
                maxn = None
            else:
                try:
                    maxn = int(entries[1])
                except ValueError:
                    raise ValueError("maxn does not specify a valid integer: %s" % entries[1])
                if maxn < 1:
                    raise ValueError("max must be at least one")
        elif entries[0].strip() == 'legendloc':
            legendloc = entries[1].strip().lower()
            if legendloc not in ['bottom', 'right']:
                raise ValueError("legendloc must be either bottom or right, got: %s" % legendloc)
        elif entries[0].strip() == 'writecounts':
            if entries[1].strip().upper() == 'FALSE':
                writecounts = False
        else:
            if len(entries) < 2:
                raise ValueError("Line must contain at least two entries:\n%s" % line)
            name = entries[0].strip()
            counts = []
            for codoncountfile in entries[1 : ]:
                codoncountfile = codoncountfile.strip()
                if not os.path.isfile(codoncountfile):
                    raise IOError("Failed to find specified codon counts file of %s" % codoncountfile)
                print "Reading codon counts for %s from %s" % (name, codoncountfile)
                counts.append(mapmuts.io.ReadCodonCounts(open(codoncountfile)))
            (all_counts, multi_nt_all_counts, syn_counts, multi_nt_syn_counts) = mapmuts.sequtils.TallyCodonCounts(counts)
            if not all_counts:
                raise ValueError("No counts for %s" % name)
            if not (max(all_counts) >= max(syn_counts) >= 0):
                raise ValueError("Count minima don't make sense for %s" % name)
            samples.append((name, all_counts, multi_nt_all_counts, syn_counts, multi_nt_syn_counts))
    samples.sort()
    if not plotfileprefix:
        raise ValueError("Failed to parse a value for plotfileprefix")
    if not samples:
        raise ValueError("Failed to find any samples.")

    # get max occurrences of any mutation if maxn not specified
    if maxn == None:
        maxn = max(samples[0][1])
        for (name, all_counts, multi_nt_all_counts, syn_counts, multi_nt_syn_counts) in samples[1 : ]:
            maxn = max(maxn, max(all_counts))

    # now make cumul_samples: entries are (name, cumul_all, cumul_all_tot, cumul_multi_nt_all, cumul_multi_nt_all_tot, cumul_syn, cumul_syn_tot, cumul_multi_nt_syn, cumul_multi_nt_syn_tot)
    # where cumul_all[n] is the fraction that have >= n occurrences
    # where cumul_all_tot is total number of mutations in this category
    cumul_samples = []
    for (name, all_counts, multi_nt_all_counts, syn_counts, multi_nt_syn_counts) in samples:
        this_tuple = [name]
        for i_list in [all_counts, multi_nt_all_counts, syn_counts, multi_nt_syn_counts]:
            i_cumul = []
            i_d = {}
            for x in i_list:
                if x in i_d:
                    i_d[x] += 1
                else:
                    i_d[x] = 1
            ntot = float(len(i_list))
            nge = ntot
            for n in range(0, maxn):
                i_cumul.append(nge / ntot)
                if n in i_d:
                    nge -= i_d[n]
            i_cumul.append(nge / ntot)
            this_tuple.append(i_cumul)
            this_tuple.append(ntot)
            assert len(i_cumul) == maxn + 1, "len(i_cumul) = %d, maxn = %d" % (len(i_cumul), maxn)
        cumul_samples.append(tuple(this_tuple))

    # write the text files and make the plots
    cumulfracs = {}
    counts = {}
    for (tuple_index, mut_type) in [(1, 'all'), (3, 'multi-nt-all'), (5, 'syn'), (7, 'multi-nt-syn')]:
        # the text file
        fname = '%s_%scodonmutcounts.txt' % (plotfileprefix, mut_type)
        print "Now writing file %s..." % fname
        f = open(fname, 'w')
        f.write('# File listing the fraction of %s mutations that are found greater than or equal to n times.\n' % mut_type)
        f.write('# There are %d total %s mutations\n' % (cumul_samples[0][tuple_index + 1], mut_type))
        f.write('#n\t%s\n' % ('\t'.join([tup[0] for tup in cumul_samples])))
        for n in range(0, maxn + 1):
            f.write('%d' % n)
            for tup in cumul_samples:
                f.write('\t%.4f' % tup[tuple_index][n])
            f.write('\n')
        f.close()
        counts[mut_type] = cumul_samples[0][tuple_index + 1]
        names = [tup[0] for tup in cumul_samples]
        cumulfracs[mut_type] = [tup[tuple_index] for tup in cumul_samples]
    plotfile = "%s_multi-nt-codonmutcounts.pdf" % plotfileprefix
    print "Now creating plot %s..." % plotfile
    mapmuts.plot.PlotMutCountFracs(plotfile, 'Multi-nucleotide codon mutations', names, cumulfracs['multi-nt-all'], cumulfracs['multi-nt-syn'], counts['multi-nt-all'], counts['multi-nt-syn'], legendloc, writecounts=writecounts)
    plotfile = "%s_codonmutcounts.pdf" % plotfileprefix
    print "Now creating plot %s..." % plotfile
    mapmuts.plot.PlotMutCountFracs(plotfile, 'Codon mutations', names, cumulfracs['all'], cumulfracs['syn'], counts['all'], counts['syn'], legendloc)


    print "Script complete."


if __name__ == '__main__':
    main() # run the script
