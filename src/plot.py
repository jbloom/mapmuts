"""Module for performing plotting for ``mapmuts`` package.

This module uses ``pylab`` and ``matplotlib`` to make plots. These plots will
fail if ``pylab`` and ``matplotlib`` are not available for importation. Before
running any function in this module, you can run the `PylabAvailable`
function to determine if ``pylab`` and ``matplotlib`` are available. Otherwise,
calling any other function will raise an Exception if thise modules are
not available. The ``pdf`` backend is used for ``matplotlib`` / ``pylab``. This means
that plots must be created as PDF files.

The *PlotCorrelation* function can use ``scipy`` to calculate correlation
coefficients if this package is available.

Written by Jesse Bloom.


List of functions
-------------------
`PylabAvailable` : tests whether ``pylab`` and ``matplotlib`` can be imported.

`PlotAlignmentStatistics`

`PlotInsertLengths`

`PlotReadMismatches`

`PlotCodonDepth`

`PlotCodonDist`

`PlotMutFracs`

`PlotPairedMutFracs`

`PlotNTMutFracs`

`PlotAAFracs`

`PlotMutCountFracs` : plots fraction mutations with at least a certain # of counts.

`PlotEnrichmentRatios` : plots enrichment ratios on log scale

`PlotEquilibriumFreqs` : plots equilibrium frequencies of amino acids

`EquilibriumFreqsHeatMap` : plots heat map of equilibrium amino-acid frequencies.

`EquilibriumFreqsLogo` : sequence logo of amino-acid frequencies - DOES NOT WORK!

`LogoOverlay` : color maps of RSA and SS to be overlayed on weblogo sequence logo

`PlotTraces` : plots traces of posterior probabilities

`PlotCorrelation` : plots correlation between two variables.

`PlotLinearDensity` : plots one or more variables along primary sequence.

`Base10Formatter`

`SplitLabel`

`KyteDoolittleColorMapping` : maps Kyte-Doolittle scale to colors.


Documentation for functions
-----------------------------
Documentation for individual functions is provided in their definitions below.

"""


import os
import sys
import math
import warnings
import mapmuts.io
import mapmuts.sequtils


# global variable _pylabavailable indicates if pylab/matplotlib present
try:
    import matplotlib
    matplotlib.use('pdf')
    import pylab
    _pylabavailable = True
except ImportError:
    _pylabavailable = False


# global variable _scipyavailable indicates if scipy is available
try:
    import scipy
    _scipyavailable = True
except ImportError:
    _scipyavailable = False


def PylabAvailable():
    """Returns `True` if ``pylab`` / ``matplotlib`` available, `False` otherwise.
    
    You should call this function to test for the availability of the
    ``pylab`` and ``matplotlib`` plotting modules before using other functions in
    this module.
    """
    return _pylabavailable


def PlotAlignmentStatistics(infiles, names, plotfile):
    """Makes plot of alignment statistics for one or more samples.

    `infiles` : a list of file name(s) containing the alignment
    statistics. These should be in the format written by
    `mapmuts.io.WriteAlignmentStatistics`.

    `names` : a list of strings giving the names of the samples
    corresponding to each entry in `infiles`.

    `plotfile` : the name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    This function uses ``pylab`` / ``matplotlib`` to plot summary statistics
    showing how the reads fall into different categories based on
    whether they could be aligned. The data is taken from the files
    in `infiles`. If there is just one infile (one sample), then a pie
    chart is plotted. If there are multiple infiles (multiple
    samples) then a stacked bar graph is plotted. The categories of
    reads are:

        `filtered` : equal to `nfiltered`

        `low Q` : equal to `nlowq`

        `excess N` : equal to `nexcessn`

        `unpaired` : equal to `nattempted - npaired`

        `unaligned` : equal to `npaired - nalignedfullgene`

        `outside gene` : equal to `nalignedfullgene - nalignedgene`

        `aligned` : equal to `nalignedgene`

    This function requires ``pylab`` / ``matplotlib``. It will raise an exception
    if these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    if not (isinstance(infiles, list)) and infiles:
        raise ValueError("infiles does not specify a non-empty list")
    if len(infiles) != len(names):
        raise ValueError("infiles and names differ in length")
    names = [name.replace('_', '\_') for name in names]
    keys = ['aligned', 'outside gene', 'unaligned',\
            'unpaired', 'excess N', 'low Q', 'filtered']
    d = dict([(key, []) for key in keys])
    for infile in infiles:
        if not os.path.isfile(infile):
            raise IOError("Cannot find file:\n%s" % infile)
        s = mapmuts.io.ReadAlignmentStatistics(open(infile))
        d['aligned'].append(s['nalignedgene'])
        d['outside gene'].append(s['nalignedfullgene'] -\
                s['nalignedgene'])
        d['unaligned'].append(s['npaired'] -\
                s['nalignedfullgene'])
        d['unpaired'].append(s['nattempted'] - s['npaired'])
        d['excess N'].append(s['nexcessn'])
        d['low Q'].append(s['nlowq'])
        d['filtered'].append(s['nfiltered'])
    nsamples = len(names)
    colors = list('bcmgwyr')
    assert len(colors) == len(keys)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('legend', fontsize=10)
    matplotlib.rc('font', size=10)
    matplotlib.rc('patch', linewidth=0.75)
    if nsamples == 1:
        # make a pie chart
        fig = pylab.figure(figsize=(4, 2))
        ax = pylab.axes([0.46, -0.08, 0.58, 1.16])
        areas = [d[key][0] for key in keys]
        (wedges, texts) = pylab.pie(areas, colors=colors)
        labels = []
        for (key, area) in zip(keys, areas):
            labels.append('%s ($%s$)' % (key, Base10Formatter(area, 3, 1,\
                    1)))
        pylab.legend(wedges, labels, loc='center left', ncol=1,\
                handlelength=1.3, columnspacing=1, labelspacing=0.8,\
                handletextpad=0.7,\
                title='{\\bf \underline{number of read pairs}}',\
                bbox_to_anchor=(-0.78, 0.5), frameon=False)
    elif nsamples > 1:
        # make stacked bar graph
        fig = pylab.figure(figsize=(6, 4.2))
        (lmargin, rmargin, bmargin, tmargin) = (0.08, 0.01, 0.41, 0.12)
        ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
                bmargin - tmargin])
        bottoms = [0] * nsamples
        bars = []
        barwidth = 0.5
        indices = [i - barwidth / 2. for i in range(1, nsamples + 1)]
        for (key, color) in zip(keys, colors):
            b = pylab.bar(indices, d[key], width=barwidth, bottom=bottoms,\
                color=color)
            bars.append(b)
            for i in range(nsamples):
                bottoms[i] += d[key][i]
        ymax = max(bottoms)
        pylab.gca().set_ylim([0, 1.04 * ymax])
        pylab.xticks([i + barwidth / 2. for i in indices], names,\
                rotation=90)
        pylab.gca().set_xlim([0.5, nsamples + 0.5])
        yformatter = pylab.ScalarFormatter(useMathText=True)
        yformatter.set_powerlimits((-3, 3))
        pylab.gca().yaxis.set_major_formatter(yformatter)
        pylab.legend([b[0] for b in bars], keys, handlelength=1.4,\
                columnspacing=1.1, bbox_to_anchor=(0.53, 1.26), loc=\
                'upper center', ncol=4)
        pylab.ylabel('number of read pairs')
    else:
        raise ValueError("No samples specified.")
    if plotfile:
        pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotInsertLengths(infile, plotfile, normed=True):
    """Plots distribution of insert lengths.

    `infile` : name of input file containing the distribution of insert
    lengths, in the format written by `mapmuts.io.WriteInsertLengths`.

    `plotfile` : name of the output plot file created by this method
    (such as 'plot.pdf'). Must end in the extension ``.pdf``.

    `normed` : Boolean switch specifying that we plot the fractional
    distribution of insert lengths (so that all bars sum to one).
    This is done if this switch has its default value of `True`. If
    set to `False`, instead plot the actual number of reads of each
    length.

    This function uses ``pylab`` / ``matplotlib`` to make a histogram showing
    the distribution of insert lengths.

    This function requires ``pylab`` / ``matplotlib``. It will raise an exception
    if these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    if not os.path.isfile(infile):
        raise IOError("Cannot find infile:\n%s" % infile)
    f = open(infile)
    d = mapmuts.io.ReadInsertLengths(f)
    f.close()
    if not d: 
        raise ValueError("Empty data for insert lengths.")
    (minl, maxl) = (min(d.keys()), max(d.keys()))
    lengths = []
    counts = []
    sum = 0
    for i in range(minl, maxl + 1):
        lengths.append(i)
        if i in d:
            counts.append(d[i])
            sum += d[i]
        else:
            counts.append(0)
    if not sum:
        raise ValueError("No insert lengths specified.")
    if normed:
        counts = [x / float(sum) for x in counts]
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('legend', fontsize=10)
    matplotlib.rc('font', size=10)
    matplotlib.rc('patch', linewidth=0.6)
    figure = pylab.figure(figsize=(4, 2))
    if normed:
        tmargin = 0.035
    else:
        tmargin = 0.09
    (lmargin, rmargin, bmargin, tmargin) = (0.14, 0.01, 0.19, tmargin)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    barwidth = 1.0
    pylab.bar([x - barwidth / 2.0 for x in lengths], counts,\
            width=barwidth)
    yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-3, 3))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    pylab.gca().set_xlim([minl - 1, maxl + 1])
    pylab.xlabel('insert length')
    if normed:
        pylab.ylabel('fraction of reads')
    else:
        pylab.ylabel('number of reads')
    if plotfile:
        pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotReadMismatches(infile, plotfile, logy, title=None):
    """Plots distribution of read mismatches.

    `infile` : name of input file containing the distribution of
    mismatches along the length of a read, in the format written
    by `mapmuts.io.WriteReadMismatches`.

    `plotfile` : name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    `logy` : Boolean switch specifying that we use a log scale for the
    y-axis (fraction of mismatches). We use a log scale if and only
    if this switch is `True`.

    `title` : string specifying a title placed above the plot. Can be
    set to `None` if you do not want a title.

    This function uses ``pylab`` / ``matplotlib`` to plot the distribution
    of mismatches along the length of the gene. The fraction of reads with
    a mutation at a given site are shown in 1, 2, ... numbering from the
    start of the read.

    This function requires ``pylab`` / ``matplotlib``. It will raise an exception
    if these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if not os.path.isfile(infile):
        raise IOError("Cannot find infile:\n%s" % infile)
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError('plotfile must end in .pdf: %s' % plotfile)
    f = open(infile)
    d = mapmuts.io.ReadReadMismatches(f)
    f.close()
    maxi = max(d.keys())
    positions = []
    fracmm = []
    for i in range(1, maxi + 1):
        positions.append(i)
        try:
            fracmm.append(d[i])
        except KeyError:
            fracmm.append(0.0)
    pylab.figure(figsize=(3, 2.25))
    (lmargin, rmargin, bmargin, tmargin) = (0.19, 0.01, 0.17, 0.1)
    pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 - tmargin -\
            bmargin])
    pylab.plot(positions, fracmm, 'd', markersize=5.5)
    if logy:
        pylab.gca().set_yscale('log')
    else:
        yticker = matplotlib.ticker.MaxNLocator(5)
        pylab.gca().yaxis.set_major_locator(yticker)
    pylab.gca().set_xlim([0.1, maxi + 0.9])
    pylab.xlabel('read position')
    pylab.ylabel('fraction with mismatch')
    if title:
        pylab.title(title.replace('_', '\_'), fontsize=10)
    if plotfile:
        pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotCodonDepth(infile, plotfile):
    """Plots codon read depth as a function of gene position.

    `infile` : name of input file of the type written by 
    `mapmuts.io.WriteCodonCounts`.

    `plotfile` : name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    This function uses ``pylab`` / ``matplotlib`` to plot the read depth
    as a function of codon position along a gene. Codon
    positions are numbered 1, 2, ... from the start of the gene.
    Plots number of times each codon is called by both paired,
    a single paired read, and is called to an ambiguous identity.

    This function requires ``pylab`` / ``matplotlib``. It will raise an exception
    these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if not os.path.isfile(infile):
        raise IOError("Cannot find infile:\n%s" % infile)
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError('plotfile must end in .pdf: %s' % plotfile)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)
    matplotlib.rc('legend', fontsize=10)
    codon_counts = mapmuts.io.ReadCodonCounts(open(infile))
    mapmuts.sequtils.ClassifyCodonCounts(codon_counts)
    maxcodon = max([i for i in codon_counts.keys() if isinstance(i, int)])
    xs = [i for i in range(1, maxcodon + 1)]
    ys = {}
    keys = ['PAIRED_COUNTS', 'SINGLE_COUNTS', 'AMBIGUOUS_COUNTS']
    labels = ['paired reads', 'single read', 'ambiguous']
    styles = {
              'PAIRED_COUNTS':'b-',
              'SINGLE_COUNTS':'r-',
              'AMBIGUOUS_COUNTS':'g-',
             }
    for key in keys:
        ys[key] = [codon_counts[i][key] for i in xs]
    fig = pylab.figure(figsize=(5.5, 2.4))
    (lmargin, rmargin, bmargin, tmargin) = (0.09, 0.02, 0.16, 0.11)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    lines = []
    for key in keys:
        plotline = pylab.plot(xs, ys[key], styles[key], lw=1.2)
        lines.append(plotline[0])
    pylab.xlabel('codon number')
    pylab.ylabel('number of reads')
    yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    pylab.gca().set_xlim([1, max(xs)])
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-3, 3))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    pylab.legend(lines, labels, handlelength=2,\
            bbox_to_anchor=(0.53, 1.15), loc=\
            'upper center', ncol=3)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotCodonDist(infile, plotfile, dist_type):
    """Plots distribution of types of codon mutations across gene.

    Takes as input a file (*infile*) specifying the number of counts
    for each codon along the gene. Creates a PDF plot (*plotfile*)
    showing the frequency of different types of codon mutations
    along the gene, as specified by *dist_type*. Only codon counts
    for which all three nucleotides are definitively called (upper
    case nucleotides in *infile*) are counted for the plots.

    *infile* : name of input file of the type written by 
    *mapmuts.io.WriteCodonCounts*. This file specifies the counts
    for all of the codon mutations at each position.

    *plotfile* : name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    *dist_type* : the type of distribution that is plotted along
    the gene. It is a string with the following possible values:
    
        * *syn_ns_stop* specifies that we plot the distribution
          of synonymous, nonsynonymous, and stop codons.

        * *codon_nmuts* specifies that we plot the distribution
          of the number of nucleotide mutations per codon for 
          all codon mutations (1, 2, or 3).

    This function requires ``pylab`` / ``matplotlib``. It will raise an exception
    these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if not os.path.isfile(infile):
        raise IOError("Cannot find infile:\n%s" % infile)
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError('plotfile must end in .pdf: %s' % plotfile)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)
    matplotlib.rc('legend', fontsize=10)
    codon_counts = mapmuts.io.ReadCodonCounts(open(infile))
    mapmuts.sequtils.ClassifyCodonCounts(codon_counts)
    maxcodon = max([i for i in codon_counts.keys() if isinstance(i, int)])
    xs = [i for i in range(1, maxcodon + 1)]
    ys = {}
    denom = float(codon_counts['TOTAL_COUNTS'])
    if not denom:
        raise ValueError('no counts')
    if dist_type == 'syn_ns_stop':
        keys = ['nonsynonymous', 'synonymous', 'stop codon']
        fracs = [codon_counts['TOTAL_NS'] / denom, 
                codon_counts['TOTAL_SYN'] / denom,
                codon_counts['TOTAL_STOP'] / denom]
        key_d = {'nonsynonymous':'N_NS', 'synonymous':'N_SYN', \
                'stop codon':'N_STOP'}
        styles = {
                'nonsynonymous':'b-',
                'synonymous':'r-',
                'stop codon':'g-',
                }
    elif dist_type == 'codon_nmuts':
        keys = ['1 nucleotide mutation', '2 nucleotide mutations', '3 nucleotide mutations']
        fracs = [codon_counts['TOTAL_N_1MUT'] / denom,
                 codon_counts['TOTAL_N_2MUT'] / denom,
                 codon_counts['TOTAL_N_3MUT'] / denom]
        key_d = {'1 nucleotide mutation':'N_1MUT',
                 '2 nucleotide mutations':'N_2MUT',
                 '3 nucleotide mutations':'N_3MUT'}
        styles = {'1 nucleotide mutation':'b-',
                  '2 nucleotide mutations':'r-',
                  '3 nucleotide mutations':'g-'}
    else:
        raise ValueError("Invalid dist_type of %s" % dist_type)
    for key in keys:
        ys[key] = []
        for i in xs:
            denom = codon_counts[i]['PAIRED_COUNTS']
            if not denom:
                ys[key].append(0)
            else:
                numerator = codon_counts[i][key_d[key]]
                ys[key].append(numerator / float(denom))
    fig = pylab.figure(figsize=(5.75, 3.1))
    (lmargin, rmargin, bmargin, tmargin) = (0.09, 0.01, 0.13, 0.23)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    lines = []
    for key in keys:
        plotline = pylab.plot(xs, ys[key], styles[key], lw=1.2)
        lines.append(plotline[0])
    pylab.xlabel('codon number')
    pylab.ylabel('fraction')
    yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    pylab.gca().set_xlim([1, max(xs)])
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-2, 2))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    labels = ["%s (overall fraction $%s$)" % (key, Base10Formatter(frac, 2, 1, 1)) for (key, frac) in zip(keys, fracs)]
    pylab.legend(lines, labels, handlelength=1.6,\
            bbox_to_anchor=(0.53, 1.36), loc=\
            'upper center', ncol=1)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotPairedMutFracs(infiles, names, plotfile):
    """Makes paired bar graph of mutation fractions per codon.
    
    For each sample, calculates the fraction of all codon counts
    (called by both reads) that represent synonymous / nonsynoymous /
    stop codons, and one- / two- / three-nucleotide mutations. Makes
    paired stacked bar graphs showing these fractions. The fractions
    for all sites in the gene are aggregated together. 

    `infiles` : a list of file names containing the codon counts for
    each sample, in the format written by `mapmuts.io.WriteCodonCounts`.
    Note that there is also an option for each entry in this list to
    be a dictionary rather than a file -- in this case, the dictionary 
    should already contain the information that would normally be
    read from the file. Specifically, it must contain a key
    for each of the strings in *keys* giving the fraction of sites
    with that type of mutation as a value.

    `names` : a list of strings giving the names of the samples
    corresponding to each entry in `infiles`.

    `plotfile` : name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    Requires ``pylab`` / ``matplotlib``, will raise an exception
    if these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    names = [name.replace('_', '\_') for name in names]
    if not (isinstance(infiles, list)) and infiles:
        raise ValueError("infiles does not specify a non-empty list")
    if len(infiles) != len(names):
        raise ValueError("infiles and names differ in length")
    bar1 = ['synonymous', 'nonsynonymous', 'stop codon']
    bar2 = ['1 nucleotide mutation', '2 nucleotide mutations', '3 nucleotide mutations']
    d = dict([(key, []) for key in bar1 + bar2])
    for infile in infiles:
        if isinstance(infile, str):
            if not os.path.isfile(infile):
                raise IOError("Cannot find file:\n%s" % infile)
            codon_counts = mapmuts.io.ReadCodonCounts(open(infile))
            mapmuts.sequtils.ClassifyCodonCounts(codon_counts)
            maxcodon = max([i for i in codon_counts.keys() if isinstance(i, int)])
            denom = float(codon_counts['TOTAL_COUNTS'])
            if not denom:
                raise ValueError("no counts for %s" % infile)
            for key in bar1 + bar2:
                if key == 'nonsynonymous':
                    d[key].append(codon_counts['TOTAL_NS'] / denom)
                elif key == 'synonymous':
                    d[key].append(codon_counts['TOTAL_SYN'] / denom)
                elif key == 'stop codon':
                    d[key].append(codon_counts['TOTAL_STOP'] / denom)
                elif key == '1 nucleotide mutation':
                    d[key].append(codon_counts['TOTAL_N_1MUT'] / denom)
                elif key == '2 nucleotide mutations':
                    d[key].append(codon_counts['TOTAL_N_2MUT'] / denom)
                elif key == '3 nucleotide mutations':
                    d[key].append(codon_counts['TOTAL_N_3MUT'] / denom)
                else:
                    raise ValueError("Invalid key of %s" % key)
        elif isinstance(infile, dict):
            for key in bar1 + bar2:
                if key not in infile:
                    raise ValueError("Did not find key for %s" % key)
                d[key].append(infile[key])
        else:
            raise ValueError("infile is not a string or dict")
    nsamples = len(names)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)
    matplotlib.rc('legend', fontsize=10)
    matplotlib.rc('xtick', labelsize=10)
    matplotlib.rc('patch', linewidth=0.5)
    # make stacked bar graph
    fig = pylab.figure(figsize=(6.5, 3.75))
    (lmargin, rmargin, bmargin, tmargin) = (0.06, 0.01, 0.43, 0.13)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    bars = []
    for (ibar, keys, colors) in [(0, bar1, 'brg'), (1, bar2, 'myc')]:
        bottoms = [0] * nsamples
        barwidth = 0.35
        indices = [i - barwidth + ibar * barwidth for i in range(1, nsamples + 1)]
        totalheights = [0 for i in range(nsamples)]
        for (key, color) in zip(keys, colors):
            totalheights = [totalheights[i] + d[key][i] for i in range(nsamples)]
            b = pylab.bar(indices, d[key], width=barwidth, bottom=bottoms, color=color)
            bars.append(b)
            for i in range(nsamples):
                bottoms[i] += d[key][i]
    ymax = max(bottoms)
    pylab.gca().set_ylim([0, 1.08 * ymax])
    pylab.xticks([i for i in indices], names, rotation=90)
    pylab.gca().set_xlim([0.4, nsamples + 0.6])
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-3, 3))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    barmarks = [b[0] for b in bars]
    barlabels = bar1 + bar2
    # reorder bar labels since pylab puts them down columns
    barmarks = [barmarks[0], barmarks[3], barmarks[1], barmarks[4], barmarks[2], barmarks[5]]
    barlabels = [barlabels[0], barlabels[3], barlabels[1], barlabels[4], barlabels[2], barlabels[5]]
    pylab.legend(barmarks, barlabels, handlelength=1.2,\
            bbox_to_anchor=(0.54, 1.3), loc='upper center', ncol=3)
    pylab.ylabel('fraction', size=10)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotMutFracs(infiles, names, plotfile, keys, writefracs, formatting={}, textwritefracs=None):
    """Makes plot of mutation fractions per codon for several samples.

    For each sample, calculates the fraction of all codon counts
    (called by both reads) that represent the types of mutations 
    specified by *keys*. Makes a stacked bar graph showing these
    counts. The fractions are for all sites in the gene
    aggregated together.

    `infiles` : a list of file names containing the codon counts for
    each sample, in the format written by `mapmuts.io.WriteCodonCounts`.
    Note that there is also an option for each entry in this list to
    be a dictionary rather than a file -- in this case, the dictionary 
    should already contain the information that would normally be
    read from the file. Specifically, it must contain a key
    for each of the strings in *keys* giving the fraction of sites
    with that type of mutation as a value.

    `names` : a list of strings giving the names of the samples
    corresponding to each entry in `infiles`.

    `plotfile` : name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    `keys` : a list specifying the order (from bottom to top) that we
    stack the bars. It can contain any of the following entries:
   
        * 'synonymous' : fraction of synonymous mutations.

        * 'nonsynonymous' : fraction of nonsynonymous mutations.

        * 'stop codon' : fraction of stop codon mutations.

        * '1 nucleotide mutation' : fraction of mutations with one
          nucleotide change.

        * '2 nucleotide mutations' : fraction of mutations with two
          nucleotide changes.

        * '3 nucleotide mutations' : fraction of mutations with three
          nucleotide changes.

    Most typically you would want to set *keys* to either
    *['synonymous', 'nonsynonymous', 'stop codon']* or to
    *['1 nucleotide mutation', '2 nucleotide mutations', '3 nucleotide mutations']*.

    `writefracs` is a Boolean switch specifying whether we report the
    numerical fraction of mutations in each category above the bar. If
    *True* then we report these fractions; if *False* then we do not.

    *formatting* is a dictionary in which you can use string keys and their
    values to specify additional formatting. Specifically:

        * *figwidth* specifies the width of the figure to be this
          in inches rather than the function's default width of 5.75.

    * *textwritefracs* is an option that allows you to write the fraction
      of mutations in each category to a file in text format. By default,
      this option is *None*, meaning that no such writing is done. If you
      want to do such writing, set it to a string representing a filename.
      This file is then created and will list the fractions of mutations
      in each category in *keys* for each of the samples in *names*.

    This function requires ``pylab`` / ``matplotlib``. It will raise an exception
    if these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    names = [name.replace('_', '\_') for name in names]
    if not (isinstance(infiles, list)) and infiles:
        raise ValueError("infiles does not specify a non-empty list")
    if len(infiles) != len(names):
        raise ValueError("infiles and names differ in length")
    d = dict([(key, []) for key in keys])
    for infile in infiles:
        if isinstance(infile, str):
            if not os.path.isfile(infile):
                raise IOError("Cannot find file:\n%s" % infile)
            codon_counts = mapmuts.io.ReadCodonCounts(open(infile))
            mapmuts.sequtils.ClassifyCodonCounts(codon_counts)
            maxcodon = max([i for i in codon_counts.keys() if isinstance(i, int)])
            denom = float(codon_counts['TOTAL_COUNTS'])
            if not denom:
                raise ValueError("no counts for %s" % infile)
            for key in keys:
                if key == 'nonsynonymous':
                    d[key].append(codon_counts['TOTAL_NS'] / denom)
                elif key == 'synonymous':
                    d[key].append(codon_counts['TOTAL_SYN'] / denom)
                elif key == 'stop codon':
                    d[key].append(codon_counts['TOTAL_STOP'] / denom)
                elif key == '1 nucleotide mutation':
                    d[key].append(codon_counts['TOTAL_N_1MUT'] / denom)
                elif key == '2 nucleotide mutations':
                    d[key].append(codon_counts['TOTAL_N_2MUT'] / denom)
                elif key == '3 nucleotide mutations':
                    d[key].append(codon_counts['TOTAL_N_3MUT'] / denom)
                else:
                    raise ValueError("Invalid key of %s" % key)
        elif isinstance(infile, dict):
            for key in keys:
                if key not in infile:
                    raise ValueError("Did not find key for %s" % key)
                d[key].append(infile[key])
        else:
            raise ValueError("infile is not a string or dict")
    nsamples = len(names)
    colors = list('brgmyk')
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)
    matplotlib.rc('legend', fontsize=11)
    matplotlib.rc('xtick', labelsize=11)
    matplotlib.rc('patch', linewidth=0.75)
    # make stacked bar graph
    ncol = 2
    if 'figwidth' in formatting:
        figwidth = formatting['figwidth']
        if figwidth > 7:
            ncol = 3
        if writefracs:
            legendy = 1.16
        else:
            legendy = 1.18
    else:
        figwidth = 5.75
        if writefracs:
            legendy = 1.24
        else:
            legendy = 1.3
    if writefracs:
        fig = pylab.figure(figsize=(figwidth, 4.5))
        (lmargin, rmargin, bmargin, tmargin) = (0.1, 0.01, 0.41, 0.11)
    else:
        fig = pylab.figure(figsize=(figwidth, 4.0))
        (lmargin, rmargin, bmargin, tmargin) = (0.1, 0.01, 0.43, 0.13)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    bottoms = [0] * nsamples
    bars = []
    barwidth = 0.5
    indices = [i - barwidth / 2. for i in range(1, nsamples + 1)]
    totalheights = [0 for i in range(nsamples)]
    for (key, color) in zip(keys, colors):
        totalheights = [totalheights[i] + d[key][i] for i in range(nsamples)]
        b = pylab.bar(indices, d[key], width=barwidth, bottom=bottoms,\
            color=color)
        bars.append(b)
        for i in range(nsamples):
            bottoms[i] += d[key][i]
    ymax = max(bottoms)
    if writefracs:
        pylab.gca().set_ylim([0, 1.25 * ymax])
        yspace = 0.02 * ymax
        for (key, color) in zip(keys, colors):
            for (index, frac, totalheight) in zip(indices, d[key], totalheights):
                x = index + barwidth / 2.0
                text = Base10Formatter(frac, 2, 1, 2)
                text = '{\setlength{\medmuskip}{0mu} $%s$}' % text
                pylab.text(x, totalheight + yspace, text, color=color, horizontalalignment='center', size=7)
            yspace += 0.07 * ymax
    else:
        pylab.gca().set_ylim([0, 1.08 * ymax])
    if textwritefracs:
        f = open(textwritefracs, 'w')
        f.write("#SAMPLE\tFRAC_%s\tFRAC_%s\tFRAC_%s\n" % tuple(keys))
        i = 0
        for name in names:
            f.write("%s" % name)
            for key in keys:
                f.write("\t%.5G" % (d[key][i]))
            f.write('\n')
            i += 1
        f.close()
    pylab.xticks([i + barwidth / 2. for i in indices], names,\
            rotation=90)
    pylab.gca().set_xlim([0.4, nsamples + 0.6])
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-3, 3))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    pylab.legend([b[0] for b in bars], keys, handlelength=1.4,\
            bbox_to_anchor=(0.53, legendy), loc='upper center', ncol=ncol)
    pylab.ylabel('fraction', size=11)
    if plotfile:
        pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotNTMutFracs(infiles, names, plotfile):
    """Makes plot of nucleotide mutation types for several samples.

    For each sample, calculates the fraction of each nucleotide call
    (called by both reads) that represents each type of nucleotide
    substitution (i.e. C->G, G->C; A->G, C->T); etc). Makes a stacked
    bar graph showing the fraction of sites with these mutations
    over aggregated over all positions.

    `infiles` : a list of file names containing the nucleotide counts for
    each sample, in the format written by `mapmuts.io.WriteNTCounts`.
    Note that there is also an option for each entry in this list to
    be a dictionary rather than a file -- in this case, the dictionary 
    should already contain the information that would normally be
    read from the file. Specifically, it must contain a key
    for each entry in `sequtils.NTMutTypes()` giving the fraction of sites
    with that type of mutation as a value.

    `names` : a list of strings giving the names of the samples
    corresponding to each entry in `infiles`.

    `plotfile` : the name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    This function uses ``pylab`` / ``matplotlib``. It will raise an exception
    if these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    if not (isinstance(infiles, list)) and infiles:
        raise ValueError("infiles does not specify a non-empty list")
    if len(infiles) != len(names):
        raise ValueError("infiles and names differ in length")
    keys = mapmuts.sequtils.NTMutTypes()
    d = dict([(key, []) for key in keys])
    names = [name.replace('_', '\_') for name in names]
    for infile in infiles:
        if isinstance(infile, str):
            if not os.path.isfile(infile):
                raise IOError("Cannot find file:\n%s" % infile)
            nt_counts = mapmuts.io.ReadNTCounts(open(infile))
            mapmuts.sequtils.ClassifyNTMuts(nt_counts)
            maxnt = max([i for i in nt_counts.keys() if isinstance(i, int)])
            denom = float(nt_counts['PAIRED_COUNTS'])
            if not denom:
                raise ValueError("no counts for %s" % infile)
            for key in keys:
                d[key].append(nt_counts[key] / denom)
        elif isinstance(infile, dict):
            for key in keys:
                if not key in infile:
                    raise ValueError("dict in infiles lacks key %s" % str(key))
                d[key].append(infile[key])
        else:
            raise ValueError("infile is not a string or dict")
    nsamples = len(names)
    colors = list('bgrcmy')
    assert len(colors) == len(keys)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('legend', fontsize=10)
    matplotlib.rc('font', size=10)
    matplotlib.rc('patch', linewidth=0.75)
    # make stacked bar graph
    fig = pylab.figure(figsize=(5.75, 4.0))
    (lmargin, rmargin, bmargin, tmargin) = (0.09, 0.01, 0.43, 0.12)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    bottoms = [0] * nsamples
    bars = []
    barwidth = 0.5
    indices = [i - barwidth / 2. for i in range(1, nsamples + 1)]
    for (key, color) in zip(keys, colors):
        b = pylab.bar(indices, d[key], width=barwidth, bottom=bottoms,\
            color=color)
        bars.append(b)
        for i in range(nsamples):
            bottoms[i] += d[key][i]
    ymax = max(bottoms)
    pylab.gca().set_ylim([0, 1.08 * ymax])
    pylab.xticks([i + barwidth / 2. for i in indices], names,\
            rotation=90)
    pylab.gca().set_xlim([0.5, nsamples + 0.5])
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-3, 3))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    barnames = []
    for key in keys:
        barname = "%s $\\rightarrow$ %s, %s $\\rightarrow$ %s" % \
                (key[0][0], key[0][1], key[1][0], key[1][1])
        barnames.append(barname)
    pylab.legend([b[0] for b in bars], barnames, handlelength=1.4,\
            bbox_to_anchor=(0.53, 1.28), loc='upper center', ncol=3)
    pylab.ylabel('fraction')
    if plotfile:
        pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotMutCountFracs(plotfile, title, names, all_cumulfracs, syn_cumulfracs, all_counts, syn_counts, legendloc):
    """Plots fraction of mutations with >= a given number of counts.

    Does this for both all mutations and synonymous mutations. The plots
    are placed side by side.

    The plots show the fractions of mutations
    that are found >= *n* times for a range of values of *n*.

    Requires *pylab* / *matplotlib*, and will raise an exception
    if not available for import.

    CALLING VARIABLES:

    * *plotfile* : string giving name of the created plot file. Must end in
      the extension ``.pdf``.

    * *title* : string giving the title placed about the plot.

    * *names* : a list of strings giving the names of the samples to
      plot.

    * *all_cumulfracs* : a list of the same length as *names* giving
      the cumulative fractions for all mutations. Each entry should be a list, 
      and all of these lists should be the same length. *cumulfracs[i][n]*
      gives the fraction of mutations to sample *names[i]* that
      are found >= *n* times. The x-axis of the created plot
      will go from 0 to *len(all_cumulfracs[0]) - 1*.

    * *syn_cumulfracs* : a list like *all_cumulfracs* except for synonymous
      mutations.

    * *all_counts* : integer counts of all mutations (the total number of mutations
      used for *all_cumulfracs*), used to create plot title.

    * *syn_counts* : like *all_counts* but for synonymous mutations.

    * *legendloc* : specifies the location of the legend. Should be a string.
      Valid values are:

        - *bottom* : put legend at the bottom of the plot.

        - *right* : put legend at the right of the plot.
    """

    # some basic checks 
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError('plotfile must end in .pdf: %s' % plotfile)

    # plot setup stuff
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=11)
    matplotlib.rc('legend', fontsize=11)
    if legendloc == 'bottom':
        ncol = 4 # number of columns
        legendrowheight = 0.2
        nrows = math.ceil(len(names) / float(ncol))
        xlegendmargin = 0.0
    elif legendloc == 'right':
        ncol = 1
        nrows = 0 # we specify by columns
        legendrowheight = 0
        xlegendmargin = 1.45
    else:
        raise ValueError("Invalid legendloc of %s" % legendloc)
    (xsize, ysize) = (4.75 + xlegendmargin, legendrowheight * nrows + 2.4)
    styles = ['k-.', 'y-.', 'b-', 'r:', 'g:', 'c--', 'm--']
    lmargin = 0.11 * (xsize - xlegendmargin) / xsize # left margin for plot
    rmargin = 0.01 + xlegendmargin / xsize # right margin for plot
    centermargin = 0.02 # horizontal space between plots
    tmargin = 0.16 # top margin
    bmargin = 0.22 + (legendrowheight * nrows) / ysize # bottom margin
    plotwidth = (1.0 - lmargin - rmargin - centermargin) / 2.0

    assert 0 < len(syn_cumulfracs) == len(all_cumulfracs) == len(names) <= len(styles), "Must specify equal numbers of fractions and names, and no more than %d" % len(styles)

    # start plotting
    fig = pylab.figure(figsize=(xsize, ysize))
    ax_all = pylab.axes([lmargin, bmargin, plotwidth, 1 - bmargin - tmargin])
    ax_syn = pylab.axes([lmargin + plotwidth + centermargin, bmargin, plotwidth, 1 - bmargin - tmargin])
    for (cumulfracs, ax, write_ylabel, ax_title) in [(syn_cumulfracs, ax_syn, False, 'synonymous (%d total)' % syn_counts), (all_cumulfracs, ax_all, True, 'all (%d total)' % all_counts)]:
        pylab.axes(ax)
        nmax = len(cumulfracs[0])
        assert nmax, "Length of entries in cumulfracs must be >= 1"
        for xlist in cumulfracs:
            if len(xlist) != nmax:
                raise ValueError("Not all entries in cumulfracs of the same length")
        lines = []
        xs = [n for n in range(0, nmax)]
        for i in range(len(names)):
            plotline = pylab.plot(xs, cumulfracs[i], styles[i], lw=1.5)
            lines.append(plotline[0])
            pylab.xlabel("Mutation counts", size=11)
            if write_ylabel:
                pylab.ylabel("Frac. $\ge$ this many counts", size=11)
            else:
                yformatter = matplotlib.ticker.NullFormatter()
                pylab.gca().yaxis.set_major_formatter(yformatter)
        pylab.gca().set_ylim([-0.02, 1.02])
        yticker = matplotlib.ticker.FixedLocator([0.0, 0.5, 1.0])
        pylab.gca().yaxis.set_major_locator(yticker)
        pylab.gca().set_xlim([0, nmax - 1])
        xticker = matplotlib.ticker.MaxNLocator(4)
        pylab.gca().xaxis.set_major_locator(xticker)
        pylab.title(ax_title, size=11)
    pylab.suptitle("{\\bf %s}" % title, size=11)
    if legendloc == 'bottom':
        fig.legend(lines, names, handlelength=2.25, handletextpad=0.2, columnspacing=0.8, ncol=ncol, bbox_to_anchor=(0.5, -0.01), loc='lower center')
    elif legendloc == 'right':
        fig.legend(lines, names, handlelength=2.25, handletextpad=0.2, columnspacing=0.8, ncol=ncol, bbox_to_anchor=(1.0, 0.52), loc='center right')
    else:
        raise ValueError("Invalid legendloc of %s" % legendloc)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()



def PlotTraces(traces, plotfile, xlabel, ylabel, title, trace_labels=None):
    """Plots traces from an MCMC run.

    Requires ``pylab`` and will raise an exception of *PylabAvailable()*
    is *False*.

    CALLING VARIABLES:

    * *traces* should be a list of one or more lists or arrays. All
      entries must be the same length. These should give the traces
      of the posterior to be plotted (it is up to you if you want
      to have applied burn-ins or thinning before the traces). In the
      legend, if *trace_labels* is *None*, then the different traces 
      are labeled as #1, #2, etc. Otherwise *trace_labels* should be 
      a list of strings of the same length as *traces* that gives
      the labels applied to the traces.

    * *plotfile* : name of the PDF file created by this method.
      Should have the extension ``.pdf``.

    * *xlabel* : string giving the label on the X-axis.

    * *ylabel* : string giving the label on the Y-axis.

    * *title* : string giving the title placed above the plot.

    * *trace_labels* : an optional argument that is *None* by default.
      If set to some other value, it should be a list of strings of the
      same length as *traces*. The strings are used to label the corresponding
      traces in *traces* in the legend. Furthermore, if an identical label
      appears multiple times in *trace_labels*, it is only shown on the 
      legend once. So if you have multiple runs of the same sample in
      *traces*, just give them the same label in *trace_labels* and then they 
      will be shown with the same linestyle and only one legend label.

    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError('plotfile must end in .pdf: %s' % plotfile)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)
    matplotlib.rc('legend', fontsize=10)
    ntraces = len(traces)
    colors = ['b', 'r', 'g', 'm', 'y', 'k', 'c']
    linestyles = ['-', ':', '--', '-.']
    styles = []
    for ls in linestyles:
        for c in colors:
            styles.append('%s%s' % (c, ls))
    assert ntraces, "No traces"
    xmax = len(traces[0])
    xs = [x for x in range(1, xmax + 1)]
    for itrace in range(ntraces):
        if xmax != len(traces[itrace]):
            raise ValueError("traces differ in length")
    fig = pylab.figure(figsize=(8, 4))
    (lmargin, rmargin, bmargin, tmargin) = (0.08, 0.02, 0.1, 0.27)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    lines = []
    labels = []
    if trace_labels and len(trace_labels) != ntraces:
        raise ValueError("Invalid number of trace_labels")
    label_d = {}
    for itrace in range(ntraces):
        if not trace_labels:
            labels.append("\# %d" % (itrace + 1))
            style = styles[itrace % len(styles)]
            plotline = pylab.plot(xs, traces[itrace], style)
            lines.append(plotline[0])
        else:
            label = trace_labels[itrace]
            if label in label_d:
                style = label_d[label]
                plotline = pylab.plot(xs, traces[itrace], style)
            else:
                style = styles[len(label_d) % len(styles)]
                label_d[label] = style
                plotline = pylab.plot(xs, traces[itrace], style)
                labels.append(label)
                lines.append(plotline[0])
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    pylab.gca().set_xlim([1, max(xs)])
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-3, 3))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    pylab.legend(lines, labels, handlelength=3, bbox_to_anchor=(0.5, 1.4),
            loc='upper center', ncol=7)
    pylab.title(title, size=10)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def EquilibriumFreqsHeatMap(sites, pi_d, aa_ordering, plotfile, otherprops=[], interpolation='none'):
    """Constructs a heat map show amino acid equilibrium frequencies.

    This function uses ``pylab``. It will raise an exception if
    *PylabAvailable()* is not *True*.

    The heat map has 20 rows corresponding to each of the amino
    acids, and as many columns as there are entries in *sites*
    (one for each site in the protein). The heat map is colored
    according to the equilibrium frequency of that site for
    that amino acid, as specified by *pi_d*.

    The heat map is generated using the *imshow* function
    of ``pylab``.

    CALLING VARIABLES:

    * *sites* is a list of integers giving all of the sites that
      are being included in the heat map. They must be consecutive,
      so that *sites = [i for i in range(sites[0], sites[-1] + 1)]*.

    * *pi_d* is a dictionary that has a key for every integer in
      *sites*. The value of *pi_d[isite]* is itself a dictionary,
      which has keys 'PI_A', 'PI_C', 'PI_D', etc for all 20
      one-letter upper-case amino acid codes. The values
      for these keys are the equilibrium frequencies of that
      amino acid at that site. So *pi_d[isite]['PI_M']* is
      the equilibrium frequency for methionine at site *isite*.

    * *aa_ordering* is a string specifying how the amino acids
      are ordered in the rows from top to bottom. Typically you
      might want them to roughly be ordered by physical property
      for a more appealing layout. *aa_ordering* is a string that
      can be passed to *mapmuts.sequtils.OrderedAminoAcids* to
      get an ordered list of the amino acids. For example, valid 
      values are 'icelogo' and 'hydrophobicity'.

    * *plotfile* is a string giving the name of the PDF file containing
      the heat map plot. It must end in the extension ``.pdf``.

    * *otherprops* is an optional argument that can be used to specify
      additional bars that are plotted above the main heat map
      of the site preferences. By default, *otherprops* is an empty
      list, which means that no additional bars are plotted. If you
      want to plot additonal bars, then add tuple entries to the list.
      The format of these tuples differs depending on whether the
      data is scalar (such as RSA or site entropy) or categorical
      (such as secondary structure):
      
        - Scalar data: each tuple should be *(propname, prop_d)* where *propname* 
          is a string giving the name of the property and *prop_d* is a dictionary
          keyed by the site number and with *prop_d[isite]* being the numerical
          value of property *propname* for site *isite* for *isite* in *sites*.

        - Categorical data: each tuple should be *(propname, prop_d, categories)*
          where *propname* is a string giving the name of the property,
          *prop_d* is a dictionary keyed by the site number and with *prop_d[isite]*
          being the value of property *propname* for site *isite* for *isite*
          in *sites, and with *categories* being a listing of all of the categories
          of values found in *prop_d*.

      Unlike for *pi_d*, it is not required that *prop_d* have a key for 
      all values of *isite*. Sites for which no value is specified
      are shown in white.

      Currently you are only allowed to specify a maximum of three
      entries in *otherprops*.

    * *interpolation* specifies how the colors are interpolated between
      points on the heat map. You can pass any value that is acceptable
      for the parameter of the same name in the ``pylab`` *imshow*
      function. The default value is 'none' which corresponds to
      no interpolation.
    """
    if not _pylabavailable:
        raise ValueError("Cannot import pylab")
    if sites != [i for i in range(sites[0], sites[-1] + 1)]:
        raise ValueError("sites does not specify consecutive numbers")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    # maxcolorbars is the maximum number of total maps
    maxcolorbars = 4
    if len(otherprops) + 1 > maxcolorbars:
        raise ValueError("Too many otherprops specified: only %d are allowed" % (maxcolorbars - 1))
    # some general properties of the plot
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('xtick', labelsize=8)
    matplotlib.rc('xtick', direction='out')
    matplotlib.rc('ytick', direction='out')
    matplotlib.rc('axes', linewidth=0.5)
    # some dimensions in inches
    plotwidth = 6.0 # width of plot
    rowheight = 0.09 # height of rows
    lmargin = 0.5 # left margin
    rmargin = 0.05 # right margin
    tmargin = 0.2 # top margin
    bmargin = 0.4 # bottom margin
    betweenbars = 0.05 # distance between bars
    colorbarmargin = 0.35 # margin between top of plot and color bars
    colorbarheight = 0.09 # total height of color bar axes
    colorbarspacing = 0.25 # horizontal spacing between color bars    
    # set up the figure and axes
    nbars = len(otherprops) # number of additional bars
    naas = 20 # number of amino acids
    figwidth = float(plotwidth + rmargin + lmargin)
    colorbarwidth = (figwidth - maxcolorbars * colorbarspacing) / float(maxcolorbars)
    figheight = float(naas * rowheight + bmargin + tmargin + nbars * (betweenbars + rowheight) + colorbarheight + colorbarmargin)  
    fig = pylab.figure(figsize=(figwidth, figheight))
    heatmap_ax = pylab.axes([lmargin / figwidth, bmargin / figheight, 1 - (lmargin + rmargin) / figwidth, 1 - (bmargin + tmargin + nbars * (betweenbars + rowheight) + colorbarmargin + colorbarheight) / figheight])
    otherprops_ax = {}
    otherprops_image = {}
    iprop = 0
    for tup in otherprops:
        propname = tup[0]
        iprop += 1
        otherprops_ax[propname] = pylab.axes([lmargin / figwidth, (bmargin + naas * rowheight + iprop * betweenbars + (iprop - 1) * rowheight) / figheight, 1 - (lmargin + rmargin) / figwidth, rowheight / figheight])
    colorbar_axs = []
    leftpad = (maxcolorbars - len(otherprops) - 1) * colorbarwidth * 0.5
    for ibar in range(maxcolorbars):
        colorbar_ax = pylab.axes([(leftpad + 0.5 * colorbarspacing + ibar * (colorbarspacing + colorbarwidth)) / figwidth, 1.0 - (tmargin + colorbarheight) / figheight, colorbarwidth / figwidth, colorbarheight / figheight], frameon=False)
        colorbar_ax.xaxis.set_ticks_position('none')
        colorbar_ax.yaxis.set_ticks_position('none')
        pylab.xticks([])
        pylab.yticks([])
        colorbar_axs.append(colorbar_ax)
    # now make the main heat map plot after organizing the data
    aas = mapmuts.sequtils.OrderedAminoAcids(aa_ordering)
    assert len(aas) == naas
    data = pylab.empty(shape=(len(aas), len(sites)))
    iaa = 0
    for aa in aas:
        isite = 0
        for site in sites:
            data[(iaa, isite)] = pi_d[site]['PI_%s' % aa]
            isite += 1
        iaa += 1
    pylab.axes(heatmap_ax)
    heatmap_image = pylab.imshow(data, interpolation=interpolation, aspect='auto', extent=[sites[0], sites[-1], len(aas) - 0.5, -0.5], cmap=pylab.get_cmap(None))
    pylab.yticks([iaa for iaa in range(len(aas))], aas, size=7, horizontalalignment='center')
    (xticklocs, xticklabels) = pylab.xticks()
    (xmin, xmax) = pylab.xlim()
    pylab.xlabel('residue number', size=9)
    pylab.ylabel('amino acid preference', size=9)
    heatmap_ax.xaxis.set_ticks_position('bottom') # use 'none' if no ticks desired
    heatmap_ax.yaxis.set_ticks_position('left') # use 'none' if no ticks desired
    # now make the plots for otherprops
    for tup in otherprops:
        if len(tup) == 2:
            proptype = 'scalar'
        elif len(tup) == 3:
            proptype = 'categorical'
            categories = tup[2]
        else:
            raise ValueError("Entry in otherprops is tuple of invalid size")
        (propname, prop_d) = (tup[0], tup[1])
        pylab.axes(otherprops_ax[propname])
        propdata = pylab.zeros(shape=(1, len(sites)))
        propdata[ : ] = pylab.nan # set to nan for all entries
        isite = 0
        for site in sites:
            if site in prop_d:
                if proptype == 'scalar':
                    propdata[(0, isite)] = prop_d[site]
                elif proptype == 'categorical':
                    if prop_d[site] not in categories:
                        raise ValueError("Invalid categorical property: %s" % str(prop_d[site]))
                    propdata[(0, isite)] = categories.index(prop_d[site])
                else:
                    raise ValueError("proptype unknown")
            isite += 1
        otherprops_image[propname] = pylab.imshow(propdata, interpolation=interpolation, aspect='auto', extent=[sites[0], sites[-1], 0.5, -0.5], cmap=pylab.get_cmap(None))
        pylab.xticks(xticklocs, ['' for x in xticklocs])
        pylab.xlim((xmin, xmax))
        pylab.yticks([0], [propname], size=8)
        otherprops_ax[propname].xaxis.set_ticks_position('bottom') 
        otherprops_ax[propname].yaxis.set_ticks_position('left') 
    # set the color bars
    for (colorbar_ax, (propname, image, proptup)) in zip(colorbar_axs, [('amino acid preference', heatmap_image, (None, None))] + [(tup[0], otherprops_image[tup[0]], tup) for tup in otherprops]):
        pylab.axes(colorbar_ax)
        colorbar_ax.xaxis.set_ticks_position('bottom')
        pylab.title(propname, size=9)
        if len(proptup) == 2: # scalar property
            ticker = matplotlib.ticker.MaxNLocator(4)
            cb = pylab.colorbar(image, cax=colorbar_ax, orientation='horizontal', ticks=ticker)
            if propname in ['RSA']:
                cb.set_ticks([0, 0.5, 1])
                cb.set_ticklabels(['0', '0.5', '1'])
        elif len(proptup) == 3: # categorical property
            categories = proptup[2]
            cb = pylab.colorbar(image, cax=colorbar_ax, orientation='horizontal', boundaries=[i for i in range(len(categories) + 1)], values=[i for i in range(len(categories))])
            cb.set_ticks([i + 0.5 for i in range(len(categories))])
            cb.set_ticklabels(categories)
    # save the plot
    pylab.savefig(plotfile)



def EquilibriumFreqsLogo(sites, pi_d, aa_ordering, plotfile, otherprops=[], plotwidth=6.25, logoheight=0.7, nperline=71, tickevery=5):
    """Constructs a sequence logo showing amino acid equilibrium frequencies.

    THIS FUNCTION DOES NOT WORK. YOU WILL GET AN EXCEPTION IF YOU TRY TO RUN IT.

    This function uses ``pylab``. It will raise an exception if
    *PylabAvailable()* is not *True*.

    CALLING VARIABLES:

    * *sites* is a list of integers giving all of the sites that
      are being included in the heat map. They must be consecutive,
      so that *sites = [i for i in range(sites[0], sites[-1] + 1)]*.

    * *pi_d* is a dictionary that has a key for every integer in
      *sites*. The value of *pi_d[isite]* is itself a dictionary,
      which has keys 'PI_A', 'PI_C', 'PI_D', etc for all 20
      one-letter upper-case amino acid codes. The values
      for these keys are the equilibrium frequencies of that
      amino acid at that site. So *pi_d[isite]['PI_M']* is
      the equilibrium frequency for methionine at site *isite*.

    * *aa_ordering* is a string specifying how the amino acids
      are ordered in in their color choices. Typically you
      might want them to roughly be ordered by physical property
      for a more appealing layout. *aa_ordering* is a string that
      can be passed to *mapmuts.sequtils.OrderedAminoAcids* to
      get an ordered list of the amino acids. For example, valid 
      values are 'icelogo' and 'hydrophobicity'.

    * *plotfile* is a string giving the name of the PDF file containing
      the heat map plot. It must end in the extension ``.pdf``.

    * *otherprops* is an optional argument that can be used to specify
      additional bars that are plotted above the main heat map
      of the site preferences. By default, *otherprops* is an empty
      list, which means that no additional bars are plotted. If you
      want to plot additonal bars, then add tuple entries to the list.
      The format of these tuples differs depending on whether the
      data is scalar (such as RSA or site entropy) or categorical
      (such as secondary structure):
      
        - Scalar data: each tuple should be *(propname, prop_d)* where *propname* 
          is a string giving the name of the property and *prop_d* is a dictionary
          keyed by the site number and with *prop_d[isite]* being the numerical
          value of property *propname* for site *isite* for *isite* in *sites*.

        - Categorical data: each tuple should be *(propname, prop_d, categories)*
          where *propname* is a string giving the name of the property,
          *prop_d* is a dictionary keyed by the site number and with *prop_d[isite]*
          being the value of property *propname* for site *isite* for *isite*
          in *sites, and with *categories* being a listing of all of the categories
          of values found in *prop_d*.

      Unlike for *pi_d*, it is not required that *prop_d* have a key for 
      all values of *isite*. Sites for which no value is specified
      are shown in white.

      Currently you are only allowed to specify a maximum of three
      entries in *otherprops*.

    * *plotwidth* is the width of the plot in inches, is 6.25 by default.

    * *logoheight* is the height of the logo in inches, 0.7 by default.

    * *nperline* is the number of sites per line, 71 by default.

    * *tickevery* specifies where the xtick label (residue number labels) go.
      They are put at multiples of this number.
    """
    raise ValueError("FUNCTION NOT IMPLEMENTED")
    if not _pylabavailable:
        raise ValueError("Cannot import pylab")
    if sites != [i for i in range(sites[0], sites[-1] + 1)]:
        raise ValueError("sites does not specify consecutive numbers")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    # maxcolorbars is the maximum number of total maps
    maxcolorbars = 4
    if len(otherprops) + 1 > maxcolorbars:
        raise ValueError("Too many otherprops specified: only %d are allowed" % (maxcolorbars - 1))
    # some general properties of the plot
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('xtick', labelsize=8)
    matplotlib.rc('xtick', direction='out')
    matplotlib.rc('ytick', direction='out')
    matplotlib.rc('axes', linewidth=0.5)
    matplotlib.rc('ytick.major', size=3)
    matplotlib.rc('xtick.major', size=3)
    matplotlib.rc('xtick.major', pad=2.5)
    matplotlib.rc('ytick.major', pad=2.5)
    matplotlib.rc('xtick.minor', size=1.75)
    # some dimensions in inches
    lmargin = 0.4 # left margin
    rmargin = 0.1 # right margin
    tmargin = 0.1 # top margin
    bmargin = 0.0 # bottom margin
    betweenbars = 0.04 # vertical distance between bars
    barheight = 0.1 # height of color bars
    colorbarmargin = 0.35 # margin between top of plot and color bar legends
    colorbarspacing = 0.25 # horizontal spacing between color bars legends
    logobmargin = 0.2 # margin above logo
    logotmargin = betweenbars # margin below logo
    # set up the figure and axes
    sitewidth = plotwidth / float(nperline)
    nlogos = int(math.ceil(len(sites) / float(nperline)))
    nbars = len(otherprops) # number of additional property bars
    naas = 20 # number of amino acids
    figwidth = float(plotwidth + rmargin + lmargin)
    colorbarwidth = (figwidth - maxcolorbars * colorbarspacing) / float(maxcolorbars)
    figheight = float(nlogos * (logoheight + logobmargin + logotmargin + nbars * (betweenbars + barheight)) + bmargin + tmargin + barheight + colorbarmargin)  
    fig = pylab.figure(figsize=(figwidth, figheight))
    # loop over each line of the multi-lined plot
    for ilogo in range(nlogos):
        xmin = sites[ilogo * nperline] 
        xmax = min(max(sites), xmin + nperline - 1)
        xlength = (xmax - xmin + 1) * sitewidth
        logo_ax = pylab.axes([lmargin / figwidth, (bmargin + logobmargin + (nlogos - ilogo - 1) * (logoheight + logobmargin + logotmargin + nbars * (betweenbars + barheight))) / figheight, xlength / figwidth, logoheight / figheight], frameon=True)
        logo_ax.yaxis.set_ticks_position('none')
        logo_ax.xaxis.set_ticks_position('bottom')
        pylab.yticks([])
        pylab.xlim(xmin - 0.5, xmax + 0.5)
        logo_ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([x for x in range(xmin, xmax + 1) if x % tickevery == 0]))
        logo_ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([x for x in range(xmin, xmax + 1) if x % tickevery != 0]))
        for iprop in range(nbars):
            tup = otherprops[iprop]
            propname = tup[0]
            otherprops_ax = pylab.axes([lmargin / figwidth, (bmargin + (nlogos - ilogo - 1) * (logoheight + logobmargin + logotmargin + nbars * (betweenbars + barheight)) + logoheight + logobmargin + logotmargin + iprop * (betweenbars + barheight)) / figheight, xlength / figwidth, barheight / figheight], frameon=True)
            otherprops_ax.xaxis.set_ticks_position('none')
            otherprops_ax.yaxis.set_ticks_position('left')
            pylab.xticks([])
            pylab.yticks([0], [propname], size=8)
            pylab.ylim(-0.5, 0.5)
    pylab.savefig(plotfile)



def LogoOverlay(sites, plotfile, rsa_d, ss_d, nperline, sitewidth, rmargin, logoheight, barheight, barspacing):
    """Overlay for *weblogo.EquilibriumFreqsLogo* showing RSA and SS.

    This function creates colored bars showing the relative
    solvent accessibility (RSA) and secondary structure (SS) as 
    well as legends for these and residue hydrophobicity that
    can be overlayed onto the sequence logos created by
    *weblogo.EquilibriumFreqsLogo*.

    The trick of this function is to create the bars the right
    size so they align when they overlay. This function works for
    at least the cases on which it has been tested, but there
    is some chance that the overlay will not work well for
    other plots.

    This function uses ``pylab``. It will raise an exception if
    *PylabAvailable()* is not *True*.

    CALLING VARIABLES:

    * *sites* is a list of integers giving all of the sites that
      are being included in the logo. They must be consecutive,
      so that *sites = [i for i in range(sites[0], sites[-1] + 1)]*.

    * *plotfile* is a string giving the name of the PDF file containing
      the overlay. It must end in the extension ``.pdf``.

    * *otherprops* is an optional argument that can be used to specify
      additional bars that are plotted above the main heat map
      of the site preferences. By default, *otherprops* is an empty
      list, which means that no additional bars are plotted. If you
      want to plot additonal bars, then add tuple entries to the list.
      The format of these tuples differs depending on whether the
      data is scalar (such as RSA or site entropy) or categorical
      (such as secondary structure):
      
    * *rsa_d* is a dictionary giving the RSA values. It is keyed by the
      integer site numbers: *rsa_d[site]* gives a float with the RSA
      for residue *site*. If *site* is not in RSA d, the value is
      shown in white.

    * *ss_d* is like *rsa_d* except that it gives the codes for the 
      secondary structures. Allowed values are 'helix', 'strand', 
      and 'loop'.

    * *nperline* is the number of sites per line.

    * *sitewidth* is the width of each site in points.

    * *rmargin* is the right margin in points.

    * *logoheight* is the total height of each logo row in points.

    * *barheight* is the total height of each bar in points.

    * *barspacing* is the vertical spacing between bars in points.

    * *cmap* is a pylab *LinearSegmentedColorMap* used for the bars.
    """
    if not _pylabavailable:
        raise ValueError("Cannot import pylab")
    if sites != [i for i in range(sites[0], sites[-1] + 1)]:
        raise ValueError("sites does not specify consecutive numbers")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    (cmap, mapping_d, mapper) = KyteDoolittleColorMapping()
    ss_categories = ['strand', 'helix', 'loop']
    nprops = 2 # two properties, RSA and SS
    pts_per_inch = 72.0 # to convert between points and inches
    # some general properties of the plot
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('xtick', labelsize=9)
    matplotlib.rc('xtick', direction='out')
    matplotlib.rc('ytick', direction='out')
    matplotlib.rc('axes', linewidth=0.5)
    matplotlib.rc('ytick.major', size=3)
    matplotlib.rc('xtick.major', size=2.5)
    # define sizes (still in points)
    colorbar_bmargin = 20 # margin below color bars in points
    colorbar_tmargin = 15 # margin above color bars in points
    nlines = int(math.ceil(len(sites) / float(nperline)))
    lmargin = 25 # left margin in points
    barwidth = nperline * sitewidth
    figwidth = lmargin + rmargin + barwidth
    figheight = nlines * (logoheight + nprops * (barheight + barspacing)) + barheight + colorbar_bmargin + colorbar_tmargin
    # set up the figure and axes
    fig = pylab.figure(figsize=(figwidth / pts_per_inch, figheight / pts_per_inch))
    # loop over each line of the multi-lined plot
    prop_image = {}
    for iline in range(nlines):
        xmin = sites[iline * nperline] 
        xmax = min(max(sites), xmin + nperline - 1)
        isites = [x for x in range(xmin, xmax + 1)]
        assert len(isites) == xmax - xmin + 1 <= nperline
        xlength = (xmax - xmin + 1) * sitewidth
        logo_ax = pylab.axes([lmargin / figwidth, ((nlines - iline - 1) * (logoheight + nprops * (barspacing + barheight))) / figheight, xlength / figwidth, logoheight / figheight], frameon=False)
        logo_ax.yaxis.set_ticks_position('none')
        logo_ax.xaxis.set_ticks_position('none')
        pylab.yticks([])
        pylab.xlim(xmin - 0.5, xmax + 0.5)
        pylab.xticks([])
        for (iprop, propname, prop_d) in [(0, 'RSA', rsa_d), (1, 'SS', ss_d)]:
            prop_ax = pylab.axes([lmargin / figwidth, ((nlines - iline - 1) * (logoheight + nprops * (barspacing + barheight)) + logoheight + iprop * (barspacing + barheight)) / figheight, xlength / figwidth, barheight / figheight], frameon=True)
            prop_ax.xaxis.set_ticks_position('none')
            prop_ax.yaxis.set_ticks_position('left')
            pylab.xticks([])
            pylab.yticks([0], [propname], size=8)
            pylab.xlim((xmin, xmax))
            pylab.ylim(-0.5, 0.5)
            propdata = pylab.zeros(shape=(1, len(isites)))
            propdata[ : ] = pylab.nan # set to nan for all entries
            isite = 0
            for site in isites:
                if site in prop_d:
                    if propname == 'RSA':
                        propdata[(0, isite)] = prop_d[site]
                    elif propname == 'SS':
                        if prop_d[site] not in ss_categories:
                            raise ValueError("Invalid secondary structure: %s" % str(prop_d[site]))
                        propdata[(0, isite)] = ss_categories.index(prop_d[site])
                    else:
                        raise ValueError("propname unknown")
                isite += 1
            if propname == 'RSA':
                (vmin, vmax) = (0, max(prop_d.values()))
            elif propname == 'SS':
                (vmin, vmax) = (0, len(ss_categories) - 1)
            prop_image[propname] = pylab.imshow(propdata, interpolation='nearest', aspect='auto', extent=[isites[0], isites[-1], 0.5, -0.5], cmap=cmap, vmin=vmin, vmax=vmax)
            pylab.yticks([0], [propname], size=8)
    # set up colorbar axes, then color bars
    colorbarspacingfrac = 0.4 # space between color bars is this fraction of bar width
    colorbarwidth = 1.0 / ((nprops + 1) * (1.0 + colorbarspacingfrac)) # width of color bars in fraction of figure width
    colorbarspacingwidth = colorbarwidth * colorbarspacingfrac # width of color bar spacing in fraction of figure width
    propnames = {'AA':'amino-acid hydrophobicity', 'RSA':'relative solvent accessibility (RSA)', 'SS':'secondary structure (SS)'}
    for (icolorbar, propcode) in  zip(range(nprops + 1), ['AA', 'RSA', 'SS']):
#        colorbar_ax = pylab.axes([colorbarspacingwidth * 0.5 + icolorbar * (colorbarwidth + colorbarspacingwidth), 1.0 - (colorbar_tmargin + barheight) / figwidth, colorbarwidth, barheight / figwidth], frameon=True)
        colorbar_ax = pylab.axes([colorbarspacingwidth * 0.5 + icolorbar * (colorbarwidth + colorbarspacingwidth), 1.0 - (colorbar_tmargin + barheight) / figheight, colorbarwidth, barheight / figheight], frameon=True)
        colorbar_ax.xaxis.set_ticks_position('bottom')
        colorbar_ax.yaxis.set_ticks_position('none')
        pylab.xticks([])
        pylab.yticks([])
        pylab.title(propnames[propcode], size=9)
        if propcode == 'AA':
            hydrophobicities = [mapmuts.sequtils.KyteDoolittle(aa) for aa in mapmuts.sequtils.AminoAcids()]
            mapper.set_array([hydrophobicities])
            cb = pylab.colorbar(mapper, cax=colorbar_ax, orientation='horizontal')
            cb.set_ticks([-4, -2, 0, 2, 4])
            cb.set_ticklabels(['4', '2', '0', '-2', '-4'])
        elif propcode == 'RSA':
            cb = pylab.colorbar(prop_image[propcode], cax=colorbar_ax, orientation='horizontal')
            cb.set_ticks([0, 0.5, 1])
            cb.set_ticklabels(['0', '0.5', '1'])
        elif propcode == 'SS':
            cb = pylab.colorbar(prop_image[propcode], cax=colorbar_ax, orientation='horizontal', boundaries=[i for i in range(len(ss_categories) + 1)], values=[i for i in range(len(ss_categories))])
            cb.set_ticks([i + 0.5 for i in range(len(ss_categories))])
            cb.set_ticklabels(ss_categories)
    # save the plot
    pylab.savefig(plotfile, transparent=True)



def PlotEquilibriumFreqs(pis, plotfile, title='', pi_errs=None):
    """Plots equilibrium frequencies of different amino acids.
    
    *pis* : a dictionary keyed by amino acid codes with values giving
    the equilibrium frequency.

    *plotfile* : the name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    *title* is an optional argument giving the string for the title affixed
    to the top of the plot. Is an empty string ('') by default.

    *pi_errs* is an optional argument that is *None* by default. If it is
    set to another value, it should specify the error bars for the plot.
    In this case, for each key in *pis*, there should be a similar key
    in *pi_errs* with the value a list or tuple such that *pi_errs[key][0]*
    gives the lower bound and *pi_errs[key][1]* gives the upper bound.

    This function uses ``pylab`` / ``matplotlib``. It will raise an exception
    if these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('legend', fontsize=10)
    matplotlib.rc('font', size=10)
    matplotlib.rc('xtick', labelsize=10)
    fig = pylab.figure(figsize=(4.5, 2.5))
    (lmargin, rmargin, bmargin, tmargin) = (0.12, 0.02, 0.17, 0.11)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    pis_list = [tup for tup in pis.iteritems()]
    pis_list.sort()
    if pis_list[0][0] == '*':
        pis_list = pis_list[1 : ] + [pis_list[0]]
    xticks = [tup[0] for tup in pis_list]
    ys = [tup[1] for tup in pis_list]
    xindices = [i for i in range(len(pis))]
    (ymin, ymax) = (-0.03, 1.03)
    if pi_errs:
        errs = pylab.ndarray(shape=(2, len(xticks)))
        i = 0
        for x in xticks:
            if x not in pi_errs:
                raise ValueError('pi_errs does not specify a value for %s' % x)
            errs[0][i] = pis[x] - pi_errs[x][0]
            errs[1][i] = pi_errs[x][1] - pis[x]
            i += 1
        pylab.errorbar(xindices, ys, yerr=errs, fmt='o', markersize=6)
    else:
        pylab.errorbar(xindices, ys, fmt='o', markersize=6)
    pylab.gca().set_ylim([ymin, ymax])
    pylab.xticks(xindices, xticks, rotation=0)
    pylab.gca().set_xlim([-0.5, len(pis) - 0.5])
    pylab.ylabel('equilibrium preference', size=10)
    pylab.xlabel('amino acid', size=10)
    pylab.title(title, size=10)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()



def PlotEnrichmentRatios(mutations, ratios, ratio_errs, plotfile):
    """Plots enrichment ratios for different mutations on a log scale.
    
    *mutations* : a list of strings giving the names of the mutations,
    such as *['M1A', 'M1C']*.
    
    *ratios* : a list of the same lengths as *mutations* giving
    the enrichment ratios.

    *ratio_errs* : a list containing two lists each 
    of the same length as mutations giving
    the upper and lower error bar bounds at 2-tuple entries. Set this 
    to *None* if you do not want to plot any error bars.

    *plotfile* : the name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    This function uses ``pylab`` / ``matplotlib``. It will raise an exception
    if these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if ratio_errs:
        if not (len(mutations) == len(ratios) == len(ratio_errs[0]) == len(ratio_errs[1]) >= 1):
            raise ValueError("Lists not all of the same length >= 1")
    else:
        if not (len(mutations) == len(ratios) >= 1):
            raise ValueError("Lists not all of the same length >= 1")
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('legend', fontsize=10)
    matplotlib.rc('font', size=10)
    matplotlib.rc('xtick', labelsize=10)
    fig = pylab.figure(figsize=(5, 3))
    (lmargin, rmargin, bmargin, tmargin) = (0.12, 0.02, 0.2, 0.04)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    xindices = [i for i in range(len(mutations))]
    pylab.gca().set_yscale('log')
    if ratio_errs:
        ymin = min(ratio_errs[0])
        ymax = max(ratio_errs[1])
        ratio_errs[0] = [ratios[i] - ratio_errs[0][i] for i in range(len(ratios))]
        ratio_errs[1] = [-ratios[i] + ratio_errs[1][i] for i in range(len(ratios))]
        pylab.errorbar(xindices, ratios, yerr=ratio_errs, fmt='o', markersize=6)
    else:
        ymin = min(ratios)
        ymax = max(ratios)
        pylab.errorbar(xindices, ratios, fmt='o', markersize=6)
    pylab.gca().set_ylim([ymin / 2., ymax * 2.])
    pylab.xticks(xindices, mutations, rotation=90)
    pylab.gca().set_xlim([-0.5, len(mutations) - 0.5])
    pylab.ylabel('enrichment ratio', size=10)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotAAFracs(counts, names, icodon, plotfile, aas, title):
    """Plots frequency of different amino acid identities at a position.

    For a given position in a protein, plots the fraction of all
    definitively called codons (called by both reads) for a set
    of samples.

    `counts` : a list of the codon_counts dictionaries of the type
    returned by `mapmuts.io.ReadCodonCounts`.

    `names` : a list of strings giving the names of the samples
    corresponding to each entry in counts. These names are shown
    in the legend, augmented with numbers giving the total
    number of reads at this codon for that sample.

    `icodon` : the number of the codon (or residue) for which we
    are plotting the amino acid frequencies. Based on 1, 2, ...
    numbering. Each gene in `counts`  must have a codon with
    this number.

    `plotfile` : the name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    `aas` : a list of the amino acids for which we are plotting the
    the fraction. Must give one-letter amino acid codes, or 'STOP'
    for stop codons. Typically you might want this to be something
    like all non-wildtype amino acids at the position. In the legend,
    a '*' character is used to represent stop codons.

    `title` : string placed above the plot as a title.

    This function uses ``pylab`` / ``matplotlib``. It will raise an exception
    if these modules cannot be imported (if `PylabAvailable() == False`).
    """
    if len(counts) < 1:
        raise ValueError("No samples specified.")
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    if not (isinstance(counts, list) and counts):
        raise ValueError("counts does not specify a non-empty list")
    names = [name.replace('_', '\_') for name in names]
    if len(counts) != len(names):
        raise ValueError("counts and names differ in length")
    if len(aas) < 1:
        raise ValueError('aas must specify at least one amino acid.')
    if not (isinstance(icodon, int) and icodon >= 1):
        raise ValueError('icodon not an integer >= 1')
    d = dict([(sample, []) for sample in names])
    if len(d) != len(names):
        raise ValueError("Duplicate sample name in names")
    labels = []
    for (codon_counts, sample) in zip(counts, names):
        if not isinstance(codon_counts, dict):
            raise ValueError("non dictionary entry in counts")
        if icodon not in codon_counts:
            raise ValueError("No codon %d for %s" % (icodon, infile))
        icodon_counts = codon_counts[icodon]
        mapmuts.sequtils.ClassifyCodonCounts(codon_counts)
        denom = float(icodon_counts['PAIRED_COUNTS'])
        labels.append("%s ($%s$)" % (sample, \
                Base10Formatter(denom, 3, 1, 0)))
        for aa in aas:
            try:
                naa = icodon_counts['N_%s' % aa]
            except KeyError:
                raise ValueError("Failed to find a key for aa %s" % aa)
            if denom:
                d[sample].append(naa / denom)
            else:
                d[sample].append(0)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('legend', fontsize=10)
    matplotlib.rc('font', size=10)
    fig = pylab.figure(figsize=(5.6, 3.4))
    (lmargin, rmargin, bmargin, tmargin) = (0.09, 0.01, 0.35, 0.09)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    symbols = ['bo', 'gv', 'r^', 'c<', 'y>', 'ms', 'kp', 'bx', 'go', 'rv', 'c^', 'y<', 'm>', 'ks', 'bp']
    if len(names) > len(symbols):
        raise ValueError("Currently not enough symbols defined for %d names"\
                % len(names))
    xindices = [i for i in range(len(aas))]
    lines = []
    i = 0
    ymax = ymin = None
    for sample in names:
        assert len(xindices) == len(d[sample]), "%s\n%s" %\
                (str(xindices), str(d[sample]))
        line = pylab.plot(xindices, d[sample], symbols[i])
        if ymax == None:
            ymax = max(d[sample])
            ymin = min(d[sample])
        else:
            ymax = max(ymax, max(d[sample]))
            ymin = min(ymin, min(d[sample]))
        lines.append(line[0])
        i += 1
    pylab.gca().set_ylim([0, 1.08 * ymax])
    xlabels = []
    for aa in aas:
        if aa == 'STOP':
            xlabels.append('*')
        else:
            xlabels.append(aa)
    pylab.xticks(xindices, xlabels, rotation=0)
    pylab.gca().set_xlim([-0.5, len(aas) - 0.5])
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-2, 2))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    pylab.legend(lines, labels, numpoints=1, handletextpad=0.2,\
            bbox_to_anchor=(0.46, -0.1), loc='upper center', ncol=2,\
            columnspacing=1.2)
    pylab.ylabel('fraction')
    pylab.title(title.replace('_', '\_'), fontsize=11)
    if plotfile:
        pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def Base10Formatter(number, exp_cutoff, exp_decimal_digits, decimal_digits):
    """Converts a number into Latex formatting with scientific notation.

    Takes a number and converts it to a string that can be shown
    in LaTex using math mode. It is converted to scientific notation
    if the criteria specified by `exp_cutoff` is met.


    `number` : the number to be formatted, should be a float or integer.
    Currently only works for `number >= 0`.

    `exp_cutoff` : convert to scientific notation if ``abs(math.log10(number)) 
    >= exp_cutoff``.

    `exp_decimal_digits` : show this many digits after the decimal if number
    is converted to scientific notation.

    `decimal_digits` : show this many digits after the decimal if number
    is NOT converted to scientific notation.

    The returned value is a LaTex formatted string. If the number is zero, the
    returned string is simply '0'.

    EXAMPLES:

    >>> Base10Formatter(103, 3, 1, 1)
    '103.0'

    >>> Base10Formatter(103.0, 2, 1, 1)
    '1.0 \\\\times 10^{2}'

    >>> Base10Formatter(103.0, 2, 2, 1)
    '1.03 \\\\times 10^{2}'

    >>> Base10Formatter(2892.3, 3, 1, 1) 
    '2.9 \\\\times 10^{3}'

    >>> Base10Formatter(0.0, 3, 1, 1) 
    '0'

    >>> Base10Formatter(0.012, 2, 1, 1)
    '1.2 \\\\times 10^{-2}'
  
    >>> Base10Formatter(-0.1, 3, 1, 1)
    Traceback (most recent call last):
        ...
    ValueError: number must be >= 0
    """
    if number < 0:
        raise ValueError('number must be >= 0')
    if number == 0:
        return '0'
    exponent = int(math.log10(number))
    if math.log10(number) < exponent and number < 1:
        exponent -= 1
    if abs(exponent) >= exp_cutoff:
        x = number / (10.**exponent)
        formatstr = '%.' + '%d' % exp_decimal_digits + 'f \\times 10^{%d}'
        return formatstr % (x, exponent)
    else:
        formatstr = '%.' + '%d' % decimal_digits + 'f'
        return formatstr % number


def SplitLabel(label, splitlen, splitchar):
    """Splits a string with a return if it exceeds a certain length.

    `label` : a string giving the label we might split.

    `splitlen` : the maximum length of a label before we attempt to 
    split it.

    `splitchar` : the character added when splitting a label.

    If `len(label) > splitlen`, we attempt to split the label in the
    middle by adding `splitchar`. The label is split as close to the
    middle as possible while splitting at a space.

    EXAMPLES:

    No splitting as label length less than `splitlen`:

    >>> SplitLabel('WT virus 1', 10, '\\n')
    'WT virus 1'

    Splitting of this label:

    >>> SplitLabel('WT plasmid 1', 10, '\\n')
    'WT\\nplasmid 1'

    Splitting of this label:

    >>> SplitLabel('mutated WT plasmid 1', 10, '\\n')
    'mutated WT\\nplasmid 1'

    """
    if len(label) <= splitlen:
        return label
    else:
        j = 0
        imid = len(label) // 2
        index = None
        while 0 <= imid - j <= imid + j < len(label):
            if label[imid - j].isspace():
                return "%s%s%s" % (label[ : imid - j], splitchar, label[imid - j + 1 : ])
            elif label[imid + j].isspace():
                return "%s%s%s" % (label[ : imid + j], splitchar, label[imid + j + 1 : ])
            j += 1
        else:
            return label # no white space to split


def PlotCorrelation(xs, ys, plotfile, xlabel, ylabel, logx=False, logy=False,\
        corr=None, title=False, alpha=1.0, symmetrize=False, fixaxes=False):
    """Plots the correlation between two variables as a scatter plot.
    
    The data is plotted as a scatter plot.

    This function uses ``pylab`` / ``matplotlib``. It will raise an Exception if
    these modules cannot be imported (if *PylabAvailable() == False*).

    The calling variables use LaTex format for strings. So for
    example, '$10^5$' will print the LaTex equivalent of this 
    string. Similarly, certain raw text strings (such as those
    including underscores) will cause problems if you do not
    escape the LaTex format meaning. For instance, 'x_label'
    will cause a problem since underscore is not valid
    outside of math mode in LaTex, so you would need to use
    'x\_label' to escape the underscore.

    CALLING VARIABLES:

    * *xs* and *ys* are lists of numbers, with the lists
      being of the same length. Entry *xs[i]* is plotted
      on the x-axis agains entrie *ys[i]* on the y-axis.

    * *plotfile* is a string giving the name of the plot PDF file
      that we create. It should end in the extension ``.pdf``.
      If this plot already exists, it is overwritten.

    * *xlabel* is a string giving the label placed on the x-axis.

    * *ylabel* is a string giving the label placed on the y-axis.

    * *logx* specifies that we log transform the data in *xs*.
      This is *False* by default; set to *True* if you want
      to log transform (base 10 logarithms) the data.

    * *logy* is like *logx*, but for the data in *ys*.

    * *corr* specifies if we calculate and include a correlation
      coefficient on the plot. If it is *None*, then no
      correlation is computed. Otherwise, the coefficient
      is calculated using ``scipy`` (so this requires
      *mapmuts.bayesian.ScipyAvailable() == True*). In this case, *corr* should
      be set to the string *Pearson* (to calculate the 
      Pearson linear correlation coefficient) or to the string
      *Spearman* (to calculate Spearman's rho rank-order correlation).
      If *logx* and / or *logy* are true, then the data are log
      transformed BEFORE calculating the correlation, so the coefficient
      is for correlation between the log-transformed data.
      In both cases, the correlations are reported along with the
      two-tailed P-values. They are written on the plot.
      If you set this to a value other than *None* and ``scipy`` is
      not available, then an exception is raised.

    * *title* is a string giving the title placed above the plot. 
      It can be *False* if no title is to be used. Otherwise, it should
      be the title string (using LaTex formatting, spaces are allowed).
      Is *False* by default.

    * *alpha* is the transparency of the plotted points. By default
      it is one, which means that there is no transparency. If you make
      the value closer to zero, then the points become partially transparent.
      The rationale is that if you have many overlapping points, making
      them partially transparent helps you better see the density.
      At any position on the plot, the intensity will be saturated
      when there are 1.0 / alpha points plotted. So a reasonable value
      of *alpha* might be something like 0.1.

    * *symmetrize* is an optional argument that is *False* by default.
      If *True*, we make the X and Y limits on the plot the same.

    * *fixaxes* is an optional argument that is *False* by default.
      If *True*, we fix both the X and Y axes to go from 0 to 1, 
      with ticks at 0, 0.5, and 1. If you set this option to *True*,
      then you must set *logx* and *logy* to *False*.

    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if corr and not _scipyavailable:
        raise ImportError("Cannot use corr option since scipy is not available")
    if not os.path.splitext(plotfile)[1].upper() == '.PDF':
        raise ValueError("plotfile does not end in PDF extension: %s " % plotfile)
    if not (len(xs) == len(ys) >= 2):
        raise ValueError("xs and ys do not specify lists of the same length with >= 2 entries")
    if fixaxes and (logy or logx):
        raise ValueError("Cannot use fixaxes with logx or logy")
    (bigmargin, smallmargin) = (0.24, 0.05)
    (lmargin, rmargin, bmargin, tmargin) = (bigmargin, smallmargin, bigmargin, smallmargin)
    titlemargin = 0.09
    plotmargin = 0.02 # add this much above and below the last data point
    logplotmargin = 2 # scale limits by this much if log scale
    xsize = 2.0
    if title:
        tmargin += titlemargin
    ysize = xsize * (1.0 - lmargin - rmargin) / (1.0 - tmargin - bmargin)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)
    matplotlib.rc('legend', fontsize=10)
    figure = pylab.figure(figsize=(xsize, ysize), facecolor='white')
    ax = pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
    pylab.plot(xs, ys, 'b.', markersize=4, alpha=alpha)
    (xmin, xmax, ymin, ymax) = (min(xs), max(xs), min(ys), max(ys))
    if fixaxes:
        xmin = ymin = 0.0
        xmax = ymax = 1.0
    elif symmetrize:
        xmin = ymin = min(xmin, ymin)
        xmax = ymax = max(xmax, xmin)
    if logy:
        pylab.gca().set_yscale('log')
        ax.set_ylim([ymin / logplotmargin, ymax * logplotmargin])
        ys = [math.log(y) for y in ys]
    else:
        ymargin = plotmargin * (ymax - ymin)
        ax.set_ylim([ymin - ymargin, ymax + ymargin])
    if logx:
        pylab.gca().set_xscale('log')
        ax.set_xlim([xmin / logplotmargin, xmax * logplotmargin])
        xs = [math.log(x) for x in xs]
    else:
        xmargin = plotmargin * (xmax - xmin)
        ax.set_xlim([xmin - xmargin, xmax + xmargin])
    pylab.xlabel(xlabel, size=10)
    pylab.ylabel(ylabel, size=10)
    if title:
        pylab.title(title, size=10)
    if corr:
        if corr == 'Pearson':
            (r, p) = scipy.stats.pearsonr(xs, ys)
            r = '$R = %.2f$' % r
        elif corr == 'Spearman':
            (r, p) = scipy.stats.spearmanr(xs, ys)
            r = '$\\rho = %.2f$' % r
        else:
            raise ValueError("Invalid value of %s for corr" % corr)
        if p < 1e-10:
            p = '$P < 10^{-10}$'
        else:
            p = '$P = %s$' % Base10Formatter(p, 2, 1, 2)
        n = '$N = %d$' % len(xs)
# don't include N        text = '%s\n%s\n%s' % (r, p, n)
        text = '%s\n%s' % (r, p)
        pylab.text(0.05, 0.96, text, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, size=10)
    if logy:
        yticker = matplotlib.ticker.LogLocator(numticks=5)
    elif fixaxes:
        yticker = matplotlib.ticker.FixedLocator([0, 0.5, 1])
    else:
        yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    if logx:
        xticker = matplotlib.ticker.LogLocator(numticks=5)
    elif fixaxes:
        xticker = matplotlib.ticker.FixedLocator([0, 0.5, 1])
    else:
        xticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().xaxis.set_major_locator(xticker)
    pylab.gca().get_xaxis().tick_bottom()
    pylab.gca().get_yaxis().tick_left()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()



def PlotLinearDensity(datalist, plotfile, xlabel, ylabel, title=False, fixymax=False):
    """Plots linear density of variable as a function of primary sequence.

    This function is designed to plot some variable (such as the number
    of epitopes as a function of the primary sequence position). It
    creates an output PDF plot *plotfile*. 

    The data is plotted as lines. If there is more than one
    data series to be plotted, a legend is included. 

    This function uses pylab / matplotlib. It will raise an Exception if
    these modules cannot be imported (if *PylabAvailable() == False*).

    The calling variables use LaTex format for strings. So for
    example, '$10^5$' will print the LaTex equivalent of this 
    string. Similarly, certain raw text strings (such as those
    including underscores) will cause problems if you do not
    escape the LaTex format meaning. For instance, 'x_label'
    will cause a problem since underscore is not valid
    outside of math mode in LaTex, so you would need to use
    'x\_label' to escape the underscore.

    CALLING VARIABLES:

    * *datalist*  is a list specifying the data to plot. It should
      be a list of one or more 2-tuples of the form 
      *(label, data)* where *label* is a string label used in
      the legend, and data is a list of 2-tuples *(x, y)*
      specifying the points to be plotted. 

    * *plotfile* is a string giving the name of the plot PDF file
      that we create. It should end in the extension ``.pdf``.
      If this plot already exists, it is overwritten.

    * *xlabel* is a string giving the label placed on the x-axis.

    * *ylabel* is a string giving the label placed on the y-axis.

    * *title* is a string giving the title placed above the plot. 
      It can be *False* if no title is to be used. Otherwise, it should
      be the title string (using LaTex formatting, spaces are allowed).
      Is *False* by default.

    * *fixymax* means that we fix the y-maximum to the specified value.
      This may be useful if you are making multiple plots for comparisons
      between them, and want them all to have the same y-maximum.
      Note that the value specified here is taken to be the data
      maximum -- the actually maximum of the y-axis is somewhat
      higher to provide some padding space. Is *False* by default.
    
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if not os.path.splitext(plotfile)[1].upper() == '.PDF':
        raise ValueError("plotfile does not end in PDF extension: %s " % plotfile)
    if not (isinstance(datalist, list) and len(datalist) >= 1):
        raise ValueError("Invalid datalist")
    linestyles = ['b-', 'r:', 'g--', 'm-', 'c:', 'y--']
    if len(datalist) > len(linestyles):
        raise ValueError("Number of data entries exceeds linestyles available for plotting")
    (lmargin, rmargin, tmargin, bmargin) = (0.1, 0.02, 0.03, 0.17)
    if title:
        tmargin = 0.11
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=9)
    matplotlib.rc('legend', fontsize=10)
    figure = pylab.figure(figsize=(5, 2), facecolor='white')
    ax = pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
    pylab.xlabel(xlabel, size=10)
    pylab.ylabel(ylabel, size=10)
    xmin = xmax = ymin = ymax = None
    for ((label, data), style) in zip(datalist, linestyles):
        xs = [tup[0] for tup in data]
        if xmin == None:
            (xmin, xmax) = (min(xs), max(xs))
        else:
            (xmin, xmax) = (min(xmin, min(xs)), max(xmax, max(xs)))
        ys = [tup[1] for tup in data]
        if ymin == None:
            (ymin, ymax) = (min(ys), max(ys))
        else:
            (ymin, ymax) = (min(ymin, min(ys)), max(ymax, max(ys)))
        pylab.plot(xs, ys, style, label=label)
    if title:
        pylab.title(title, size=10)
    ax.set_xlim([xmin, xmax])
    if fixymax != False:
        ymax = fixymax
    if len(datalist) > 1:
        ax.set_ylim([0, 1.25 * ymax]) # higher value to allow space for legend
        pylab.legend(loc='upper center', ncol=3, borderaxespad=0)
    else:
        ax.set_ylim([0, 1.1 * ymax]) # higher value to allow space for legend
    yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def KyteDoolittleColorMapping(maptype='jet', reverse=True):
    """Maps amino-acid hydrophobicities to colors.

    Uses the Kyte-Doolittle hydrophobicity scale defined by 
    *mapmuts.sequtils.KyteDoolittle*.

    The returned variable is the 3-tuple *(cmap, mapping_d, mapper)*:

        * *cmap* is a ``pylab`` *LinearSegmentedColorMap* object.

        * *mapping_d* is a dictionary keyed by the one-letter amino-acid
          codes. The values are the colors in CSS2 format (e.g. #FF0000
          for red) for that amino acid. The value for a stop codon 
          (denoted by a * character) is black (#000000).

        * *mapper* is the actual *ScalarMappable* object.

    The optional calling argument *maptype* should specify a valid
    ``pylab`` color map. Is *jet* by default.

    The optional calling argument *reverse* specifies that we set up the color
    map so that the most hydrophobic residue comes first (in the Kyte-Doolittle
    scale the most hydrophobic comes last as it has the largest value). This option
    is *True* by default as it seems more intuitive to have charged residues red
    and hydrophobic ones blue.
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    hydrophobicities = [mapmuts.sequtils.KyteDoolittle(aa) for aa in mapmuts.sequtils.AminoAcids()]
    if reverse:
        hydrophobicities = [-1 * x for x in hydrophobicities]
    mapper = pylab.cm.ScalarMappable(cmap=maptype)
    mapper.set_clim(min(hydrophobicities), max(hydrophobicities))
    mapping_d = {'*':'#000000'}
    for (aa, h) in zip(mapmuts.sequtils.AminoAcids(), hydrophobicities):
        tup = mapper.to_rgba(h, bytes=True)
        (red, green, blue, alpha) = tup
        mapping_d[aa] = '#%02x%02x%02x' % (red, green, blue)
        assert len(mapping_d[aa]) == 7
    cmap = mapper.get_cmap()
    return (cmap, mapping_d, mapper)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
