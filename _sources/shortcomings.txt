========================================
Shortcomings
========================================

Here are some known potential shortcomings of `mapmuts`_. These shortcomings are not currently believed to pose significant problems in practice for the applications in which this program has actually been used -- but it is still good to be aware of them.

* In terms of computational speed, the reading of *.gz* files is relatively slow. This is because the `Python gzip`_ module currently has poor performance, while a bug in the `Python subprocess`_ module causes problems piping directly. See http://codebright.wordpress.com/2011/03/25/139/ and http://www.macaronikazoo.com/?p=607. It is possible that you may be able to improve this performance with future versions of `Python gzip`_ or `Python subprocess`_.

* The alignment algorithms just searches for substrings rather than global best alignments. This is not a problem for genes with little redundancy, but could pose a problem on longer templates or highly redundant ones. This will definitely become highly inefficient for much longer reads and/or templates, in which case some sort of hashing approach might be preferable.

* There is currently no mechanism for handling gaps in the alignments. In general, this just causes gapped reads to be discarded. However, it could cause incorrect calling of mutations near the end of reads when there is actually a gap. Consider two reads that actually align like this::

    ATGGGACCATTT
    ATGGGACC-TTT

  It is possible that they will actually be aligned like this::

    ATGGGACCATTT
    ATGGGACCTTT

  which implies a mismatch. This does not seem to be a problem in practice, but you never know...

* As a related issue, all alignments are performed at the nucleotide level. In principle, it might be conceptually preferable to perform the alignments at the codon level for codon mutant libraries.

* The inference of the preferences uses MCMC implemented in `pymc`_. This is definitely the slowest part of the package. It is possible that this could be accelerated. For the inference of the equilbrium preferences, the script currently uses the default `pymc`_ step method for Dirichlet variables rather than the delta-exchange approach used by ``BEAST`` and some other MCMC implementations -- that may adversely affect performance? Overall, it appears that the current approach can lead to good convergence, but further optimizations certainly might be possible.

.. include:: weblinks.txt
