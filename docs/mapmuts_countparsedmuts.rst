.. _mapmuts_countparsedmuts.py:

=======================================
mapmuts_countparsedmuts.py
=======================================

Determine the fraction of all and all syonymous mutations that are counted at least a certain number of times.

This script is useful if you want to summarize how extensively you have sampled all mutations in a library.

Specifically, for all codon mutations for all sites, determines the number of times the mutations are observed. Then determines the fraction of mutations that are observed at least *n* times for all *n* up the maximum occurrence of any mutation (or to the *maxn* value specified in the input file). The script does this for the following types of mutations:

    * All codon mutations.
    
    * All codon mutations that involve more than one nucleotide change to the codon.
    
    * All synonymous codon mutations.
    
    * All synonymous codon mutations that involve more than one nucleotide change to the codon.

The two categories for mutations involving multiple nucleotide changes to the same codon are useful for distinguishing true codon mutations from sequencing or reverse-transcription errors, as these only rarely induce multi-nucleotide codon changes.

The distinction between all and synonymous is useful when looking at libraries that have been subjected to selection, since synonymous mutations are generally under much weaker selection than the typical mutations.

This script is designed to be run after you have already run 
:ref:`mapmuts_parsecounts.py` on one or more sample. Each run of that program
generates a ``*_codoncounts.txt`` file. This script gathers its information about mutation counts from these ``*_codoncounts.txt`` files. Note that only definitely called codons (those indicated in all uppercase without any *N* nucleotides) are considered here.

This script also makes plots of the data shown in the aforementioned output files. Making these plots require that you have `matplotlib`_ installed.

To run this script from the prompt, first create a text infile of the
format described below. Then simply type :ref:`mapmuts_countparsedmuts.py`
followed by the infile name. For example, if the name is ``infile.txt``,
type::

    mapmuts_countparsedmuts.py infile.txt



Input file format
-------------------
In the input file, blank lines or lines that begin with # are ignored (i.e. as comment lines).

Other than these comments, the input file should have the following format:

* The first line should have the key *plotfileprefix*, and be followed by a string giving the prefix that is preprended to the `Output files`_. This prefix can contain a directory path (i.e. *my_directory/my_prefix*) or just contain the file name prefix (*my_prefix*) as you prefer.

* The next line should have the key *maxn* and be followed either by an integer >= 1 or *None*. For each value of *n* ranging from 0 to *maxn*, the script reports the fraction of mutations that occur at least that many times (occurrences :math:`\ge n`). If you specify *maxn*, then the script looks at values of *n* ranging from 0 up to and including *maxn*. If you assign *maxn* a value of *None*, then the script sets *maxn* to the maximum occurrence of any mutation (when looking at all mutations) or any synonymous mutation (when looking at synonymous mutations).

* The next line should have the key *legendloc* and be followed by either *bottom* or *right*. If it is *bottom*, then the legend is at the bottom of the plot. If it is *right*, then the legend is at the right of the plot.

* The next line is optional. If specified, it should have a key of *writecounts* and a value of *True* or *False*. If this line is not included or is *True*, then the total counts for each type of mutation is written above each plot. If this line is included and has a value of *False*, then the total counts for each type of mutation is not written above each plot.

* All of the subsequent lines should list the samples for which we count the mutations. These are listed as follows: 

    - First, provide the name that you want to assign to the sample. This name can **not** contain any spaces.

    - Next, list one or more ``*_codoncounts.txt`` files of the type created by :ref:`mapmuts_parsecounts.py`. These file names can **not** contain spaces.

  You have at least one line listing such samples, and you can provide more than one such line if you want to analyze multiple samples.

  If multiple ``*_codoncounts.txt`` files are listed, then looks at each possible synonymous mutation **only at sites where all files specify a sequence with the same wildtype codon.** Then computes the total number of times that the mutation is found over all of the files listed.


Example input file
---------------------
Here is an example input file::

    # input file for running script mapmuts_countparsedmuts.py for countparsedmuts
    plotfileprefix countparsedmuts
    maxn 50
    legendloc right
    writecounts False
    DNA /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-1/DNA/replicate_A_WT-1_DNA_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-2/DNA/replicate_A_WT-2_DNA_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-1/DNA/replicate_A_N334H-1_DNA_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-2/DNA/replicate_A_N334H-2_DNA_codoncounts.txt
    RNA /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-1/RNA/replicate_A_WT-1_RNA_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-2/RNA/replicate_A_WT-2_RNA_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-1/RNA/replicate_A_N334H-1_RNA_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-2/RNA/replicate_A_N334H-2_RNA_codoncounts.txt
    mutDNA /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-1/mutDNA/replicate_A_WT-1_mutDNA_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-2/mutDNA/replicate_A_WT-2_mutDNA_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-1/mutDNA/replicate_A_N334H-1_mutDNA_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-2/mutDNA/replicate_A_N334H-2_mutDNA_codoncounts.txt
    virus-p1 /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-1/virus-p1/replicate_A_WT-1_virus-p1_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-2/virus-p1/replicate_A_WT-2_virus-p1_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-1/virus-p1/replicate_A_N334H-1_virus-p1_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-2/virus-p1/replicate_A_N334H-2_virus-p1_codoncounts.txt
    virus-p2 /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-1/virus-p2/replicate_A_WT-1_virus-p2_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-2/virus-p2/replicate_A_WT-2_virus-p2_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-1/virus-p2/replicate_A_N334H-1_virus-p2_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-2/virus-p2/replicate_A_N334H-2_virus-p2_codoncounts.txt
    mutvirus-p1 /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-1/mutvirus-p1/replicate_A_WT-1_mutvirus-p1_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-2/mutvirus-p1/replicate_A_WT-2_mutvirus-p1_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-1/mutvirus-p1/replicate_A_N334H-1_mutvirus-p1_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-2/mutvirus-p1/replicate_A_N334H-2_mutvirus-p1_codoncounts.txt
    mutvirus-p2 /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-1/mutvirus-p2/replicate_A_WT-1_mutvirus-p2_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/WT-2/mutvirus-p2/replicate_A_WT-2_mutvirus-p2_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-1/mutvirus-p2/replicate_A_N334H-1_mutvirus-p2_codoncounts.txt /home/jbloom/mapmuts/examples/2013Analysis_Influenza_NP_Aichi68/replicate_A/N334H-2/mutvirus-p2/replicate_A_N334H-2_mutvirus-p2_codoncounts.txt

Output files
---------------
The following output files are created:

* A text file giving the fraction of all of the multi-nucleotide codon mutations to a protein of length *L* that are found :math:`\ge n` times for :math:`0 \le n \le maxn` (these are codon mutations that involve more than one nucleotide change to the same codon). The name of this file begins with the prefix specified by *plotfileprefix* and is followed by the suffix ``_multi-nt-allcodonmutcounts.txt``. 

  The first two lines are headers beginning with a *#* character. After that, the first column gives *n* and the remaining columns (there are as many remaining columns are there are samples listed in the input file) give the fraction of mutations that are found :math:`\ge` that many times for each sample. The sample columns are in the same order that they are listed in the input file. The columns are tab delimited. Here are a few example lines::

    # File listing the fraction of multi-nt-all mutations that are found greater than or equal to n times.
    # There are 26892 total multi-nt-all mutations
    #n  DNA RNA mutDNA  mutvirus-p1 mutvirus-p2 virus-p1    virus-p2
    0   1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
    1   0.1936  0.2220  0.9981  0.8801  0.8510  0.2259  0.3395
    2   0.0750  0.0950  0.9972  0.7986  0.7448  0.0926  0.1451
    3   0.0468  0.0627  0.9964  0.7317  0.6689  0.0590  0.0862
    4   0.0340  0.0489  0.9955  0.6772  0.6082  0.0451  0.0624
    5   0.0272  0.0396  0.9948  0.6302  0.5591  0.0357  0.0483
    6   0.0223  0.0333  0.9939  0.5891  0.5173  0.0297  0.0391
    7   0.0193  0.0294  0.9927  0.5540  0.4836  0.0255  0.0329
    8   0.0161  0.0259  0.9922  0.5244  0.4544  0.0217  0.0287
    9   0.0141  0.0231  0.9911  0.4980  0.4299  0.0195  0.0248
    10  0.0124  0.0206  0.9901  0.4734  0.4087  0.0173  0.0219



* A text file like the ``_multi-nt-allcodonmutcounts.txt``, but for all codon mutations (single and multi-nucleotide codon mutations). This file begins with the prefix specified by *plotfileprefix* and is followed by the suffix ``_allcodonmutcounts.txt``.

* A text file like ``*_allcodonmutcounts.txt`` but just for synonymous mutations. This file begins with the prefix specified by *plotfileprefix* and is followed by the suffix ``_syncodonmutcounts.txt``. 

* A text file like ``*_multi-nt-allcodonmutcounts.txt`` but just for synonymous mutations that involve more than one nucleotide change in the codon mutation. This file begins with the prefix specified by *plotfileprefix* and is followed by the suffix ``_multi-nt-syncodonmutcounts.txt``. 

* A plot beginning with *plotfileprefix* and followed by the suffix ``_multi-nt-codonmutcounts.pdf``. Here is an example of this plot:

    .. figure:: example_2013Analysis_Influenza_NP_Aichi68_replicate_A_countparsedmuts_multi-nt-codonmutcounts.jpg
       :width: 80%
       :align: center
       :alt: example_2013Analysis_Influenza_NP_Aichi68_replicate_A_countparsedmuts_multi-nt-codonmutcounts.jpg

* A plot like the ``_multi-nt-codonmutcounts.pdf`` but for all codon mutations regardless of whether or not they have multi-nucleotide changes. This plot has the suffix ``_codonmutcounts.pdf``.

.. include:: weblinks.txt
