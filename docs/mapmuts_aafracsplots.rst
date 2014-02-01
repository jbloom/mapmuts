.. _mapmuts_aafracsplots.py:

===============================================
mapmuts_aafracsplots.py
===============================================


Makes plots summarizing fractions of each different mutant amino acid at each site. 
This script is designed to be run after you have already run 
:ref:`mapmuts_parsecounts.py` on several samples. Each run of that program
generates a ``*_codoncounts.txt`` file, which is used as input for this
script.

To run this script from the prompt, first create a text infile of the
format described below. Then simply type :ref:`mapmuts_aafracsplots.py`
followed by the infile name. For example, if the name is ``infile.txt``,
type::

    mapmuts_aafracsplots.py infile.txt

This script will only work if `pylab`_ / `matplotlib`_ are available.


Input file format
-------------------
The input file is a text file that should have the following entries. Blank lines or lines beginning with # are ignored (i.e. as comment lines).

The first line should specify the key *plotfileprefix* followed by the prefix given to the plot for each site. For example, to place the plots in the directory *aafracsplots* with the prefix *WT-1_aafracs* you would use::

    plotfileprefix aafracsplots/WT-1_aafracs

This would then create files ``aafracsplots/WT-1_aafracs_M1.pdf``, etc for all sites. Note that if *plotfileprefix* specifies a directory as part of the prefix, that directory must already exist (so in the above example, ``aafracsplots/`` must already exist).

There should then be a list of lines specifying the locations of the ``*_codoncounts.txt`` files created by :ref:`mapmuts_parsecounts.py`. Each line should give the name / location of the file (no spaces allowed), followed by the name given to that sample in the plots.

Finally, there should be a line with the key *printprogress* and the value being either *True* or *False*. If it is *True*, a progress line is printed as the script runs, indicating how many plots have been made already.

Here is an example input file::

    # input file to mapmuts_aafracsplots.py
    plotfileprefix aafracsplots/WT-1_aafracs
    DNA/WT-1_DNA_codoncounts.txt WT-1 DNA
    mutDNA/WT-1_mutDNA_codoncounts.txt WT-1 mutDNA
    RNA/WT-1_RNA_codoncounts.txt WT-1 RNA
    virus-p1/WT-1_virus-p1_codoncounts.txt WT-1 virus-p1
    virus-p2/WT-1_virus-p2_codoncounts.txt WT-1 virus-p2
    mutvirus-p1/WT-1_mutvirus-p1_codoncounts.txt WT-1 mutvirus-p1
    mutvirus-p2/WT-1_mutvirus-p2_codoncounts.txt WT-1 mutvirus-p2
    printprogress True

Output files
--------------
This script will generate an output plot for each site in the protein. These sites will be named with the prefix given by *plotfileprefix* and the suffix being the wildtype identity at that site and then the site number. For example, assuming the first residue in methionine (M), the first plot for the above example input file would be ``aafracsplots/WT-1_aafracs_M1A.pdf``.

In general, there are going to be quite a few of these plots (for example, there will be 500 for a 500-residue protein). This is why you may want to specify a new directory with *plotfileprefix*.

Here is an example plot:

.. figure:: aafracsplot.jpg
   :width: 60%
   :align: center
   :alt: aafracsplot.jpg


.. include:: weblinks.txt
