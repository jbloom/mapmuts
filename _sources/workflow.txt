===========================================
Workflow
===========================================

`mapmuts`_ can be run directly using
scripts that should be executable from the command line if you have installed
`mapmuts`_ into a directory that is part of the search path. You therefore
do not have to interface directly with the `Python`_ modules, but can just run these
executable scripts. 

Typically you begin with a target gene and some FASTQ files from Illumina paired-end reads. This workflow progressively takes you through alignment of the reads, calling of identities and analysis of mutation frequencies, and finally Bayesian inference of the amino-acid preferences.

Each script takes as input a single text file which specifies various input
files and parameters. You will need to create these input files manually.

1) First, you will want to align the paired reads to the target gene 
   using :ref:`mapmuts_makealignments.py`.  This script aligns paired Illumina reads to each other and then to a template gene according to various user specified parameters. It takes as input the FASTQ files giving the Illumina reads. If you are analyzing several samples, you will want to run this script on each of them.

2) After you have made alignments for several different samples, you may want to 
   compare the alignment statistics for them using :ref:`mapmuts_alignmentsummaryplot.py`.
   This script will create a summary
   plot comparing the alignment statistics for the different samples. 

3) You will then want to parse the alignments into text files that give the
   number of nucleotide and codon mutations at each position. You can do this
   using :ref:`mapmuts_parsecounts.py`. This script will give you a summary of the identities called definitively by both pairs of an Illumina read pair. 

4) If you parse alignments for several different samples, you may want to compare
   the mutation frequencies among these samples using :ref:`mapmuts_parsesummaryplots.py`.

5) If you want to check the completeness with which mutations are sampled, you can use :ref:`mapmuts_countparsedmuts.py`.

6) You may then want to infer the preference of each site for each of the amino acids using :ref:`mapmuts_inferpreferences.py`. This is done using MCMC to account for sampling errors and sequencing errors / mutations. This MCMC will take a while to run (up to a few days perhaps depending on your computer and the size of the data set). 

7) After you have inferred the preferences for several samples, you may want to examine the correlations among the inferred preferences using :ref:`mapmuts_preferencescorrelate.py`.

8) You may want to average inferences for several samples using :ref:`mapmuts_preferencemeans.py`.

9) You may want to summarize the preferences for the whole gene and create a sequence logo plot summarizing the results. This can be done using :ref:`mapmuts_siteprofileplots.py`.

.. include:: weblinks.txt
