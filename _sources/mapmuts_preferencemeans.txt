.. _mapmuts_preferencemeans.py:

======================================
mapmuts_preferencemeans.py
======================================

This is a simple script that computes the mean of the amino acid preferences or differential preferences inferred by several runs of :ref:`mapmuts_inferpreferences.py` or :ref:`mapmuts_inferdifferentialpreferences.py`.

Preferences
~~~~~~~~~~~~~~~
If using amino-acid preferences (*preferencefiles* argument in `Input file`_), then we are working with the following data.

For each residue *r* in a protein, the script :ref:`mapmuts_inferpreferences.py` computes the preference of that site for each amino acid *a*, which is denoted as :math:`\pi_{r,a}`. These preferences sum to one, so :math:`1 = \sum_a \pi_{r,a}`. It is possible that you have run several replicates of :ref:`mapmuts_inferpreferences.py`, such as for several different replicates of an experiment. You then might want to compute the average preference for each amino acid taken over the replicates. This script does that, and writes the output to a file.

Specifically, let :math:`\pi_{r,a}^i` denote the preference for amino acid *a* at site *r* computed for replicate *i*, where :math:`i = 1, \ldots N` where *N* is the number of replicates. Then the mean preference of site *r* for amino acid *a* is :math:`\langle \pi_{r,a} \rangle = \frac{1}{N} \sum_i \pi_{r,a}^i`. Note that because the preferences sum to one for each individual library, it will also be the case that :math:`1 = \sum_a \langle \pi_{r,a} \rangle`.

Differential preferences
~~~~~~~~~~~~~~~~~~~~~~~~~~
If using differential amino-acid preferences (*differentialpreferencefiles* argument in `Input file`_), then we are working with the following data.

For each residue :math:`r` in a protein, the script :ref:`mapmuts_inferdifferentialpreferences.py` computes the differential preference :math:`\Delta\pi_{r,a}` of that site for amino acid :math:`a`. These differential preferences sum to zero, so :math:`0 = \sum_a \Delta\pi_{r,a}`. This script is designed for the case when you have performed several runs of :ref:`mapmuts_inferdifferentialpreferences.py` on different samples, and want to compute the average.

Specifically, then :math:`\Delta\pi_{r,a}^i` denote the differential preferences for replicate :math:`i`, where :math:`i = 1, \ldots N` where :math:`N` is the number of replicates. Then the mean differential preference of site :math:`r` for amino acid :math:`a` is :math:`\langle \Delta\pi_{r,a} \rangle = \frac{1}{N} \sum_i \Delta\pi_{r,a}^i`. Note that because the differential preferences sum to zero for each individual replicate, it will also be the case that the average differential preferences sum to zero.

Dependencies
--------------
This script has no dependencies outside of `mapmuts`_ and the standard `Python`_ library.


Running the script
--------------------
To run this script from the prompt, first create a text infile of the
format described below. Then simply run the script with the single text file as the argument, as in::

    mapmuts_preferencemeans.py infile.txt


Input file
--------------
The input file is a text file with a series of *key* / *value* pairs. The required keys are indicated below. The values should not include spaces.

Lines beginning with # and empty lines are ignored.

Keys for the input file:

* *preferencefiles* : You should specify **either** this argument or *differentialpreferencefiles*, but not both. If using *preferencefiles*, this argument should be followed by a list of the ``*_equilibriumpreferences.txt`` files created by :ref:`mapmuts_inferpreferences.py`. You must specify at least two such files (if there were just one, you wouldn't need to run this script). So this key should be followed by a list of two or more ``*_equilibriumpreferences.txt`` files. 

  Note that it is optional whether a preference is specified for stop codons (indicated by a * character) in the last column. But all files specified here must be consistent -- they either all must specify a preference for a stop codon, or all must not specify a preference for a stop codon.

* *differentialpreferencefiles* : You should specify **either** this argument or *preferencefiles*, but not both. If using *differentialpreferencefiles*, this argument should be followed by a list of the ``differentialpreferences_selection_*.txt`` files created by :ref:`mapmuts_inferdifferentialpreferences.py`. You must specify at least two such files (if there were just one, you wouldn't need to run this script).

  Note that it is optional whether a preference is specified for stop codons (indicated by a * character) in the last column. But all files specified here must be consistent -- they either all must specify a preference for a stop codon, or all must not specify a preference for a stop codon.
* *outfile* is the name of the created output file. This output file is in the same format as the input files in *preferencefiles*, but now contains the average preference.

* *includestop* specifies whether we include stop codons (denoted by a * character) as possible amino acids in the generated *outfile*. If the input *preferencefiles* / *differentialpreferencefiles* do not specify stop codons as an amino acid, then no stop codons are ever included in *outfile* and this option is meaningless. But the input files do specify stop codon preferences, the *includestop* is meaningful. In this case:

    * If *includestop* is *True* then *outfile* also contains these stop codons as possible amino acids. 
    
    * If *includestop* is *False*, then stop codons are not included as possible amino acids. In this case, the :math:`\pi_{r,a}` values for the 20 non-stop amino acids are normalized so that they sum to one.

Example input file
---------------------
Here is an example input file::

    # Input file for mapmuts_preferencemeans.py
    preferencefiles WT-1_equilibriumpreferences.txt WT-2_equilibriumpreferences.txt N334H-1_equilibriumpreferences.txt N334H-2_equilibriumpreferences.txt
    outfile mean_equilibriumpreferences.txt
    includestop False


Output
--------

The script writes some brief output to standard output.

The file specified by *outfile* is created (or overwritten if it already exists). This file has the exact same format as the input files listed in *preferencefiles* or *differentialpreferencfiles*, except the indicated preferences / differential preferences are now average values. The entropies or RMS differential preferences (third columns of these files) are recomputed. 

The entries for the wildtype amino acid (second column, *WT_AA*) show either the single wildtype amino acid (if identical in all files) or a comma delimited list of all wildtype amino acids in the different files (such as *N,N,H,H* if the *preferencefiles* have wildtype amino acids *N*, *N*, *H*, and *H* for four such files).

The *outfile* created by this script is valid input to :ref:`mapmuts_siteprofileplots.py`.


.. include:: weblinks.txt
