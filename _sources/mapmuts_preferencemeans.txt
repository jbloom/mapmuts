.. _mapmuts_preferencemeans.py:

======================================
mapmuts_preferencemeans.py
======================================

This is a simple script that computes the mean of the amino acid preferences inferred by several runs of :ref:`mapmuts_inferpreferences.py`. 

Briefly, for each residue *r* in a protein, the script :ref:`mapmuts_inferpreferences.py` computes the preference of that site for each amino acid *a*, which is denoted as :math:`\pi_{r,a}`. These preferences sum to one, so :math:`1 = \sum_a \pi_{r,a}`. It is possible that you have run several replicates of :ref:`mapmuts_inferpreferences.py`, such as for several different replicates of an experiment. You then might want to compute the average preference for each amino acid taken over the replicates. This script does that, and writes the output to a file.

Specifically, let :math:`\pi_{r,a}^i` denote the preference for amino acid *a* at site *r* computed for replicate *i*, where :math:`i = 1, \ldots N` where *N* is the number of replicates. Then the mean preference of site *r* for amino acid *a* is :math:`\langle \pi_{r,a} \rangle = \frac{1}{N} \sum_i \pi_{r,a}^i`. Note that because the preferences sum to one for each individual library, it will also be the case that :math:`1 = \sum_a \langle \pi_{r,a} \rangle`.


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

* *preferencefiles* : The input data for this script are the ``*_equilibriumpreferences.txt`` files created by :ref:`mapmuts_inferpreferences.py`. You must specify at least two such files (if there were just one, you wouldn't need to run this script). So this key should be followed by a list of two or more ``*_equilibriumpreferences.txt`` files. These files will have a format as follows (and as detailed in the documentation for :ref:`mapmuts_inferpreferences.py`::

            #SITE   WT_AA   SITE_ENTROPY    PI_A    PI_C    PI_D    PI_E    PI_F    PI_G    PI_H    PI_I    PI_K    PI_L    PI_M    PI_N    PI_P    PI_Q    PI_R    PI_S    PI_T    PI_V    PI_W    PI_Y    PI_*        
            1   M   0.885694    0.000680047 0.00185968  8.48048e-06 0.00128217  0.0289946   0.000654589 0.0212144   0.0205306   0.0021597   0.00056726  0.876551    0.00130602  0.000789781 0.00144568  0.000414457 0.00071164  0.00126055  0.00167084  0.0352283   0.00176486  0.000905383
            2   A   1.1671  0.709575    0.00177747  3.79503e-05 0.00298854  0.00390524  0.00131246  0.00379599  0.000913006 0.00121497  0.000535209 0.00276335  0.00140061  0.000694694 0.0032134   0.00464972  0.000463924 0.000908368 0.254241    0.00335749  0.00137437  0.000876985        
            3   N   3.2593  0.129992    0.000540239 4.2514e-06  0.00171887  0.0212062   0.0006086   0.102743    0.0431489   0.00524637  0.0771062   0.021714    0.0641921   0.00338082  0.0149898   0.0106522   0.13886 0.0247841   0.00607319  0.0500692   0.282389    0.000581555

 Note that it is optional whether a preference is specified for stop codons (indicated by a * character) in the last column. But all files specified here must be consistent -- they either all must specify a preference for a stop codon, or all must not specify a preference for a stop codon.

* *outfile* is the name of the created output file. This output file is in the same format as the input files in *preferencefiles*, but now contains the average preference.

* *includestop* specifies whether we include stop codons (denoted by a * character) as possible amino acids in the generated *outfile*. If the input *preferencesfiles* do not specify stop codons as an amino acid, then no stop codons are ever included in *outfile* and this option is meaningless. But the input *preferencesfiles* do specify stop codon preferences, the *includestop* is meaningful. In this case:

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

The file specified by *outfile* is created (or overwritten if it already exists). This file has the exact same format as the input files listed in *preferencefiles*, except the indicated preferences are now the :math:`\langle \pi_{r,a} \rangle` average preferences (this format matches that for the ``*_equilibriumpreferences.txt`` files described in the documentation for :ref:`mapmuts_inferpreferences.py`. The entropies (third column) are recomputed based on the average preferences as :math:`h_r = \sum_a \langle \pi_{r,a} \rangle \log_2 \langle \pi_{r,a}\rangle`. The entries for the wildtype amino acid (second column, *WT_AA*) show either the single wildtype amino acid (if identical in all *preferencefiles*) or a comma delimited list of all wildtype amino acids in the different *preferencefiles* (such as *N,N,H,H* if the *preferencefiles* have wildtype amino acids *N*, *N*, *H*, and *H* for four such files).

The *outfile* created by this script is valid input to :ref:`mapmuts_siteprofileplots.py`.


.. include:: weblinks.txt
