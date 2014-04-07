"""Module for making sequence logos with ``weblogo``.

This module uses ``weblogo`` (http://weblogo.threeplusone.com/). It has been
tested with ``weblogo`` version 3.3. It is not guaranteed to work with
other version of ``weblogo``.

Before running any function in this module, you can run *WebLogoAvailable()*
to test if the ``weblogo`` executable can be found. If this function returns
*False*, then you need to get the weblogo executable in the current search
path before running other functions in the module.

If you are using an overlay on the sequence logo, the package ``pyPdf``
(http://pybrary.net/pyPdf/) is also required. You can test this with
*PyPdfAvailable()*. 

Written by Jesse Bloom.


List of functions
-------------------
*WebLogoAvailable* : tests whether ``weblogo`` can be run from the current search path.

*ConvertAvailable* : tests whether ``convert`` utility from ImageMagick is available.

*PyPdfAvailable* : tests whether *pyPdf* package is available.

*EquilibriumFreqsLogo* : plots a sequence logo of equilibrium amino-acid preferences.


Documentation for functions
-----------------------------
Documentation for individual functions is provided in their definitions below.

"""


import os
import shutil
import subprocess
import tempfile
try:
    import pyPdf
    _pypdf_available = True
except ImportError:
    _pypdf_available = False
import mapmuts.sequtils


def PyPdfAvailable():
    """Returns *True* if and only if ``pyPdf`` is available for import."""
    return _pypdf_available


def WebLogoAvailable():
    """Tests if ``weblogo`` is in the current search path.
    
    Returns *False* if ``weblogo`` cannot be run using the *os* system commands
    for running external programs. Otherwise returns a non-empty string
    giving the current ``weblogo`` version in the search path. This string
    evaluates to *True* using *bool(WebLogoAvailable())*.
    """
    sout = tempfile.TemporaryFile()
    serror = tempfile.TemporaryFile()
    code = subprocess.call(['weblogo', '--version'], stdout=sout, stderr=serror)
    sout.seek(0)
    version = sout.read()
    sout.close()
    serror.close()
    if code == 0:
        if version:
            return version
        else:
            return 'version not available'
    else:
        return False


def ConvertAvailable():
    """Returns *True* if and only if ``convert`` is in the current search path.
    
    Returns *True* if ``convert`` can be run using the *os* system commands
    for running external programs, and *False* otherwise.

    ``convert`` is the command-line access to the ImageMagick suite.
    """
    sout = tempfile.TemporaryFile()
    serror = tempfile.TemporaryFile()
    code = subprocess.call(['convert', '--version'], stdout=sout, stderr=serror)
    sout.close()
    serror.close()
    if code == 0:
        return True
    else:
        return False



def EquilibriumFreqsLogo(sites, pi_d, plotfile, nperline, overlay, sitenumbermapping=None, numberevery=10):
    """Constructs a sequence logo showing amino-acid equilibrium preferences.

    The heights of each amino-acid letter are equal to the preference of
    that site for the amino acid.

    Note that stop codons may or may not be included in the logo
    depending on whether they are present in *pi_d* as detailed in the
    explanation of that calling parameter. So if you do not want to
    include stop codons as a potential letter in the logo, then first
    remove them by calling *mapmuts.bayesian.PreferencesRemoveStop*
    on *pi_d*. Otherwise if stop codons are present, they are plotted
    using a black *X* character -- we are not able to use the normal
    asterisk character for these codons as ``weblogo`` does not appear
    to allow that.

    This function uses ``weblogo``. It will raise an exception if
    *not WebLogoAvailable()*.

    CALLING VARIABLES:

    * *sites* is a list of integers giving all of the sites that
      are being included in the logo. They must be consecutive,
      so that *sites = [i for i in range(sites[0], sites[-1] + 1)]*.

    * *pi_d* is a dictionary that has a key for every integer in
      *sites*. The value of *pi_d[isite]* is itself a dictionary,
      which has keys 'PI_A', 'PI_C', 'PI_D', etc for all 20
      one-letter upper-case amino acid codes. The values
      for these keys are the equilibrium preference of that
      amino acid at that site. So *pi_d[isite]['PI_M']* is
      the equilibrium preference for methionine at site *isite*.
      *pi_d* is allowed to either contain or not contain stop codons.
      If it contains stop codons, then there should be a key
      'PI_*' giving the preference for a stop codon for each
      dictionary *pi_d[isite]*. However, we only check that there
      are actually stop codons by looking to see if there
      is a key 'PI_*' in *pi_d[sites[0]]* -- if there is not,
      then we don't look for stop codons at any other sites either.
      Note that even though stop codons are denoted by an asterisk
      in *pi_d*, they are plotted using an *X* character in the 
      sequence logo.

    * *plotfile* is a string giving the name of the PDF file containing
      the heat map plot. It must end in the extension ``.pdf``.

    * *nperline* is the number of sites per line. You will get a good
      appearance with values in the range from 40 to 80.

    * *overlay* specifies that we also overlay a bar showing the hydrophobicity
      scale for the sequence logo colors, as well as bars showing the relative
      solvent accessibility (RSA) and either the secondary structure (SS)
      or some other discrete custom property with between 1 or more classes. This 
      requires *mapmuts.plot.PylabAvailable()* to be *True*, and also
      requires *PyPdfAvailable()* to be *True*. If you do not want to
      do the overlaying, set *overlay* to *False*. If you do want to do
      the overlaying, then set *overlay* to a list. The list must be
      composed as follows::

        overlay = [rsa_d, ss_d]

      * *rsa_d* is a dictionary giving the RSA values. It is keyed by the
        integer site numbers: *rsa_d[site]* gives a float with the RSA
        for residue *site*. If *site* is not in RSA d, the value is
        shown in white.

      * *ss_d* is like *rsa_d* except that it gives the codes for some
        discrete property. These can be secondary structure ('helix', 
        'strand', and 'loop') or some other set of one ore more strings.

    * *sitenumbermapping* is an optional argument that is *None* by
      default. If it is set to some other value then it should be 
      a dictionary keyed by every integer in *sites*. The values are
      then the string number annotation assigned to that site.

    * *numberevery* is an option that is only meaningful if *sitenumbermapping*
      is not *None*. In this case, it should be an integer. We only label
      sites at this interval. It is 10 by default, meaning that only every 10th 
      site is labeled.
    """
    stopchar = 'X' # character for stop codon in logo plot
    if not WebLogoAvailable():
        raise ValueError("Cannot run weblogo")
    if overlay and not PyPdfAvailable():
        raise ValueError("Cannot use overlay as pyPdf is not available.")
    if overlay and not mapmuts.plot.PylabAvailable():
        raise ValueError("Cannot use overlay as pylab is not available.")
    if overlay:
        if not (len(overlay) == 2 and isinstance(overlay[0], dict) and isinstance(overlay[1], dict)):
            raise ValueError("overlay is not a list of two dictionaries.")
    if sites != [i for i in range(sites[0], sites[-1] + 1)]:
        raise ValueError("sites does not specify consecutive numbers")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    if os.path.isfile(plotfile):
        os.remove(plotfile) # remove existing plot
    #
    # Following are specifications of weblogo sizing taken from its documentation
    # or specified when weblogo is called
    stackwidth = 9.5 # stack width in points, not default size of 10.8, but set to this in weblogo call below
    barheight = 5.5 # height of bars in points if using overlay
    barspacing = 2.0 # spacing between bars in points if using overlay
    stackaspectratio = 4.4 # ratio of stack height:width, doesn't count part going over maximum value of 1
    if overlay:
        ymax = (stackaspectratio * stackwidth + len(overlay) * (barspacing + barheight)) / float(stackaspectratio * stackwidth)
        aspectratio = ymax * stackaspectratio # effective aspect ratio for full range
    else:
        ymax = 1.0
        aspectratio = stackaspectratio
    rmargin = 11.5 # right margin in points, fixed by weblogo
    stackheightmargin = 16 # margin between stacks in points, fixed by weblogo
    # End specifications of weblogo sizing taken from its documentation
    #
    assert sites, "No sites specified"
    if 'PI_*' in pi_d[sites[0]]:
        includestop = True
    else:
        includestop = False
    aas = mapmuts.sequtils.AminoAcids(includestop=includestop)
    if includestop:
        aas_for_string = aas[ : -1] + [stopchar]
    else:
        aas_for_string = aas
    try:
        # write data into transfacfile (a temporary file)
        transfacfile = tempfile.mkstemp()[1]
        f = open(transfacfile, 'w')
        f.write('ID ID\nBF BF\nP0 %s\n' % ' '.join(aas_for_string))
        for site in sites:
            f.write('%d %s\n' % (site, ' '.join([str(pi_d[site]['PI_%s' % aa]) for aa in aas])))
        f.close()
        # set up weblogo calling arguments
        aastring = ''.join(aas_for_string)
        args = ['weblogo', '-f', transfacfile, '-D', 'transfac', '-o', plotfile, '-F', 'pdf', '--weight', '0', '--alphabet', aastring, '-i', str(sites[0]), '-n', str(nperline), '--units', 'probability', '--composition', 'equiprobable', '-Y', 'NO', '--fineprint', '', '--errorbars', 'NO', '--stack-width', str(stackwidth), '--aspect-ratio', str(aspectratio), '--yaxis', str(ymax)]
        # add color mapping to weblogo calling arguments
        (cmap, colormapping, mapper) = mapmuts.plot.KyteDoolittleColorMapping()
        for (aa, aaforstring) in zip(aas, aas_for_string): 
            args.append('--color')
            args.append(colormapping[aa])
            args.append(aaforstring)
            args.append("'%s'"% aaforstring)
        # add site number mapping
        if sitenumbermapping:
            args.append('--annotate')
            annotatestring = []
            isite = 0
            for site in sites:
                if isite % numberevery == 0:
                    annotatestring.append(sitenumbermapping[site].strip())
                else:
                    annotatestring.append('')
                isite += 1
            args.append(','.join(annotatestring))
        # run weblogo
        sout = tempfile.TemporaryFile()
        serror = tempfile.TemporaryFile()
        code = subprocess.call(args, stdout=sout, stderr=serror)
        # check for errors running weblogo
        if code != 0 or not os.path.isfile(plotfile):
            sout.seek(0)
            serror.seek(0)
            raise IOError("Failed to run weblogo successfully.\nHere is standard out:\n%s\nHere is standard error:\n%s" % (sout.read(), serror.read()))
        serror.close()
        sout.close()
    finally:
        # remove temporary file
        if os.path.isfile(transfacfile):
            os.remove(transfacfile)
    # now build the overlay
    if overlay:
        # make the overlay plot
        overlayfile = '_overlay_tempfile.pdf'
        mergedfile = '_merged_tempfile.pdf'
        mapmuts.plot.LogoOverlay(sites, overlayfile, overlay[0], overlay[1], nperline, sitewidth=stackwidth, rmargin=rmargin, logoheight=stackwidth * stackaspectratio + stackheightmargin, barheight=barheight, barspacing=barspacing)
        # overlay onto plotfile using pyPdf
        plot = pyPdf.PdfFileReader(open(plotfile, 'rb')).getPage(0)
        overlay = pyPdf.PdfFileReader(open(overlayfile, 'rb')).getPage(0)
        xshift = overlay.artBox[2] - plot.artBox[2]
        overlay.mergeTranslatedPage(plot, xshift, 0)
        output = pyPdf.PdfFileWriter()
        output.addPage(overlay)
        outputstream = open(mergedfile, 'wb')
        output.write(outputstream)
        outputstream.close()
        os.rename(mergedfile, plotfile)
# COMMENTED OUT OLD CODE DOING THIS WITH CONVERT. THE PROBLEM WAS IMAGE RASTERIZATION
#        args = ['convert', '-density', '600', '-transparent', 'white', overlayfile, '-gravity', 'SouthEast', plotfile, '-transparent', 'white', '-composite', plotfile]
#        sout = tempfile.TemporaryFile()
#        serror = tempfile.TemporaryFile()
#        code = subprocess.call(args, stdout=sout, stderr=serror)
#        if code != 0:
#            sout.seek(0)
#            serror.seek(0)
#            raise IOError("Failed to run convert successfully.\nHere is standard out:\n%s\nHere is standard error:\n%s" % (sout.read(), serror.read()))
#        serror.close()
#        sout.close()
        os.remove(overlayfile)




if __name__ == '__main__':
    import doctest
    doctest.testmod()
