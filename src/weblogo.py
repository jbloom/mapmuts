"""Module for making sequence logos with ``weblogo``.

This module uses ``weblogo`` (http://weblogo.threeplusone.com/). It has been
tested with ``weblogo`` version 3.4. It is not guaranteed to work with
other version of ``weblogo``.

Before running any function in this module, you can run *WebLogoAvailable()*
to test if the ``weblogo`` executable can be found. If this function returns
*False*, then you need to get the weblogo executable in the current search
path before running other functions in the module.

This package also imports the ``weblogo`` libraries ``weblogolib`` and ``corebio.matrix``
into ``Python``.

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

*DifferentialPreferencesLogo* : plots a sequence logo of differential amino-acid preferences.


Documentation for functions
-----------------------------
Documentation for individual functions is provided in their definitions below.

"""


import os
import shutil
import subprocess
import tempfile
import string
try:
    import pyPdf
    _pypdf_available = True
except ImportError:
    _pypdf_available = False
# the following are part of the weblogo 3.4 library
import weblogolib # weblogo library
import weblogolib.colorscheme # weblogo library
import corebio.matrix # weblogo library
import corebio.utils # weblogo library

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


def DifferentialPreferencesLogo(sites, dpi_d, plotfile, nperline, overlay, sitenumbermapping=None, numberevery=10, ydatamax=1.0):
    """Creates a logo plot of differential amino-acid preferences.

    *ydatamax* is the maximum that the logo stacks extend in the positive
    and negative directions. Is 1.0 by default.
    """
    stopchar = 'X' # character for stop codon in logo plot
    firstblankchar = 'B' # character for first blank space
    lastblankchar = 'b' # character for last blank space
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
    if 'dPI_*' in dpi_d[sites[0]]:
        includestop = True
    else:
        includestop = False
    aas = mapmuts.sequtils.AminoAcids(includestop=includestop)
    if includestop:
        aas_for_string = aas[ : -1] + [stopchar]
    else:
        aas_for_string = aas
    aas_for_string = [aa for aa in aas_for_string] + [firstblankchar, lastblankchar]
    ydatamax *= 2.0 # maximum possible range of data, multiply by two for range
    try:
        # write data into transfacfile (a temporary file)
        transfacfile = tempfile.mkstemp()[1]
        f = open(transfacfile, 'w')
        f.write('ID ID\nBF BF\nP0 %s\n' % ' '.join(aas_for_string))
        ordered_alphabets = {} # keyed by site (consecutive 0-index) with values ordered lists of aas from bottom to top
        isite = 0
        for site in sites:
            positivesum = sum([dpi_d[site]['dPI_%s' % aa] for aa in aas if dpi_d[site]['dPI_%s' % aa] > 0])
            negativesum = sum([dpi_d[site]['dPI_%s' % aa] for aa in aas if dpi_d[site]['dPI_%s' % aa] < 0])
            if abs(positivesum + negativesum) > 1.0e-5:
                raise ValueError("Differential preference sums not close to zero for site %d" % site)
            f.write('%d' % site)
            dpi_aa = []
            for aa in aas:
                y = dpi_d[site]['dPI_%s' % aa]
                dpi_aa.append((y, aa))
                f.write(' %g' % (abs(y) / float(ydatamax)))
            dpi_aa.sort()
            ordered_alphabets[isite] = [firstblankchar] + [tup[1] for tup in dpi_aa] + [lastblankchar]
            isite += 1
            f.write(' %g %g' % (0.5 * (ydatamax + 2.0 * negativesum) / ydatamax, 0.5 * (ydatamax + 2.0 * negativesum) / ydatamax))
            f.write('\n')
        f.close()
        # create web logo
        aastring = ''.join(aas_for_string)
        logoprior = weblogolib.parse_prior('equiprobable', aastring, 0)
        motif = corebio.matrix.Motif.read_transfac(open(transfacfile), aastring)
        logodata = weblogolib.LogoData.from_counts(motif.alphabet, motif, logoprior)
        logo_options = weblogolib.LogoOptions()
        logo_options.fineprint = None
        logo_options.stacks_per_line = nperline
        logo_options.stack_aspect_ratio = aspectratio
        logo_options.stack_width = stackwidth
        logo_options.unit_name = 'probability'
        logo_options.show_yaxis = False
        logo_options.yaxis_scale = ymax 
        logo_options.first_index = sites[0]
        (cmap, colormapping, mapper) = mapmuts.plot.KyteDoolittleColorMapping()
        colormapping[firstblankchar] = colormapping[lastblankchar] = '#FFFFFF' # white
        color_scheme = weblogolib.colorscheme.ColorScheme()
        for (aa, aaforstring) in zip(aas + [firstblankchar, lastblankchar], aas_for_string):
            color_scheme.groups.append(weblogolib.colorscheme.ColorGroup(aaforstring, colormapping[aa], "'%s'" % aaforstring))
        logo_options.color_scheme = color_scheme
        # add site number mapping
        if sitenumbermapping:
            annotate = []
            isite = 0
            for site in sites:
                if isite % numberevery == 0:
                    annotate.append(sitenumbermapping[site].strip())
                else:
                    annotate.append('')
                isite += 1
            logo_options.annotate = annotate
        logoformat = weblogolib.LogoFormat(logodata, logo_options)
        # _my_pdf_formatter is modified from weblogo version 3.4 source code
        # to allow custom ordering of the symbols.
        pdf = _my_pdf_formatter(logodata, logoformat, ordered_alphabets) 
        open(plotfile, 'w').write(pdf)
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
        os.remove(overlayfile)

#########################################################################
# The following code is modified from weblogo (version 3.4), which
# comes with the following license:
#
# -------------------------------- WebLogo --------------------------------

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006-2011, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without 
#  modification, are permitted provided that the following conditions are met: 
#
#  (1) Redistributions of source code must retain the above copyright notice, 
#  this list of conditions and the following disclaimer. 
#
#  (2) Redistributions in binary form must reproduce the above copyright 
#  notice, this list of conditions and the following disclaimer in the 
#  documentation and or other materials provided with the distribution. 
#
#  (3) Neither the name of the University of California, Lawrence Berkeley 
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors 
#  may be used to endorse or promote products derived from this software 
#  without specific prior written permission. 
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#  POSSIBILITY OF SUCH DAMAGE. 

# Replicates README.txt

def _my_pdf_formatter(data, format, ordered_alphabets) :
    """ Generate a logo in PDF format.
    
    Modified from weblogo version 3.4 source code.
    """
    eps = _my_eps_formatter(data, format, ordered_alphabets).decode()
    gs = weblogolib.GhostscriptAPI()    
    return gs.convert('pdf', eps, format.logo_width, format.logo_height)


def _my_eps_formatter(logodata, format, ordered_alphabets) :
    """ Generate a logo in Encapsulated Postscript (EPS)
    
    Modified from weblogo version 3.4 source code. 

    *ordered_alphabets* is a dictionary keyed by zero-indexed
    consecutive sites, with values giving order of characters
    from bottom to top.
    """
    substitutions = {}
    from_format =[
        "creation_date",    "logo_width",           "logo_height",      
        "lines_per_logo",   "line_width",           "line_height",
        "line_margin_right","line_margin_left",     "line_margin_bottom",
        "line_margin_top",  "title_height",         "xaxis_label_height",
        "creator_text",     "logo_title",           "logo_margin",
        "stroke_width",     "tic_length",           
        "stacks_per_line",  "stack_margin",
        "yaxis_label",      "yaxis_tic_interval",   "yaxis_minor_tic_interval",
        "xaxis_label",      "xaxis_tic_interval",   "number_interval",
        "fineprint",        "shrink_fraction",      "errorbar_fraction",
        "errorbar_width_fraction",
        "errorbar_gray",    "small_fontsize",       "fontsize",
        "title_fontsize",   "number_fontsize",      "text_font",
        "logo_font",        "title_font",          
        "logo_label",       "yaxis_scale",          "end_type",
        "debug",            "show_title",           "show_xaxis",
        "show_xaxis_label", "show_yaxis",           "show_yaxis_label",
        "show_boxes",       "show_errorbars",       "show_fineprint",
        "rotate_numbers",   "show_ends",            "stack_height",
        "stack_width"
        ]
   
    for s in from_format :
        substitutions[s] = getattr(format,s)

    substitutions["shrink"] = str(format.show_boxes).lower()


    # --------- COLORS --------------
    def format_color(color):
        return  " ".join( ("[",str(color.red) , str(color.green), 
            str(color.blue), "]"))  

    substitutions["default_color"] = format_color(format.default_color)

    colors = []  
    for group in format.color_scheme.groups :
        cf = format_color(group.color)
        for s in group.symbols :
            colors.append( "  ("+s+") " + cf )
    substitutions["color_dict"] = "\n".join(colors)
        
    data = []
    
    # Unit conversion. 'None' for probability units
    conv_factor = None #JDB
    #JDB conv_factor = std_units[format.unit_name]
    
    data.append("StartLine")

    seq_from = format.logo_start- format.first_index
    seq_to = format.logo_end - format.first_index +1

    # seq_index : zero based index into sequence data
    # logo_index : User visible coordinate, first_index based
    # stack_index : zero based index of visible stacks
    for seq_index in range(seq_from, seq_to) :
        logo_index = seq_index + format.first_index 
        stack_index = seq_index - seq_from
        
        if stack_index!=0 and (stack_index % format.stacks_per_line) ==0 :
            data.append("")
            data.append("EndLine")
            data.append("StartLine")
            data.append("")
        
        data.append("(%s) StartStack" % format.annotate[seq_index] )

        if conv_factor: 
            stack_height = logodata.entropy[seq_index] * std_units[format.unit_name]
        else :
            stack_height = 1.0 # Probability

        # The following code modified by JDB to use ordered_alphabets
        s_d = dict(zip(logodata.alphabet, logodata.counts[seq_index]))
        s = [(s_d[aa], aa) for aa in ordered_alphabets[seq_index]]

        # Sort by frequency. If equal frequency then reverse alphabetic
        # (So sort reverse alphabetic first, then frequencty)
        # TODO: doublecheck this actual works
        #s = list(zip(logodata.counts[seq_index], logodata.alphabet))
        #s.sort(key= lambda x: x[1])
        #s.reverse()
        #s.sort(key= lambda x: x[0])
        #if not format.reverse_stacks: s.reverse()

        C = float(sum(logodata.counts[seq_index])) 
        if C > 0.0 :
            fraction_width = 1.0
            if format.scale_width :
                fraction_width = logodata.weight[seq_index] 
            # print(fraction_width, file=sys.stderr)
            for c in s:
                data.append(" %f %f (%s) ShowSymbol" % (fraction_width, c[0]*stack_height/C, c[1]) )

        # Draw error bar on top of logo. Replaced by DrawErrorbarFirst above.
        if logodata.entropy_interval is not None and conv_factor and C>0.0:

            low, high = logodata.entropy_interval[seq_index]
            center = logodata.entropy[seq_index]
            low *= conv_factor
            high *= conv_factor
            center *=conv_factor
            if high> format.yaxis_scale : high = format.yaxis_scale 

            down = (center - low) 
            up   = (high - center) 
            data.append(" %f %f DrawErrorbar" % (down, up) )
            
        data.append("EndStack")
        data.append("")
               
    data.append("EndLine")
    substitutions["logo_data"] = "\n".join(data)  


    # Create and output logo
    template = corebio.utils.resource_string( __name__, 'template.eps', __file__).decode()
    logo = string.Template(template).substitute(substitutions)

    return logo.encode()
#
# End of code modified from weblogo    
#########################################################################


if __name__ == '__main__':
    import doctest
    doctest.testmod()
#########################################################################
