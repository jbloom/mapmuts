"""Script to run mapmuts python programs on Illumina sequencing results.

Written by Bargavi Thyagarajan, 2014."""


import os

def RunScript(script, infilename, commands, subdir, output):
    """Runs a script in a subdirectory.
  
    *script* is the string name of a script.

    *infilename* is the name of the created input file.

    *commands* is a list of the lines to place in the input file.
  
    *subdir* is the name of the subdirectory in
    which the script is run. Created if it doesn't
    exist.

    *output* is the name of the expected output file in *subdir*.
    This script will raise an error if this file is not created.

    After going to subdirectory *subdir*, 
    runs *script* on *infile* using the command::
  
      python *script* *infile*
  
    """
    print "\nRunning script %s in subdirectory %s..." % (script, subdir)
    if not os.path.isdir(subdir):
        os.mkdir(subdir)
    cwd = os.getcwd()
    os.chdir(subdir)
    if os.path.isfile(output):
        print "Output %s already exists." % output 
    else:
        open(infilename, 'w').write('\n'.join(commands))
        os.system('%s %s' % (script, infilename))
        if not (os.path.isfile(output) or os.path.isfile(os.path.basename(output))):
            raise IOError("Script %s failed to create output %s in subdirectory %s" % (script, output, subdir))
    os.chdir(cwd)
    print "Script complete, created output %s." % output



def main():
    """Main body of script."""
    # list of setup variables
    replicates = ['replicate_1', 'replicate_2', 'replicate_3', 'replicate_1_repeat']
    unique_replicates = ['replicate_1', 'replicate_2', 'replicate_3']
    replicate_names = {'replicate_1':'\#1',
                       'replicate_2':'\#2',
                       'replicate_3':'\#3',
                       'replicate_1_repeat':'\#1 repeat',
                      }
    samples = ['DNA', 'mutDNA', 'virus', 'mutvirus']
    inferaakeys = ['DNA', 'RNA', 'mutDNA', 'mutvirus'] # key inputs for infer amino acids script
    inferaasamples = ['DNA', 'virus', 'mutDNA', 'mutvirus'] # sample inputs for infer amino acids script
    use_existing_output = True
    fullgenefile = ' %s/WSN-HA-amplicon.txt' % os.getcwd()
    sitesrange = (2, 565) # range of residues to consider
    generange = (33, 1727) # range of nucleotides in the amplicon to consider for alignment
    

    # make the alignments
    print "\nMaking alignments..."
    for replicate in replicates:
        if not os.path.isdir(replicate):
            os.mkdir(replicate)
        for sample in samples:
            subdir = "%s/%s" % (replicate, sample)
            outfileprefix = '%s_%s' % (replicate, sample)
            output = "%s_alignments.txt.gz" % outfileprefix
            infilename = 'makealignments_infile.txt'
            commands = [
                       'r1files %s/FASTQ_files/%s/%s/*R1*.gz' % (os.getcwd(), replicate, sample),
                       'r2files %s/FASTQ_files/%s/%s/*R2*.gz' % (os.getcwd(), replicate, sample),
                       'gzipped True',
                       'applyfilter True',
                       'minq 25',
                       'fullgenefile %s' % fullgenefile,
                       'generange %d %d' % generange,
                       'a1file %s/R1_trim3.fasta' % os.getcwd(),
                       'a2file %s/R2_trim3.fasta' % os.getcwd(),
                       'maxn 2',
                       'minoverlap 30',
                       'maxrm 1',
                       'maxa1m 1',
                       'maxa2m 1',
                       'maxgenem 6',
                       'upcase test',
                       'outfileprefix %s' % outfileprefix,
                       'samplename %s %s' % (replicate_names[replicate], sample),
                       'write_unaligned True',
                       ]
            RunScript('mapmuts_makealignments.py', infilename, commands, subdir, output)
    print "\nCompleted mapmuts_makealignments script."

    # Creating the alignmentsummaryplots
    commands = []
    for sample in samples:
        output = "alignmentsummaryplot.pdf" 
        infilename = 'alignmentsummary_infile.txt'
        for replicate in replicates:
            commands.append('%s/%s/%s_%s_alignmentstatistics.txt %s %s' % (replicate, sample, replicate, sample, replicate_names[replicate], sample))
    alignmentsummaryplot = 'alignmentsummaryplot.pdf'
    commands.append('plotfile %s' % alignmentsummaryplot)
    if os.path.isfile(alignmentsummaryplot) and use_existing_output:
        print "\nThe alignment summary plot %s already exists." % alignmentsummaryplot
    else:
        print "\nGenerating the alignment summary plot %s." % alignmentsummaryplot
        RunScript('mapmuts_alignmentsummaryplot.py', infilename, commands, './', output)

    # Parse the counts
    print "\nParsing counts..."
    for replicate in replicates:
        if not os.path.isdir(replicate):
            os.mkdir(replicate)
        for sample in samples:
            subdir = "%s/%s" % (replicate, sample)
            outfileprefix = '%s_%s' % (replicate, sample)
            output = "%s/%s_aacounts.txt" % (subdir, outfileprefix) 
            infilename = 'parsecounts_infile.txt'
            commands = [
                       'fullgenefile %s' % fullgenefile,
                       'generange %d %d' % generange,
                       'upcase test',
                       'alignmentfile %s/%s/%s/%s_%s_alignments.txt.gz' % (os.getcwd(), replicate, sample, replicate, sample),
                       'r1exclude None',
                       'r2exclude None',
                       'outfileprefix %s' % outfileprefix,
                       'samplename %s %s' % (replicate_names[replicate], sample),
                       ]
            if os.path.isfile(output) and use_existing_output:
                print "\nUsing existing parse counts output %s." % output
            else:
                print "\nParsing counts to generate %s." % output
                RunScript('mapmuts_parsecounts.py', infilename, commands, subdir, output)
    print "\nCompleted mapmuts_parsecounts script."

    # Creating the parsesummaryplots
    commands = []
    for sample in samples:
        output = "parsesummary_codonnmuts.pdf" 
        infilename = 'parsesummaryplots_infile.txt'
        for replicate in replicates:
            commands.append('%s/%s/%s_%s %s %s' % (replicate, sample, replicate, sample, replicate_names[replicate], sample))
    commands.append('plotfileprefix parsesummary')
    commands.append('writefracs False')
    commands.append('figwidth 8.5')
    commands.append('textwritefracs parsesummary')
    commands.append('pairedcodonplot True')
    if os.path.isfile(output) and use_existing_output:
        print "\nThe parsesummary plot %s already exists." % output 
    else:
        print "\nGenerating the parsesummary plot %s." % output
        RunScript('mapmuts_parsesummaryplots.py', infilename, commands, './', output)

    # Creating the countparsedmuts plots
    commands = ['plotfileprefix countparsedmuts',
                'maxn 50',
                'legendloc right',
                'sites %d %d' % sitesrange,
                'writecounts False']
    for sample in samples:
        output = "countparsedmuts_codonmutcounts.pdf" 
        infilename = 'countparsedmuts_infile.txt'
        files = []
        for replicate in unique_replicates:
            files.append('%s/%s/%s_%s_codoncounts.txt' % (replicate, sample, replicate, sample))
        commands.append('%s %s' % (sample, ' '.join(files)))
    if os.path.isfile(output) and use_existing_output:
        print "\nThe countparsedmuts plot %s already exists." % output
    else:
        print "\nGenerating the countparsedmuts plot %s." % output
        RunScript('mapmuts_countparsedmuts.py', infilename, commands, './', output)
    for replicate in replicates:
        commands = ['plotfileprefix countparsedmuts',
                'maxn 50',
                'legendloc right',
                'writecounts False',
                'sites %d %d' % sitesrange]
        for sample in samples:
            output = "%s/%s/countparsedmuts_codonmutcounts.pdf" % (os.getcwd(), replicate)
            infilename = '%s/%s/countparsedmuts_infile.txt' % (os.getcwd(), replicate)
            commands.append('%s %s/%s/%s/%s_%s_codoncounts.txt' % (sample, os.getcwd(), replicate, sample, replicate, sample))
        if os.path.isfile(output) and use_existing_output:
            print "\nThe countparsedmuts plot %s already exists." % output
        else:
            print "\nGenerating the countparsedmuts plot %s." % output
            RunScript('mapmuts_countparsedmuts.py', infilename, commands, replicate, output)

    # Infer amino acid preferences
    for replicate in replicates:
        commands = []
        print "\nInferring preferences for %s..." % replicate
        if not os.path.isdir(replicate):
            os.mkdir(replicate)
        if not os.path.isdir("%s/MCMC_traces" % replicate):
            os.mkdir('%s/MCMC_traces' % replicate)
        if not os.path.isdir("%s/preference_plots" % replicate):
            os.mkdir('%s/preference_plots' % replicate)
        subdir = "%s" % replicate        
        outfileprefix = '%s' % replicate
        output = "%s_equilibriumpreferences.txt" % outfileprefix
        infilename = '%s_inferpreferences_infile.txt' % replicate
        for (inferaakey, inferaasample) in zip(inferaakeys, inferaasamples): 
            commands.append('%s_files %s/%s/%s/%s_%s_codoncounts.txt' % (inferaakey, os.getcwd(), replicate, inferaasample, replicate, inferaasample))
        commands.append('mu_concentration 1.0')
        commands.append('MCMC_traces MCMC_traces')
        commands.append('preference_plots preference_plots')
        commands.append('stepincrease 4')
        commands.append('seed 1')
        commands.append('convergence 0.01')
        commands.append('ncpus 12')
        commands.append('nruns 3')
        commands.append('pi_concentration 1.0')
        commands.append('removeoutlier False')
        commands.append('rho_concentration 1.0')
        commands.append('epsilon_concentration 1.0')
        commands.append('minvalue 1e-7')
        commands.append('thin 200')
        commands.append('nsteps 200000')
        commands.append('sites %d %d' % sitesrange)
        commands.append('outfileprefix %s' % outfileprefix)
        if os.path.isfile(output) and use_existing_output:
            print "\nThe %s for  %s already exists." % (output, replicate)
        else:
            print "\nGenerating the  %s for %s." % (output, replicate)
            RunScript('mapmuts_inferpreferences.py', infilename, commands, subdir, output)

    # Preferences correlations
    correlationdir = './correlations/'
    if not os.path.isdir(correlationdir):
        os.mkdir(correlationdir)
    for i1 in range(len(replicates)):
        replicate1 = replicates[i1]
        for replicate2 in replicates[i1 + 1 : ]:
            outfile = '%s_vs_%s.pdf' % (replicate1, replicate2)
            if os.path.isfile(outfile) and use_existing_output:
                print "Correlation file %s already exists." % outfile
            else:
                infilename = 'preferencescorrelate_%s_vs_%s_infile.txt' % (replicate1, replicate2)
                commands = [
                            'plotdir ./',
                            'samplenames %s %s' % (replicate1, replicate2),
                            'preferencesfiles %s/%s/%s_equilibriumpreferences.txt %s/%s/%s_equilibriumpreferences.txt' % (os.getcwd(), replicate1, replicate1, os.getcwd(), replicate2, replicate2),
                            'alpha 0.1',
                            ]
                RunScript('mapmuts_preferencescorrelate.py', infilename, commands, correlationdir, outfile)

    # Preference means
    output = "average_equilibriumpreferences.txt" 
    infilename = 'preferencemeans_infile.txt'
    commands = []
    files = []
    for replicate in unique_replicates:
        files.append('%s/%s_equilibriumpreferences.txt' % (replicate, replicate))
    commands.append('preferencefiles %s' % (' '.join(files)))
    commands.append('includestop False') 
    commands.append('outfile average_equilibriumpreferences.txt')
    if os.path.isfile(output) and use_existing_output:
        print "\nThe file %s already exists." % output
    else:
        print "\nGenerating the %s." % output
        RunScript('mapmuts_preferencemeans.py', infilename, commands, './', output)

    # Plotting site profile preferences - Logo plot
    for (numberingname, numberingscheme) in [('sequentialnumbering', 'None'), ('H3numbering', './PDB_structure/sequential_to_H3.txt')]:
        infilename = 'siteprofileplots_infile_%s.txt' % numberingname
        commands = [
               'sitepreferences average_equilibriumpreferences.txt',
               'siterange %d %d' % sitesrange,
               'dsspfile ./PDB_structure/1RVX_trimer_renumbered.dssp',
               'dsspchain A',
               'add_rsa True',
               'add_ss False',
               'add_custom antigenic Caton_H1_HA_antigenic_sites.txt receptor-binding receptor_binding_residues.txt',
               'nperline 63',
               'includestop False',
               ]
        commands.append('outfileprefix %s_' % numberingname)
        commands.append('sitenumbermapping %s' % numberingscheme)
        output = '%s_site_preferences_logoplot.pdf' % numberingname
        if os.path.isfile(output) and use_existing_output:
            print "\nThe file %s already exists." % output
        else:
            print "\nGenerating the %s." % output
            RunScript('mapmuts_siteprofileplots.py', infilename, commands, './', output)

    # comparing entropy among subsets after correcting for rsa
    for (subsetname, subsetfile) in [('nearby_antigenic', 'nearby_antigenic_sites.txt'), ('antigenic', 'Caton_H1_HA_antigenic_sites.txt'), ('receptor_binding', 'receptor_binding_residues.txt'), ('allRBS', 'allRBS_residues.txt')]:
        plotfile = '%s_entropy_rsa_correlation.pdf' % subsetname
        commands = [
               'dsspfile ./PDB_structure/1RVX_trimer_renumbered.dssp',
               'dsspchain A',
               'aapreferences average_equilibriumpreferences.txt',
               'siterange 2 365', # HA1
               'plotfile %s' % plotfile,
               'linearmodelfile %s_linearmodelresults.txt' % subsetname,
               'selectedsites %s' % subsetfile,
               ]
        if os.path.isfile(plotfile) and use_existing_output:
            print "\nThe file %s already exists." % plotfile
        else:
            print "\nGenerating %s." % plotfile
            RunScript('mapmuts_entropycomparison.py', 'entropy_comparison_infile_%s.txt' % subsetname, commands, './', plotfile)

    # Compare entropy among NP sites in CTL epitopes
    plotfile = 'NP_CTL_entropy_rsa_correlation.pdf'
    if not (use_existing_output and os.path.isfile(plotfile)):
        print "\nGenerating %s" % plotfile
        open('NP_stuff/NP_sites_with_CTL_epitopes.txt', 'w').write('# NP sites with at least one CTL epitope\n%s' % '\n'.join([line.split(',')[0] for line in open('NP_stuff/NP_CTL_epitope_sitecounts.csv').readlines()[1 : ] if int(line.split(',')[1]) > 1]))
        commands = [
               'dsspfile ./NP_stuff/2IQH_monomerC.dssp',
               'dsspchain None',
               'aapreferences NP_stuff/NP_amino_acid_preferences.txt',
               'siterange all',
               'plotfile %s' % plotfile,
               'linearmodelfile NP_CTL_linearmodelresults.txt',
               'selectedsites NP_stuff/NP_sites_with_CTL_epitopes.txt',
               ]
        RunScript('mapmuts_entropycomparison.py', 'entropy_comparison_infile_NP_CTL.txt', commands, './', plotfile)
    else:
        print "\nThe file %s already exists." % plotfile

    # Correlate with values of Wu et al
    plotfile = 'correlation_with_Wu_et_al.pdf'
    if os.path.isfile(plotfile) and use_existing_output:
        print "\nThe file %s already exists." % plotfile
    else:
        print "\nGenerating %s" % plotfile
        os.system('python correlate_with_Wu_data.py')

    # convert PDFs to JPGs for files to be included in README
    print "\nConverting PDFs to JPGs..."
    to_convert = [
                    'alignmentsummaryplot.pdf',
                    'replicate_3/DNA/replicate_3_DNA_codondepth.pdf',
                    'parsesummary_codon_types_and_nmuts.pdf',
                    'countparsedmuts_multi-nt-codonmutcounts.pdf',
                    'replicate_1/countparsedmuts_multi-nt-codonmutcounts.pdf',
                    'replicate_2/countparsedmuts_multi-nt-codonmutcounts.pdf',
                    'replicate_3/countparsedmuts_multi-nt-codonmutcounts.pdf',
                    'replicate_1_repeat/countparsedmuts_multi-nt-codonmutcounts.pdf',
                    'correlations/replicate_1_vs_replicate_2.pdf',
                    'correlations/replicate_1_vs_replicate_3.pdf',
                    'correlations/replicate_2_vs_replicate_3.pdf',
                    'correlations/replicate_1_vs_replicate_1_repeat.pdf',
                    'sequentialnumbering_site_preferences_logoplot.pdf',
                    'H3numbering_site_preferences_logoplot.pdf',
                    'receptor_binding_entropy_rsa_correlation.pdf',
                    'antigenic_entropy_rsa_correlation.pdf',
                    'nearby_antigenic_entropy_rsa_correlation.pdf',
                    'allRBS_entropy_rsa_correlation.pdf',
                    'NP_CTL_entropy_rsa_correlation.pdf',
                    'correlation_with_Wu_et_al.pdf',
                 ]
    for pdfname in to_convert:
        if not os.path.isfile(pdfname):
            raise IOError("Cannot find PDF %s")
        else:
            jpgname = "%s.jpg" % os.path.splitext(pdfname)[0]
            if not os.path.isfile(jpgname):
                os.system('convert -density 200 %s %s' % (pdfname, jpgname))
            if not os.path.isfile(jpgname):
                raise ValueError("Cannot find %s" % jpgname)
            else:
                print "%s now exists." % jpgname



main() # run the script





