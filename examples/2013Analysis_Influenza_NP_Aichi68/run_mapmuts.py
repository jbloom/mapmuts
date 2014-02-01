"""This script runs a suite of ``mapmuts`` programs.

This is a master Python script that runs the ``mapmuts`` scripts to analyze
a deep sequencing library. 

Can run multiple scripts at the same time, using either direction submission
of the jobs on the current CPU or submission to a queue via ``sbatch``.

The locations and names of various files are hard-coded into the script.

Because this is my personal script for running jobs, the code is not
fully documented elsewhere, and you will have to look at the script source
code to understand exactly how it works. It has been tested 
on the FHCRC's computing core. Note however that the ``mapmuts`` scripts
run by this code are fully documented with the ``mapmuts`` package, and
do not requiring understanding the source code. This script is just running
the ``mapmuts`` scripts.

Written by Jesse Bloom.
"""


import os
import string
import time
import multiprocessing


def RunScript(rundir, run_name, script_name, commands, use_sbatch, sbatch_cpus, walltime=None):
    """Runs a ``mapmuts`` script.

    *rundir* is the directory in which we run the job. Created if it does
    not exist.

    *run_name* is the name of the run, which should be a string without
    spaces. The input file has this prefix followed by ``_infile.txt``.

    *script_name* is the name of the script that we run.

    *commands* contains the commands written to the input file. Itis a list 
    of 2-tuples giving the key / value pairs.
    Both keys and values should be strings.

    *use_sbatch* is a Boolean switch specifying whether we use ``sbatch``
    to run the script. If *False*, the script is just run with the command
    line instruction. If *True*, then ``sbatch`` is used, and the command file
    has the prefix *run_name* followed by the suffix ``.sbatch``.

    *sbatch_cpus* is an option that is only meaningful if *use_sbatch* is 
    *True*. It gives the integer number of CPUs that are claimed via
    ``sbatch`` using the option ``sbatch -c``. 

    *waltime* is an option that is only meaningful if *use_sbatch* is
    *True*. If so, it should be an integer giving the number of hours 
    to allocate for the job. If *walltime* has its default value of 
    *None*, no wall time for the job is specified.

    It is assumed that the script can be run at the command line using::

        script_name infile

    Returns *runfailed*: *True* if run failed, and *False* otherwise.
    """
    print "Running %s for %s in directory %s..." % (script_name, run_name, rundir)
    currdir = os.getcwd()
    if not os.path.isdir(rundir):
        os.mkdir(rundir)
    os.chdir(rundir)
    if (not run_name) or not all([x not in string.whitespace for x in run_name]):
        raise ValueError("Invalid run_name of %s" % run_name)
    infile = '%s_infile.txt' % run_name
    open(infile, 'w').write('# input file for running script %s for %s\n%s' % (script_name, run_name, '\n'.join(['%s %s' % (key, value) for (key, value) in commands])))
    if use_sbatch:
        sbatchfile = '%s.sbatch' % run_name # sbatch command file
        jobidfile = 'sbatch_%s_jobid' % run_name # holds sbatch job id
        jobstatusfile = 'sbatch_%s_jobstatus' % run_name # holds sbatch job status
        joberrorsfile = 'sbatch_%s_errors' % run_name # holds sbatch job errors
        sbatch_f = open(sbatchfile, 'w')
        sbatch_f.write('#!/bin/sh\n#SBATCH\n')
        if walltime:
            sbatch_f.write('#PBS -l walltime=%d:00:00\n' % walltime)
        sbatch_f.write('%s %s' % (script_name, infile))
        sbatch_f.close()
        os.system('sbatch -c %d -e %s %s > %s' % (sbatch_cpus, joberrorsfile, sbatchfile, jobidfile))
        time.sleep(1) # short 1 second delay
        jobid = int(open(jobidfile).read().split()[-1])
        while True:
            time.sleep(1) # delay 1 second
            os.system('squeue -j %d > %s' % (jobid, jobstatusfile))
            lines = open(jobstatusfile).readlines()
            if len(lines) < 2:
                break # no longer in slurm queue
        errors = open(joberrorsfile).read().strip()
    else:
        errors = os.system('%s %s' % (script_name, infile))
    os.chdir(currdir)
    if errors:
        print "ERROR running %s for %s in directory %s." % (script_name, run_name, rundir)
        return True
    else:
        print "Successfully completed running %s for %s in directory %s." % (script_name, run_name, rundir)
        return False


def RunProcesses(processes, nmultiruns):
    """Runs a list *multiprocessing.Process* processes.

    *processes* is a list of *multiprocessing.Process* objects that
    have not yet been started.

    *nmultiruns* is an integer >= 1 indicating the number of simultaneous
    processes to run.

    Runs the processes in *processes*, making sure to never have more than
    *nmultiruns* running at a time. If any of the processes fail (return
    an exitcode with a boolean value other than *False*), an exception
    is raised immediately. Otherwise, this function finishes when all
    processes have completed.
    """
    if not (nmultiruns >= 1 and isinstance(nmultiruns, int)):
        raise ValueError("nmultiruns must be an integer >= 1")
    processes_started = [False] * len(processes)
    processes_running = [False] * len(processes)
    processes_finished = [False] * len(processes)
    while not all(processes_finished):
        if (processes_running.count(True) < nmultiruns) and not all(processes_started):
            i = processes_started.index(False)
            processes[i].start()
            processes_started[i] = True
            processes_running[i] = True
        for i in range(len(processes)):
            if processes_running[i]:
                if not processes[i].is_alive():
                    processes_running[i] = False
                    processes_finished[i] = True
                    if processes[i].exitcode:
                        raise IOError("One of the processes failed to complete.")
        time.sleep(1)
        

def main():
    """Main body of script."""

    # Current directory
    basedir = os.getcwd()

    # Do we use sbatch when appropriate? Some simple operations will still
    # be run without sbatch even when True
    use_sbatch = True

    # Maximum number of CPUs to try to use at once. If not using sbatch, 
    # don't make this bigger than the number of available cores.
    max_cpus = 100 

    # Density for JPGs converted from PDFs
    jpg_density = 150

    # Specify files holding Sanger sequencing codon mutations and types for 
    # mapmuts_parsesummaryplots.py, as line that can be passed to that script.
    # Put *None* if none specified.
    sanger_codonnmuts = 'Sanger_codonnmuts.txt'
    sanger_codontypes = 'Sanger_codontypes.txt'

    # Specify the replicates, samples, and amplicons.
    replicates = ['replicate_A', 'replicate_B']
    samples = ['WT-1', 'WT-2', 'N334H-1', 'N334H-2']
    amplicons = ['DNA', 'RNA', 'mutDNA', 'virus-p1', 'mutvirus-p1', 'virus-p2', 'mutvirus-p2']

    # Now run mapmuts_makealignments.py. Each sample is run its own subdirectory
    processes = []
    convert_to_jpgs = []
    command_d = {'gzipped':'True',
                 'applyfilter':'True',
                 'minq':'25',
                 'generange':'62 1555',
                 'a1file':'%s/R1_trim3.fasta' % basedir,
                 'a2file':'%s/R2_trim3.fasta' % basedir,
                 'maxn':'2',
                 'minoverlap':'30',
                 'maxrm':'1',
                 'maxa1m':'1',
                 'maxa2m':'1',
                 'maxgenem':'6',
                 'upcase':'test',
                 'write_unaligned':'True',
                }
    for replicate in replicates:
        if not os.path.isdir(replicate):
            os.mkdir(replicate)
        for sample in samples:
            if not os.path.isdir("%s/%s" % (replicate, sample)):
                os.mkdir("%s/%s" % (replicate, sample))
            if 'WT' in sample:
                command_d['fullgenefile'] = '%s/Aichi68-NP_amplicon.fasta' % basedir
            elif 'N334H' in sample:
                command_d['fullgenefile'] = '%s/Aichi68-NP-N334H_amplicon.fasta' % basedir
            else:
                raise ValueError("Cannot assign fullgenefile to sample %s" % sample)
            for amplicon in amplicons:
                subdir = '%s/%s/%s' % (replicate, sample, amplicon)
                command_d['r1files'] = '%s/FASTQ_files/%s/Sample_%s_%s/*R1*.fastq.gz' % (basedir, replicate, sample, amplicon)
                command_d['r2files'] = '%s/FASTQ_files/%s/Sample_%s_%s/*R2*.fastq.gz' % (basedir, replicate, sample, amplicon)
                command_d['outfileprefix'] = subdir.replace('/', '_')
                command_d['samplename'] = subdir.replace('/', ', ')
                processes.append(multiprocessing.Process(target=RunScript,\
                    args=(subdir, "makealignments", 'mapmuts_makealignments.py', list(command_d.items()), use_sbatch, 1))) # use *list* to make a copy of *commands* sinces lists are mutable on subsequent runs
                convert_to_jpgs.append('%s/%s_alignmentstatistics.pdf' % (subdir, command_d['outfileprefix']))
                convert_to_jpgs.append('%s/%s_insertlengths.pdf' % (subdir, command_d['outfileprefix']))
    RunProcesses(processes, nmultiruns=max_cpus)
    for f in convert_to_jpgs:
        assert os.path.isfile(f), "Cannot find file %s" % f
        os.system('convert -density %d %s %s.jpg' % (jpg_density, f, os.path.splitext(f)[0]))

    # Now run mapmuts_alignmentsummaryplot.py to make a summary for each replicate
    convert_to_jpgs = []
    for replicate in replicates:
        commands = []
        for sample in samples:
            for amplicon in amplicons:
                commands.append(("%s/%s/%s/%s/%s_%s_%s_alignmentstatistics.txt" % (basedir, replicate, sample, amplicon, replicate, sample, amplicon), "%s %s" % (sample, amplicon)))
        plotfile = "%s/%s/alignmentsummaryplot.pdf" % (basedir, replicate)
        convert_to_jpgs.append(plotfile)
        commands.append(('plotfile', plotfile))
        RunScript(replicate, 'alignmentsummaryplot', 'mapmuts_alignmentsummaryplot.py', commands, False, 1)
    for f in convert_to_jpgs:
        assert os.path.isfile(f), "Cannot find file %s" % f
        os.system('convert -density %d %s %s.jpg' % (jpg_density, f, os.path.splitext(f)[0]))

    # Now run mapmuts_parsecounts.py to parse counts from alignments.
    processes = []
    convert_to_jpgs = []
    command_d = {'generange':'62 1555',
                 'upcase':'test',
                 'r1exclude':'None',
                 'r2exclude':'None',
                }
    for replicate in replicates:
        for sample in samples:
            if 'WT' in sample:
                command_d['fullgenefile'] = '%s/Aichi68-NP_amplicon.fasta' % basedir
            elif 'N334H' in sample:
                command_d['fullgenefile'] = '%s/Aichi68-NP-N334H_amplicon.fasta' % basedir
            else:
                raise ValueError("Cannot assign fullgenefile to sample %s" % sample)
            for amplicon in amplicons:
                subdir = '%s/%s/%s' % (replicate, sample, amplicon)
                command_d['outfileprefix'] = subdir.replace('/', '_')
                command_d['samplename'] = subdir.replace('/', ', ')
                command_d['alignmentfile'] = '%s/%s/%s_alignments.txt.gz' % (basedir, subdir, command_d['outfileprefix'])
                processes.append(multiprocessing.Process(target=RunScript,\
                    args=(subdir, 'parsecounts', 'mapmuts_parsecounts.py', list(command_d.items()), use_sbatch, 1)))
                convert_to_jpgs.append('%s/%s_codondepth.pdf' % (subdir, command_d['outfileprefix']))
    RunProcesses(processes, nmultiruns=max_cpus)
    for f in convert_to_jpgs:
        assert os.path.isfile(f), "Cannot find file %s" % f
        os.system('convert -density %d %s %s.jpg' % (jpg_density, f, os.path.splitext(f)[0]))

    # Now run mapmuts_parsesummaryplots.py to make a summary for each replicate
    for replicate in replicates:
        commands = []
        for amplicon in amplicons:
            for sample in samples:
                commands.append(("%s/%s/%s/%s/%s_%s_%s" % (basedir, replicate, sample, amplicon, replicate, sample, amplicon), "%s %s" % (sample, amplicon)))
        commands.append(open(sanger_codontypes).read().strip().split(None, 1))
        commands.append(open(sanger_codonnmuts).read().strip().split(None, 1))
        plotfileprefix = "%s/%s/parsesummary" % (basedir, replicate)
        commands.append(('plotfileprefix', plotfileprefix))
        commands.append(('writefracs', 'True'))
        commands.append(('textwritefracs', plotfileprefix))
        commands.append(('pairedcodonplot', 'True'))
        convert_to_jpgs.append("%s_codon_types_and_nmuts.pdf" % plotfileprefix)
        RunScript(replicate, 'parsesummaryplots', 'mapmuts_parsesummaryplots.py', commands, False, 1)
    for f in convert_to_jpgs:
        assert os.path.isfile(f), "Cannot find file %s" % f
        os.system('convert -density %d %s %s.jpg' % (jpg_density, f, os.path.splitext(f)[0]))

    # Now run mapmuts_countparsedmuts.py for each replicate
    convert_to_jpgs = []
    for replicate in replicates:
        commands = [('plotfileprefix', 'countparsedmuts'), ('maxn', '50'), ('legendloc', 'right')]
        for amplicon in ['DNA', 'RNA', 'mutDNA', 'virus-p1', 'virus-p2', 'mutvirus-p1', 'mutvirus-p2']:
            files = ['%s/%s/%s/%s/%s_%s_%s_codoncounts.txt' % (basedir, replicate, sample, amplicon, replicate, sample, amplicon) for sample in samples]
            commands.append((amplicon, ' '.join(files)))
        RunScript(replicate, 'countparsedmuts', 'mapmuts_countparsedmuts.py', commands, False, 1)
        convert_to_jpgs.append('%s/%s/countparsedmuts_multi-nt-codonmutcounts.pdf' % (basedir, replicate))
    for f in convert_to_jpgs:
        assert os.path.isfile(f), "Cannot find file %s" % f
        os.system('convert -density %d %s %s.jpg' % (jpg_density, f, os.path.splitext(f)[0]))

    # Now run mapmuts_inferpreferences.py for each sample in each replicate
    command_d = {'epsilon_concentration':'1.0',
                 'mu_concentration':'1.0',
                 'rho_concentration':'1.0',
                 'pi_concentration':'1.0',
                 'minvalue':'1e-7',
                 'seed':'1',
                 'nruns':'3',
                 'nsteps':'200000',
                 'thin':'200',
                 'stepincrease':'4',
                 'convergence':'0.01',
                 'MCMC_traces':'None',
                 'preference_plots':'None',
                 'ncpus':'6',
                 'removeoutlier':'False',
                }
    processes = []
    for replicate in replicates:
        for sample in samples:
            subdir = "%s/%s/%s" % (basedir, replicate, sample)
            command_d['DNA_files'] = '%s/DNA/%s_%s_DNA_codoncounts.txt' % (subdir, replicate, sample)
            command_d['RNA_files'] = '%s/RNA/%s_%s_RNA_codoncounts.txt' % (subdir, replicate, sample)
            command_d['mutDNA_files'] = '%s/mutDNA/%s_%s_mutDNA_codoncounts.txt' % (subdir, replicate, sample)
            for passage in ['p1', 'p2']:
                command_d['outfileprefix'] = passage
                command_d['mutvirus_files'] = '%s/mutvirus-%s/%s_%s_mutvirus-%s_codoncounts.txt' % (subdir, passage, replicate, sample, passage)
                processes.append(multiprocessing.Process(target=RunScript,\
                    args=(subdir, "%s_inferpreferences" % passage, 'mapmuts_inferpreferences.py', list(command_d.items()), use_sbatch, 6))) # must use *list* to make a copy of *commands* sinces lists are mutable on subsequent runs
    RunProcesses(processes, nmultiruns=max_cpus)

    # Now run mapmuts_preferencemeans.py for each replicate
    for replicate in replicates:
        for p in ['p1', 'p2']:
            preferencefiles = ["%s/%s_equilibriumpreferences.txt" % (sample, p) for sample in samples]
            commands = [('preferencefiles', ' '.join(preferencefiles)), ('outfile', '%s_%s_equilibriumpreferences.txt' % (replicate, p)), ('includestop', 'False')]
            RunScript(replicate, 'preferencemeans', 'mapmuts_preferencemeans.py', commands, False, 1)

    # Now run mapmuts_preferencemeans.py to combine the replicates
    for p in ['p1', 'p2']:
        preferencefiles = ['%s/%s_%s_equilibriumpreferences.txt' % (replicate, replicate, p) for replicate in replicates]
        commands = [('preferencefiles', ' '.join(preferencefiles)), ('outfile', '%s_equilibriumpreferences.txt' % p), ('includestop', 'False')]
        RunScript('./', 'preferencemeans', 'mapmuts_preferencemeans.py', commands, False, 1)

    # Now run mapmuts_preferencescorrelate.py
    correlationdir = './correlations/'
    convert_to_jpgs = []
    if not os.path.isdir(correlationdir):
        os.mkdir(correlationdir)
    # correlations between p1 and p2
    for replicate in replicates:
        commands = [('plotdir', './')]
        commands.append(('samplenames', '%s_p1 %s_p2' % (replicate, replicate)))
        commands.append(('preferencesfiles', '%s/%s/%s_p1_equilibriumpreferences.txt %s/%s/%s_p2_equilibriumpreferences.txt' % (basedir, replicate, replicate, basedir, replicate, replicate)))
        commands.append(('alpha', '0.1'))
        RunScript(correlationdir, 'preferencescorrelate_%s_passages' % replicate, 'mapmuts_preferencescorrelate.py', commands, False, 1)
        convert_to_jpgs.append('%s/%s_p1_vs_%s_p2.pdf' % (correlationdir, replicate, replicate))
    # correlations between replicates for p1
    for i1 in range(len(replicates)):
        replicate1 = replicates[i1]
        for replicate2 in replicates[i1 + 1 : ]:
            commands = [('plotdir', './')]
            commands.append(('samplenames', '%s_p1 %s_p1' % (replicate1, replicate2)))
            commands.append(('preferencesfiles', '%s/%s/%s_p1_equilibriumpreferences.txt %s/%s/%s_p1_equilibriumpreferences.txt' % (basedir, replicate1, replicate1, basedir, replicate2, replicate2)))
            commands.append(('alpha', '0.1'))
            RunScript(correlationdir, 'preferencescorrelate_%s_vs_%s' % (replicate1, replicate2), 'mapmuts_preferencescorrelate.py', commands, False, 1)
            convert_to_jpgs.append('%s/%s_p1_vs_%s_p1.pdf' % (correlationdir, replicate1, replicate2))
    for f in convert_to_jpgs:
        assert os.path.isfile(f), "Cannot find file %s" % f
        os.system('convert -density %d %s %s.jpg' % (jpg_density, f, os.path.splitext(f)[0]))

    # Now run mapmuts_siteprofileplots.py for p1
    logoplot = 'p1_site_preferences_logoplot.pdf'
    if os.path.isfile(logoplot):
        os.remove(logoplot)
    commands = [
            ('sitepreferences', 'p1_equilibriumpreferences.txt'),
            ('outfileprefix', 'p1_'),
            ('siterange', '2 498'),
            ('dsspfile', './DSSP_analysis/2IQH_monomerC.dssp'),
            ('dsspchain', 'C'),
            ('add_rsa', 'True'),
            ('add_ss', 'True'),
            ('nperline', '63'),
            ('includestop', 'False'),
            ]
    RunScript('./', 'siteprofileplots', 'mapmuts_siteprofileplots.py', commands, False, 1)
    assert os.path.isfile(logoplot), "Cannot find file %s" % logoplot
    os.system('convert -density %d %s %s.jpg' % (jpg_density, logoplot, os.path.splitext(logoplot)[0]))



if __name__ == '__main__':
    main() # run the script

