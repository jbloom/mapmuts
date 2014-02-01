"""Profiles core alignment computations in main.ParseNTCodonCounts

This script profiles the time taken for alignment parsing,
including reading the data from gzipped FASTQ files and
processing it. The process is profiled 
with cProfile, and then analyzed to determine which computations
are consuming the most CPU time. 

Written by Jesse Bloom, 2012.
"""


import time
import sys
import cProfile
import pstats
import mapmuts.main


def ParseAlignments():
    """Runs mapmuts.main.ParseNTCodonCounts() on example data."""
    r1files = ['R1.fastq.gz']
    r2files = ['R2.fastq.gz']
    alignmentfile = 'testparse_alignments.txt.gz'
    fullgene = mapmuts.sequtils.ReadFASTA('Aichi68-NP_amplicon.fasta')[0][1]
    generange = (61, 1555)
    gene = fullgene[generange[0] : generange[1]].upper()
    r1exclude = [1, 23, 27]
    r2exclude = [32]
    outfileprefix = 'profile_parsealignments'
    logfile = 'profile_parsealignments.log'
    log = open(logfile, 'w')
    mapmuts.main.ParseNTCodonCounts(alignmentfile, outfileprefix, gene,\
            r1exclude, r2exclude, log=log)
    log.close()


# Main body of script
profilefile = 'pstats-profileparsealignments'
print "Running ParseNTCodonCounts (profiling to %s)..." % profilefile,
start = time.clock()
sys.stdout.flush()
cProfile.run('ParseAlignments()', profilefile)
t = time.clock() - start
print "completed parsing and profiling in %.3f seconds." % t
s = pstats.Stats(profilefile)
s.strip_dirs()
s.sort_stats('time')
print "Here are the profiling stats:"
s.print_stats()
