"""Profiles core alignment computations in main.MakeAlignments.

This script profiles the time taken for the complete alignment
operations on example data, including reading the data from
gzipped FASTQ files, aligning the data, and writing it to an
output file. The process is profiled
with cProfile, and then analyzed to determine which computations
are consuming the most CPU time. 

Written by Jesse Bloom, 2012.
"""


import time
import sys
import cProfile
import pstats
import mapmuts.main


def DoAlignments():
    """Runs mapmuts.main.MakeAlignments() on example data."""
    r1files = ['R1.fastq.gz']
    r2files = ['R2.fastq.gz']
    a1 = mapmuts.sequtils.ReadFASTA('R1_trim3.fasta')[0][1]
    a2 = mapmuts.sequtils.ReadFASTA('R2_trim3.fasta')[0][1]
    minoverlap = 20
    maxa1m = maxa2m = 2
    maxn = 2
    maxrm = 2
    maxgenem = 4
    minq = 25
    gzipped = True
    applyfilter = True
    fullgene = mapmuts.sequtils.ReadFASTA('Aichi68-NP_amplicon.fasta')[0][1]
    generange = (61, 1555)
    outfileprefix = 'profile_makealignments'
    mapmuts.main.MakeAlignments(r1files, r2files, gzipped, applyfilter,\
            minq, fullgene, generange, a1, a2, maxn, minoverlap, maxrm,\
            maxa1m, maxa2m, maxgenem, outfileprefix)


# Main body of script
profilefile = 'pstats-readpairtogenealigner'
print "Running MakeAlignments (profiling to %s)..." % profilefile,
start = time.clock()
sys.stdout.flush()
cProfile.run('DoAlignments()', profilefile)
t = time.clock() - start
print "completed alignments and profiling in %.3f seconds." % t
s = pstats.Stats(profilefile)
s.strip_dirs()
s.sort_stats('time')
print "Here are the profiling stats:"
s.print_stats()
