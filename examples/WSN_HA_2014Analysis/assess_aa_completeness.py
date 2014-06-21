"""Assesses how completely amino-acid substitutions are sampled.

This script assumes that you know the fraction of codon
substitutions that have been sampled. Assuming these codon
substitutions are random, the script then computes the
fraction of amino-acid substitutions that are sampled.

Jesse Bloom, 2014."""


import random
import mapmuts.sequtils


def main():
    """Main body of script."""
    fracsampled = 0.85
    naasampled = naatotal = 0 
    codingsequence = ''.join([nt for nt in mapmuts.sequtils.ReadFASTA('WSN-HA-amplicon.txt')[0][1] if nt.upper() == nt])[3 : -3]
    ncodons = len(codingsequence) // 3
    assert ncodons * 3 == len(codingsequence)
    for icodon in range(ncodons):
        codon = codingsequence[3 * icodon : 3 * icodon + 3]
        aa = mapmuts.sequtils.Translate([('head', codon)])[0][1]
        naatotal += 19
        mutcodons = [c for c in mapmuts.sequtils.Codons() if c != codon]
        assert len(mutcodons) == 63
        mutaas = {}
        for c in mutcodons:
            if random.random() > fracsampled:
                continue
            mutaa = mapmuts.sequtils.Translate([('head', c)])[0][1]
            if not mutaa:
                continue # stop codon
            if mutaa == aa:
                continue # synonymous
            mutaas[mutaa] = True
        naasampled += len(mutaas)
    print "The fraction of amino acids sampled is %.3f (%d of %d)" % (naasampled / float(naatotal), naasampled, naatotal)


main() # run the script
