"""Labels PDB file with site entropies.

To run this script, open ``pymol`` in the directory where the script is found.
Then at the ``pymol`` command line type::

    run entropy_label_PDB.py

This script was written by Jesse Bloom, 2014."""


import mapmuts.io
import mapmuts.bayesian


def main():
    """Main body of script."""
    # input / output files
    aapreferencesfile = 'average_equilibriumpreferences.txt'
    pdbfile = 'PDB_structure/1RVX_trimer_renumbered.pdb'
    subsets = ['allRBS',
               'conservedRBS',
               'antigenicsites',
               'nearantigenicsites',
               ]
    files = {
            'allRBS':'allRBS_residues.txt',
            'conservedRBS':'receptor_binding_residues.txt',
            'antigenicsites':'Caton_H1_HA_antigenic_sites.txt',
            'nearantigenicsites':'nearby_antigenic_sites.txt',
            }
    imagefile = 'PDB_structure/entropy_colored_structure.png'
    dpi = 350

    # compute entropies
    print "\nComputing entropies from amino-acid preferences in %s..." % aapreferencesfile
    aapreferences = mapmuts.io.ReadEntropyAndEquilFreqs(aapreferencesfile)
    mapmuts.bayesian.PreferencesRemoveStop(aapreferences)
    sites = aapreferences.keys()
    sites.sort()
    entropies = {}
    for r in sites:
        entropies[r] = mapmuts.bayesian.SiteEntropy(dict([(aa, aapreferences[r]['PI_%s' % aa]) for aa in mapmuts.sequtils.AminoAcids()]))

    # load the pdbfile and replace residue b-factors with entropies
    cmd.delete('all')
    cmd.load(pdbfile, 'fullPDB')
    cmd.show('cartoon', 'fullPDB')
    for r in sites:
        cmd.alter("fullPDB and resi %d" % r, "b = %g" % entropies[r])

    # color the entire HA molecule by entropy
    print "Coloring entire HA..."
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.show('surface', 'fullPDB')
    cmd.select(None)
    cmd.spectrum('b', 'blue_red', 'fullPDB')
    cmd.color('gray', 'chain B+C')
    # make image
    cmd.set_view (\
        '0.856079578,    0.142146140,   -0.496918648,' +\
        '0.514902651,   -0.151199639,    0.843809366,' +\
        '0.044809252,   -0.978229046,   -0.202629298,' +\
        '0.000000000,    0.000000000, -426.549072266,' +\
        '75.983520508,    0.051200867,   14.679933548,' +\
        '340.855773926,  512.242492676,  -20.000000000')
    cmd.ray(1024 * 2, 768 * 2)
    cmd.png(imagefile, dpi=dpi)
    print "Wrote image to %s" % imagefile

    # color the HA subsets
    for subset in subsets:
        cmd.hide('everything')
        cmd.bg_color('white')
        cmd.show('cartoon', 'fullPDB')
        cmd.color('gray', 'fullPDB')
        cmd.select('HA1_m1', 'chain A and resi 1-365')
        cmd.select(None)
        selectedsites = [line.split()[0] for line in open(files[subset]).readlines() if (not line.isspace()) and line[0] != '#']
        cmd.spectrum('b', 'blue_red', 'chain A and resi %s' % '+'.join([r for r in selectedsites]))
        cmd.show('spheres', 'chain A and resi %s' % '+'.join([r for r in selectedsites]))
        # make image
        cmd.set_view (\
            '0.856079578,    0.142146140,   -0.496918648,' +\
            '0.514902651,   -0.151199639,    0.843809366,' +\
            '0.044809252,   -0.978229046,   -0.202629298,' +\
            '0.000000000,    0.000000000, -426.549072266,' +\
            '75.983520508,    0.051200867,   14.679933548,' +\
            '340.855773926,  512.242492676,  -20.000000000')
        cmd.ray(1024 * 2, 768 * 2)
        ifile = 'PDB_structure/%s_entropy_colored_structure.png' % subset
        cmd.png(ifile, dpi=dpi)
        print "Wrote image to %s" % ifile

    print "Script complete"


main() # run the script
