"""Gets all residues within contact with each residue in PDB, using HA1 chain only.

To run this script, open ``pymol`` in the directory where the script is found.
Then at the ``pymol`` command line type::

    run get_PDB_contacts.py

This script was written by Jesse Bloom, 2014."""


def main():
    """Main body of script."""
    # input / output files
    pdbfile = 'PDB_structure/1RVX_trimer_renumbered.pdb'
    pdbchain = 'A'
    contactdistance = 6.0 # in contact if Ca-Ca distance <= this
    contactfile = 'PDB_structure/contacting_residues.txt' # created output file
    residuerange = (1, 365) # residue numbers in HA1, inclusive

    # begin analysis
    cmd.delete('all')
    cmd.load(pdbfile, 'fullPDB')
    cmd.select('HA1_m1', 'chain A and resi %d-%d' % residuerange)
    print 'Finding all residues with CA-CA distance <= %.2f from indicated residue in chain %s of PDB %s looking at residues %d to %d' % (contactdistance, pdbchain, pdbfile, residuerange[0], residuerange[1])
    f = open(contactfile, 'w')
    f.write('# All residues with CA-CA distance <= %.2f from indicated residue in chain %s of PDB %s looking at residues %d to %d\n#REFERENCE_RESIDUE\tCONTACTING_RESIDUES\n' % (contactdistance, pdbchain, pdbfile, residuerange[0], residuerange[1]))
    for r in range(residuerange[0], residuerange[1] + 1):
        cmd.select('CA%d' % r, 'HA1_m1 and resi %d and name CA' % r)
        cmd.select('near%d' % r, '(CA%d around %.2f) and name CA' % (r, contactdistance))
        stored.nearlist = []
        cmd.iterate('near%d' % r, "stored.nearlist.append(resi)")
        contacts = ','.join(stored.nearlist)
        if not contacts:
            contacts = 'None'
        f.write('%d\t%s\n' % (r, contacts))
    f.close()
    print "Script complete."

main() # run the script
