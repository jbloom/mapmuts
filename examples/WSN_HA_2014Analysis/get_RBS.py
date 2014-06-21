"""Gets all residues contacting receptor in PDB.

To run this script, open ``pymol`` in the directory where the script is found.
Then at the ``pymol`` command line type::

    run get_RBS.py

This script was written by Jesse Bloom, 2014."""


def main():
    """Main body of script."""
    # input / output files
    pdbfile1 = 'PDB_structure/1RVX_trimer_renumbered.pdb' # correctly numbered, but lacks receptor
    pdbfile2 = 'PDB_structure/1RVX_trimer.pdb' # not correctly numbered, but has receptor. Should automatically align to pdbfile1
    pdbchain = 'A'
    contactdistance = 5.0 # in contact if distance <= this
    rbsfile = 'allRBS_residues.txt' # created output file

    # begin analysis
    cmd.delete('all')
    cmd.load(pdbfile1, 'renumbered')
    cmd.select('HA1_m1', 'renumbered and chain A')
    cmd.load(pdbfile2, 'withreceptor')
    cmd.select('receptor', 'withreceptor and resn SIA+GAL and chain A')
    print 'Finding all residues with any atoms within <= %.2f angstroms from the substrate.' % contactdistance
    cmd.select('near_receptor', '(byres (receptor around %.2f)) and name CA and HA1_m1' % contactdistance)
    stored.nearlist = []
    cmd.iterate('near_receptor', "stored.nearlist.append(resi)")
    contacts = set(stored.nearlist)
    print "Found %d contacts." % len(contacts)
    contacts = '\n'.join(contacts)
    open(rbsfile, 'w').write('# All residues with any atoms within <= %.2f from the substrate\n%s' % (contactdistance, contacts))
    print "Contacting residues have been written to %s" % rbsfile
    print "Script complete."

main() # run the script
