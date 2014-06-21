"""Gets surface-exposed sites near Caton et al antigenic sites.

Written by Jesse Bloom."""


import mapmuts.dssp



def main():
    """Main body of script."""

    # files names / variables
    catonsitesfile = 'Caton_H1_HA_antigenic_sites.txt' 
    contactingsitesfile = 'PDB_structure/contacting_residues.txt'
    dsspfile = 'PDB_structure/1RVX_trimer_renumbered.dssp'
    dsspchain = 'A'
    rsa_cutoff = 0.2 # only keep sites with >= this RSA
    outfile = 'nearby_antigenic_sites.txt' # created file

    # begin analysis
    catonsites = [int(line.split()[0]) for line in open(catonsitesfile).readlines() if not (line.isspace() or line[0] == '#')]
    print "Read %d sites from %s." % (len(catonsites), catonsitesfile)
    rsas = dict([(r, d['RSA']) for (r, d) in mapmuts.dssp.ReadDSSP(dsspfile, 'Tien2013', dsspchain).iteritems()])
    print "Read %d RSAs from %s." % (len(rsas), dsspfile)
    contactingsites = dict([(int(line.split()[0]), set([int(r) for r in line.split()[1].split(',') if r != 'None'])) for line in open(contactingsitesfile).readlines() if not (line.isspace() or line[0] == '#')])
    print "Read contacting sites for %d sites from %s" % (len(contactingsites), contactingsitesfile)
    allcontacts = []
    for r in catonsites:
        contacts = [r2 for r2 in contactingsites[r] if rsas[r2] >= rsa_cutoff and r2 not in catonsites]
        allcontacts += contacts
    allcontacts = set(allcontacts)
    print "Found a total of %d new contacting residues (not already in antigenic site) that contact an existing antigenic site and have a RSA >= %.2f" % (len(allcontacts), rsa_cutoff)
    open(outfile, 'w').write(\
            '# Antigenic siites that contact antigenic sites.\n' +\
            '# This is all antigenic sites and sites that contact the antigenic sites (antigenic sites defined by %s, contacts defined by %s) and have an RSA >= %.2f (defined by %s)\n' % (catonsitesfile, contactingsitesfile, rsa_cutoff, dsspfile) +\
            '%s' % '\n'.join([str(r) for r in allcontacts] + [str(r) for r in catonsites])\
            )
    print "The antigenic sites and their surface-expose contacts have been written to %s" % outfile


if __name__ == '__main__':
    main()
