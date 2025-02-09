#!/bin/env python
import pmx 
from sys import argv
sequence = argv[1]
protein = pmx.Chain().create(sequence)
for res in protein.residues:
    res.set_phi( -180, True)
    res.set_psi( 180, True)

prolines = protein.fetch_residues('HIS')
for res in prolines:
    res.set_phi( -75, True)
    res.set_psi( 145, True)
    res.set_omega(180)

protein.write('extended.pdb')