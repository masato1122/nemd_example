# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
import ase.io

def  check_nemd_structure(atoms, index, filename='nemd.xyz'):
    """ make a xyz file of NEMD structure
    Parameters
    ------------
    atoms : Atoms object of NEMD simulation
    index : dictionary of integer array
        index[name] : range of atomic index at the 'name' region
        name : 'left', 'center', or 'right'
    """
    atoms_out = atoms.copy()
    
    ## change chemical symbols depending on the region
    symbols = ['O', 'C', 'Si']
    for j, name in enumerate(['left', 'center', 'right']):
        for ia in range(index[name][0], index[name][1]+1):
            atoms_out[ia].symbol = symbols[j]
    
    ## output file: allowed format are 'POSCAR' or 'xyz'
    if filename == 'POSCAR':
        ase.io.write(filename, atoms_out, format='vasp',
                direct=True, vasp5=True, sort=False)
    else:
        ase.io.write(filename, atoms_out)
    
    print(" Output", filename)


