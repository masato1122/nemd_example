# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
from ase import Atoms
import ase.io

from nemd.build import build_graphite, write_lammps_data
from nemd.check import check_nemd_structure

def _build_nemd_structure(atoms, nthermo=2, ncenter=10):
    """ Make a structure for NEMD simulation with simply repeating the unit
    structure. The heat flux must be along z-axis.
    
    Parameters
    ------------
    atoms : Atoms object
    nthermo : integer
        number of repetition at the thermostats
    ncenter : integer
        number of repetition at the center region
    
    Return
    --------
    atoms_nemd : Atoms object
        structure for NEMD simulation
    index_nemd : dictionary of integer array
        index_nemd[name] : range of atomic index in the 'name' region
        name : 'left', 'center', or 'right'
    """
    natoms = atoms.get_number_of_atoms()
    
    ## set cell size
    atoms_nemd = Atoms()
    atoms_nemd.cell = np.copy(atoms.cell)
    atoms_nemd.cell[2,2] = atoms.cell[2,2] * (nthermo*2 + ncenter)
    
    ## set atomic positions
    trans = np.zeros(3)
    trans[2] -= atoms.cell[2,2]
    index_nemd = {}
    names = ['left', 'center', 'right']
    ncells = [nthermo, ncenter, nthermo]
    count = 0
    for j in range(3):
        index_nemd[names[j]] = [count, count + natoms * ncells[j] - 1]
        count = index_nemd[names[j]][1] + 1
        for ic in range(ncells[j]):
            trans[2] += atoms.cell[2,2]
            for ia in range(atoms.get_number_of_atoms()):
                atoms_nemd.append(atoms[ia])
                atoms_nemd.positions[-1] += trans
    return index_nemd, atoms_nemd

def main(options):
    
    ## make a conventional cell of graphite
    graphite = build_graphite(n=options.n, m=options.m, 
            distance=options.distance)
    
    ## make a structure for NEMD simulation
    index_nemd, structure_nemd = \
            _build_nemd_structure(graphite, nthermo=1, ncenter=3)
    
    ## output xyz file to visualize the structure in which thermostats and
    ## center region are separated.
    check_nemd_structure(structure_nemd, index_nemd, filename='nemd_check.xyz')
    
    ## output LAMMPS data file
    write_lammps_data('data.lammps', structure_nemd)

if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-n", dest="n", type="int",
            default=4, help="chiral number 1")
    parser.add_option("-m", dest="m", type="int",
            default=4, help="chiral number 2")
    
    parser.add_option("--distance", dest="distance", type="float",
            default=3.35, help="layer distance [A] (default: 3.35)")

    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

