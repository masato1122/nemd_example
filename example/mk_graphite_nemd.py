# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
from ase import Atoms
import ase.io

from nemd.build import build_graphite, write_lammps_data
from nemd.check import check_nemd_structure
from nemd.input import write_nemd_input

def _build_nemd_structure(atoms, nfix=1, nthermo=2, ncenter=10):
    """ Make a structure for NEMD simulation with simply repeating the unit
    structure. The heat flux must be along z-axis.
    
    Parameters
    ------------
    atoms : Atoms object
        rectangular unit cell
    nfix : integer
        number of unit cells at the adiabatic region
    nthermo : integer
        number of unit cells at the thermostats
    ncenter : integer
        number of unit cells at the center region
    
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
    atoms_nemd.cell[2,2] = atoms.cell[2,2] * (nfix*2 + nthermo*2 + ncenter)
    
    ## set atomic positions
    trans = np.zeros(3)
    trans[2] -= atoms.cell[2,2]
    index_nemd = {}
    names = ['bottom', 'hot', 'middle', 'cold', 'top']
    ncells = [nfix, nthermo, ncenter, nthermo, nfix]
    count = 0
    for j in range(len(names)):
        index_nemd[names[j]] = [count, count + natoms * ncells[j] - 1]
        count = index_nemd[names[j]][1] + 1
        for ic in range(ncells[j]):
            trans[2] += atoms.cell[2,2]
            for ia in range(atoms.get_number_of_atoms()):
                atoms_nemd.append(atoms[ia])
                atoms_nemd.positions[-1] += trans
    return index_nemd, atoms_nemd

def write_lmp_input4nemd(structure, index, datafile='data.lammps'):
    """ Make scripts for NEMD simulation
    Unit 
    -------
    temperature : K
    time : ps
    """
    
    tave = 300.
    tcold = 290.
    thot = 310.
    
    time_npt = 10.
    time_nemd = 10.

    ## make a script for NPT simulation
    restartfile0 = 'restart0.nemd'
    write_nemd_input(structure, index, type='npt',
            md_time=time_npt,
            read_data=datafile, output='nemd0.in',
            dumpfile='npt.dump', restartfile=restartfile0,
            taverage=tave, tcold=tcold, thot=thot
            )
    
    restartfile1 = 'restart1.nemd'
    write_nemd_input(structure, index, type='nemd',
            md_time=time_nemd,
            read_data=restartfile0, output='nemd1.in',
            dumpfile='nemd1.dump', restartfile='restart1.nemd',
            taverage=tave, tcold=tcold, thot=thot
            )

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
    datafile = 'data.lammps'
    write_lammps_data(datafile, structure_nemd)

    ## write NEMD input file
    write_lmp_input4nemd(
            structure_nemd, index_nemd, 
            datafile=datafile)

    

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

