# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
from ase import Atoms
import ase.io

from nemd.build import build_graphite
from nemd.file  import check_nemd_structure, write_lammps_data
from nemd.input import write_nemd_inputs

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
    natoms = len(atoms)
    
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
            for ia in range(len(atoms)):
                atoms_nemd.append(atoms[ia])
                atoms_nemd.positions[-1] += trans
    return index_nemd, atoms_nemd

def main(options):
    
    ## make a conventional cell of graphite
    graphite = build_graphite(n=options.n, m=options.m, 
            distance=options.distance)
    
    ## make a structure for NEMD simulation
    index_nemd, structure_nemd = _build_nemd_structure(
            graphite, nthermo=options.nthermo, ncenter=options.ncenter)
    
    ## output xyz file to visualize the structure in which thermostats and
    ## center region are separated.
    check_nemd_structure(structure_nemd, index_nemd, filename='nemd_check.xyz')
    
    ## output LAMMPS data file
    datafile = 'data.lammps'
    write_lammps_data(datafile, structure_nemd)
    
    ## write NEMD input file
    write_nemd_inputs(
            structure_nemd, index_nemd,
            datafile=datafile,
            time_npt=options.time_npt,
            time_nemd=options.time_nemd,
            tcold=options.tcold, thot=options.thot)

if __name__ == '__main__':
    parser = OptionParser()
    
    # --- graphite
    parser.add_option("-n", dest="n", type="int",
            default=4, help="chiral number 1")
    parser.add_option("-m", dest="m", type="int",
            default=4, help="chiral number 2")
    
    # --- structure
    parser.add_option("--nthermo", dest="nthermo", type="int",
            default=1, help="number of unit cells in a thermostat (default: 1)")
    parser.add_option("--ncenter", dest="ncenter", type="int",
            default=3, 
            help="number of unit cells in the center region (default: 3)")
    
    # --- temperature
    parser.add_option("--thot", dest="thot", type="float",
            default=310., help="hot temperature (default: 300 K)")
    parser.add_option("--tcold", dest="tcold", type="float",
            default=290., help="cold temperature (default: 290 K)")
    
    # --- MD time
    parser.add_option("--time_npt", dest="time_npt", type="float",
            default=10., help="NPT time (default: 10 ps)")
    parser.add_option("--time_nemd", dest="time_nemd", type="float",
            default=100., help="NEMD time (default: 100 ps)")
    
    # --- layer distance
    parser.add_option("--distance", dest="distance", type="float",
            default=3.35, help="layer distance [A] (default: 3.35)")

    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

