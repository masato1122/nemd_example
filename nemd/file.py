# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
import pymatgen as mg
import ase, ase.build

def _check_rectangular(cell, symprec=1e-5):
    for j in range(-1,2,1):
        p1 = np.dot(cell[j,:], cell[j+1,:])
        p2 = np.linalg.norm(cell[j]) * np.linalg.norm(cell[j+1])
        angle = np.arccos(p1/p2) * 180./np.pi
        if abs(angle-90.) > symprec:
            return None
    return 0

def write_lammps_data(filename, images=None):
    """ Make a LAMMPS data file 
    filename : string
        output file name of the LAMMPS data file
    images : Atoms object of ASE
        atomic type will be represented by 'tag' of each atom
    """
    ## get the list of tags
    tags_list = list(set(images.get_tags()))
    tags_list.sort()
    
    ## get corresponding list of symbols to the tag list
    symbols_list = {}
    for ia in range(len(images)):
        symbols_list[images[ia].tag] = images[ia].symbol
    
    # --- check if the system is rectangular or not
    if _check_rectangular(images.cell) is None:
        print("Error: the system is not rectangular.")
        exit()
    
    # --- output a file
    natoms = len(images)
    ofs = open(filename, "w")
    
    ofs.write("Position data\n")
    ofs.write("\n")
    ofs.write(" %d atoms\n"%(natoms))
    ofs.write("\n")
    ofs.write(" %d atom types\n" % (len(tags_list)))
    ofs.write("\n")
    for j, direction in enumerate(['x', 'y', 'z']):
        ofs.write(" 0.0000  %15.13f  %slo %shi\n" % (images.cell[j,j],
            direction, direction))
    ofs.write("\n")
    
    # --- masses
    ofs.write("Masses\n")
    ofs.write("\n")
    for itag, tag in enumerate(tags_list):
        sym = symbols_list[tag]
        number = ase.data.atomic_numbers[sym]
        mass_each = ase.data.atomic_masses[number]
        ofs.write(" %2d %13.8f %3s\n" % (tag, mass_each, sym))
    ofs.write("\n")
    
    # --- atoms
    ofs.write("Atoms\n")
    ofs.write("\n")
    for ia, atom in enumerate(images):
        ofs.write("%d  %d  "%(ia+1, atom.tag))
        for j in range(3):
            ofs.write(" %15.10f" % (atom.position[j]))
        ofs.write("\n")
    ofs.close()
    print(" Output", filename)

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
    symbols = ['Ge', 'O', 'C', 'Si', 'Ge']
    if len(symbols) != len(index):
        print('Error')
        exit()
    for j, name in enumerate(index.keys()):
        for ia in range(index[name][0], index[name][1]+1):
            atoms_out[ia].symbol = symbols[j]

    ## output file: allowed format are 'POSCAR' or 'xyz'
    if filename == 'POSCAR':
        ase.io.write(filename, atoms_out, format='vasp',
                direct=True, vasp5=True, sort=False)
    else:
        ase.io.write(filename, atoms_out)

    print(" Output", filename)

