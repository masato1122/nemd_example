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

def write_lammps_data(filename, images=None, charge=False):
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
        ofs.write("%d  %d"%(ia+1, atom.tag))
        if charge:
            ofs.write(" %10.7f"%(images.get_initial_charges()[ia]))
        for j in range(3):
            ofs.write(" %15.10f" % (atom.position[j]))
        ofs.write("\n")
    ofs.close()
    print(" Output", filename)

def check_nemd_structure(atoms, index, filename='nemd.xyz'):
    """ make a xyz file of NEMD structure
    Parameters
    ------------
    atoms : Atoms object of NEMD simulation
    index : dictionary of integer array or dictionary of list of integer array
        index[name] : range of atomic index at the 'name' region
        name : 'left', 'center', or 'right'
        #
        "index[name][0][0]:index[name][0][1]+1"
    """
    atoms_out = atoms.copy()

    ## change chemical symbols depending on the region
    symbols = ['Ge', 'O', 'C', 'Si', 'Ge']
    if len(symbols) != len(index):
        print('Error')
        exit()
    for j, name in enumerate(index.keys()):
        if type(index[name][0]) is int:
            for ia in range(index[name][0], index[name][1]+1):
                atoms_out[ia].symbol = symbols[j]
        elif type(index[name][0]) is list:
            for idx_each in index[name]:
                for ia in range(idx_each[0], idx_each[1]+1):
                    atoms_out[ia].symbol = symbols[j]
        else:
            print(" Error in group")
            print(index[name], type(index[name][0]))
            exit()

    ## output file: allowed format are 'POSCAR' or 'xyz'
    if filename == 'POSCAR':
        ase.io.write(filename, atoms_out, format='vasp',
                direct=True, vasp5=True, sort=False)
    else:
        ase.io.write(filename, atoms_out)

    print(" Output", filename)

def read_group_ids(filename):
    """ Read group IDs from a LAMMPS input script
    """
    ifs = open(filename, 'r')
    lines = ifs.readlines()
    ifs.close()
    
    group_ids = {}
    for line in lines:
        data = line.split()
        if len(data) < 4:
            continue
        if line[0] == "#":
            continue
        if data[0] == "group" and data[2] == "id":
            ids = data[3].split(':')
            group_ids[data[1]] = [int(ids[0]), int(ids[1])]
    return group_ids

def get_averaged_temperatures_atom(
        filename, natoms=None, nskip=9, label='f_tempave',
        idx_group=None, outfile=None):
    """ Get averaged temperature from LAMMPS dump file. Note that the first
    data of temperature should be skipped.

    Return
    --------
    taverage : float, shape=(natoms)
    """
    temperatures = read_temperatures_from_dump(
            filename, natoms=natoms, nskip=nskip, label=label)
    taverage = np.zeros(natoms)
    for ia in range(natoms):
        taverage[ia] = np.average(temperatures[1:,ia])
    nstructures = len(temperatures) - 1
    if idx_group is not None and outfile is not None:
        write_averaged_temperatures_atom(outfile, taverage,
                nstructures=nstructures, idx_group=idx_group)
    return taverage, nstructures

def read_temperatures_from_dump(
        filename, natoms=None, nskip=9, label='f_tempave'):
    """
    Parameters
    --------------
    filename : string
        dump file name

    Return
    --------
    temperatures : ndarray, shape=(nstructures, natoms)
    """
    ifs = open(filename, "r")
    lines = ifs.readlines()
    nstructures = int(len(lines) / (natoms + nskip))
    ifs.close()
    
    ## read temperatures
    temperatures = np.zeros((nstructures, natoms))
    for ist in range(nstructures):
        istart = ist * (natoms + nskip) + nskip
        
        ## get label
        labels = lines[istart-1].split()
        if label not in labels:
            print(" %s does not exist in"%(label), end=" ")
            print(labels)
            exit()
            return None
        
        idx_temp = labels.index(label) - 2
        
        ## read temperature of each atom
        for ia in range(natoms):
            data = lines[istart+ia].split()
            temperatures[ist,ia] = float(data[idx_temp])
    
    return temperatures

def get_averaged_temperatures_layer(filename, lmpinput=None, 
        atoms=None, nskip=9, label='f_tempave', iaxis=2, outfile=None, gap=1.5):
    """ Atoms in 'atoms' object should be ordered along z-axis.
    atoms[i].number denotes the layer.
    
    Parameters
    --------------
    outfile : string
        output file name. If it is given, averaged positions and temperatures of
        the layers will be written.

    Return
    ---------
    zt_layers : dict of dict of narray
        zt_layers[group][keyword] : shape=(nlayers)
        keyword : 'positions', 'temperatures', or 'natoms'
    """
    from nemd.structure import get_indices_at_layers
    
    ids_group = read_group_ids(lmpinput)

    ## get averaged temperature of the atoms
    tave_atom, nstructures = get_averaged_temperatures_atom(
            filename, natoms=len(atoms), nskip=nskip, label=label)
    
    ## get atom index in each layer for 'group'
    ids_layers = {}
    for group in ids_group.keys():
         
        ia0 = ids_group[group][0] - 1
        ia1 = ids_group[group][1]
        
        z_tmp, idx_tmp = get_indices_at_layers(
                atoms[ia0:ia1], iax_out=iaxis, gap=gap)
        
        ##
        ids_layers[group] = []
        for il, ids in enumerate(idx_tmp):
            ids_layers[group].append([])
            for iatom in ids:
                ids_layers[group][il].append(iatom+ia0)
        
    ## get the averaged positions and temperatures in the layers
    zt_layers = {}
    keys = ['natoms', 'positions', 'temperatures']
    for group in ids_group.keys():
        nlayers = len(ids_layers[group])
        zt_layers[group] = np.zeros((nlayers,3))
        for il in range(nlayers):
            idx = ids_layers[group][il]
            zt_layers[group][il,0] = len(idx)
            zt_layers[group][il,1] = \
                    np.average(atoms.positions[idx,iaxis]) * 0.1
            zt_layers[group][il,2] = \
                    np.average(tave_atom[idx])
        
    ## outputfile
    if outfile is not None:
        write_nemd_temperatures(
                outfile, zt_layers, 
                nstructures=nstructures,
                unit_length='nm'
                )
    return zt_layers

def write_nemd_temperatures(outfile, temp_layers, 
        unit_length='nm', nstructures=None):
    """ Write temperature profile obtained with a NEMD simulation
    temp_layers : dict of ndarray
        temp_layers[group][n,3], (number of atoms, position[A], temperature[K])
    """
    ofs = open(outfile, 'w')
    if nstructures is not None:
        ofs.write("# number of structures : %d\n"%(nstructures))
    ofs.write("# natoms position[%s] temperature[K]\n"%(unit_length))
    for group in temp_layers.keys():
        ofs.write('## %s\n'%(group))
        nlayers = len(temp_layers[group])
        for il in range(nlayers):
            ofs.write("%3d %15.8f  %8.3f\n"%(
                temp_layers[group][il,0],
                temp_layers[group][il,1],
                temp_layers[group][il,2]
                ))
        ofs.write("\n")
    ofs.close()
    print(" Output", outfile)
    

