# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
import pubchempy as pcp
import ase, ase.io
from ase.build.supercells import make_supercell

## pymatgen
from pymatgen import MPRester
from pymatgen.symmetry.bandstructure import HighSymmKpath
import spglib

from nemd.structure.crystal import (
        get_pubchem_structure, 
        get_graphene_structure,
        get_stacking_graphene)

def get_FeCl3mols_intercalated_graphite(
        nprim=[2,2,1], ncells=[1,1,1],
        rectangular=True, iaxial=2,
        layer_distance=3.35, mol_distance=3.0): 
    """
    Parameters
    --------------
    nprim : array, shape=(3)
        number of primitive unit cells of graphene along its translational
        vectors
    ncells : array, shape=(3)
        number of unit cells
    rectangular : bool
        False returns the primitive unit cell of the intercalated structure,
        while True returns the rectangular cell.
    layer_distance : float
        distance between graphene layers
    mol_distance : float
        distance between FeCl3-layer and graphene layer
    """
    fecl3 = get_pubchem_structure(cid=24380)
    intercalated = get_gic_structure(
            nprim=nprim, molecule=fecl3,
            layer_distance=layer_distance, distance=mol_distance
            )
    if rectangular:
        intercalated = get_rectangular_structure(intercalated)
    
    ## make a supercell
    if ncells:
        intercalated = make_supercell(
                intercalated, 
                [[ncells[0],0,0], [0,ncells[1],0], [0,0,ncells[2]]]
                )
    atoms_new = get_ordered_structure(intercalated, iax=iaxial)
    
    return atoms_new

def _get_center_of_a_hexagon(graphene, ibase=0):
    """ Get the center position of a hexagonal lattice in a graphene layer
    Parameter
    -----------
    graphene : ase.Atoms object
        graphene structure
    Return 
    ----------
    center_hext : numpy array, shape=(3)
    """

    ## get ID of an atom at the last layer
    natoms = len(graphene)
    
    ## get IDs of adjacent atoms and vectors to them
    ds = graphene.get_all_distances(mic=True)[ibase]
    idx = np.where((ds < 2.0) & (ds > 1.0))[0]
    if len(idx) < 2:
        print(" Error", idx)
        exit()
    vectors = graphene.get_distances(ibase, idx, vector=True)
    
    ## get a center position of a hexagonal lattice
    r0 = graphene[ibase].position
    center_hex = graphene[idx[0]].position + graphene[idx[1]].position - r0
    return center_hex

def get_gic_structure(
        nprim=[2,2,1], molecule=None, distance=3.0,
        layer_distance=3.35, set_nemdtag=True):
    """ Create a graphene intercalation compund
    Parameters
    -------------
    distance : float, unit=[A]
        distance between the molecule and graphite layer
    """
    ## make a graphtie structure
    graphite = get_stacking_graphene(
            na=nprim[0], nb=nprim[1], nc=nprim[2], distance=layer_distance)
    
    ## center of a hexagonal lattice in the top layer
    natoms_layer = int(len(graphite) / nprim[2])
    center_hex = _get_center_of_a_hexagon(graphite[-natoms_layer:])
    
    ## get the center of the molecule
    center_mol = np.zeros(3)
    for j in range(3):
        center_mol[j] += np.average(molecule.get_positions()[:,j])
    
    ## put a molecule on the center of the hexagonal lattice
    ## Be careful for the z-position
    disp = center_hex + np.array([0, 0, distance]) - center_mol
    intercalate = graphite.copy()
    intercalate.cell[2,2] += distance * 2.0 - layer_distance
    for ia in range(len(molecule)):
        intercalate.append(
                ase.Atom(
                    symbol=molecule[ia].symbol,
                    position=molecule[ia].position + disp
                    )
                )
    ## If the number of graphene layers is even, the position of the next layers
    ## should be adjusted.
    if nprim[2]%2 == 0:
        intercalate = _superpose_next_intercalated_layer(
                intercalate, natoms_layer=natoms_layer)

    ##
    _set_element_center(intercalate, element="Fe")
    
    ## set tag for NEMD simulation
    if set_nemdtag:
        set_tags4md_structure(intercalate)
    return intercalate

def _superpose_next_intercalated_layer(block, natoms_layer=None):
    """
    Parameters
    -----------
    block : ase.Atoms object
        structure on which a molecule is put on graphite
    natoms_layer : integer
        number of C atoms in a layer
    """
    ## vector for the displacement of the next block
    ibase = 0
    center_hex = _get_center_of_a_hexagon(block[:natoms_layer], ibase=ibase)
    disp = center_hex - block[0].position + block.cell[2,:]
    
    ## prepare the new cell
    block2 = block.copy()
    block2.cell[2,2] *= 2
    
    ## displacement vector for the next block
    disp = block[natoms_layer].position - block[0].position
    disp[2] = block.cell[2,2]

    ## superpose the next block
    for ia in range(len(block)):
        block2.append(ase.Atom(
            symbol=block[ia].symbol,
            position=block[ia].position + disp
            ))
    return block2

def get_rectangular_structure(atoms, order_adjust=True, iax=2, eps=1e-7):
    
    #print(atoms.cell)
    P = [[1,1,0], [-1,1,0], [0,0,1]]
    atoms_rec = make_supercell(atoms, P)
    
    ## check cell shape
    for i1 in range(3):
        for i2 in range(3):
            if i1 == i2:
                continue
            if abs(atoms_rec.cell[i1,i2]) > eps:
                print("Error", atoms_rec.cell[i1,i2])
                exit()
            else:
                atoms_rec.cell[i1,i2] = 0.
    ##
    atoms_rec.wrap()

    ##
    if order_adjust:
        atoms_new = get_ordered_structure(atoms_rec, iax=iax)
        return atoms_new
    else:
        return atoms_rec

def get_ordered_structure(atoms, iax=2):
    
    zcoords = atoms.get_positions()[:,iax]
    isort = np.argsort(zcoords)
    atoms_new = ase.Atoms(cell=atoms.cell, pbc=[True, True, True])
    atoms_new = atoms[isort]
    return atoms_new

def _set_element_center(atoms, element='Fe'):
    """ Adjust position of the atoms
    """
    ## get center of Fe atoms
    symbols = list(atoms.get_chemical_symbols())
    ifes = [i for i, el in enumerate(symbols) if el == element]
    center_Fe = np.zeros(3)
    for ia in ifes:
        center_Fe += atoms[ia].position / len(ifes)
     
    ## get center of the cell
    center_cell = np.zeros(3)
    for j in range(3):
        center_cell += atoms.get_cell()[j,:] * 0.5
    
    ## set Fe atom to be the center at xy plane
    disp = center_cell - center_Fe
    disp[2] = 0.
    atoms.translate(disp)
    
def set_tags4md_structure(atoms, iax_out=2, gap=1.0):
    """ Add tags to Atoms object.
    Adjacent graphene layers and different molecules have different tags.
    
    tag :
        1 or 2 for carbon layers
        >= 3 for other elements
    
    Parameters
    --------------
    atoms : ase.Atoms object
    iax_out : integer
        axial index for the out-of-plane
    gap : float
        distance to judge the gap of layers    
    """
    ## get tags for each species
    sym_all = atoms.get_chemical_symbols()
    sym_list = sorted(set(sym_all), key=sym_all.index)
    tags_elements = {}
    itag_others = 3
    for ii, el in enumerate(sym_list):
        if el == "C":
            tags_elements[el] = "12"
        else:
            tags_elements[el] = itag_others
            itag_others += 1
    
    ## add tag for each layer
    zlayers, idx_layers = get_indices_at_layers(atoms, iax_out=iax_out, gap=gap)
    
    tags_all = np.zeros(len(atoms))
    count_carbon = 0
    for idx in idx_layers:
        if atoms[idx[0]].symbol == "C":
            tags_all[idx] = count_carbon%2 + 1
            count_carbon += 1
        else:
            for ia in idx:
                tags_all[ia] = tags_elements[atoms[ia].symbol]
    
    atoms.set_tags(tags_all)

def get_indices_at_layers(atoms, iax_out=2, gap=1.5):
    """
    Return
    -----------
    zlayers : array, float, shape=(number_of_layers)
        averaged z-position of the atoms in each layer
    idx_layers : array
        idx_layers[il] : 
            atomic indices at the il-th layer
    """
    zall = atoms.get_positions()[:,iax_out]
    return get_layer_information(zall, gap=gap)

def get_layer_information(zall, gap=1.5):
    isort = np.argsort(zall)
    
    ## initial set
    zbase = zall[isort[0]]
    idx_layers = []
    idx_layers.append([])
    idx_layers[-1].append(isort[0])
    
    ## analyze for all
    for i in range(1, len(isort)):
        ia = isort[i]
        if abs(zall[ia] - zbase) > gap:
            zbase = zall[ia]
            idx_layers.append([])
        else:
            zbase = max(zbase, zall[ia])
        ##
        idx_layers[-1].append(ia)
    
    ## averaged positions of the layers
    nlayers = len(idx_layers)
    zlayers = np.zeros(nlayers)
    for il in range(nlayers):
        zlayers[il] = np.average(zall[idx_layers[il]])
    ##
    return zlayers, idx_layers

