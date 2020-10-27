# -*- coding: utf-8 -*-
import numpy as np
from scipy.spatial.transform import Rotation

import ase
import ase.io
from ase.build import make_supercell
import spglib
from nemd.structure.gic import get_indices_at_layers
from nemd.structure.crystal import (
        get_standerdized_cell,
        get_stacking_graphene)

def get_FeCl3_intercalated_graphite(
        nglayers=3, rectangular=True, distance=3.0, tgra=3.35):
    """ Create FeCl3-intercalated graphite
    nglayers :
        number of graphene layers
    distance : unit=[A]
        space between a graphene layer and the top or bottom of FeCl3 layer
        Note that FeCl3 layer is not 2D, but qusi-3D.
    """
    
    ## prepare FeCl3
    fecl3_std = get_FeCl3_structure(shape="standerdized")
    fecl3 = make_supercell(fecl3_std, [[2,0,0], [0,2,0], [0,0,1]])
     
    ## prepare graphene
    graphene = get_stacking_graphene(na=5, nb=5, nc=nglayers, distance=tgra)
    
    ## set axis
    _rotate_in_2D(fecl3)
    _rotate_in_2D(graphene)
    
    ## apply strain to match two cells
    Lgra = np.linalg.norm(graphene.cell[0])
    Lfecl3 = np.linalg.norm(fecl3.cell[0])
    strain = Lgra / Lfecl3 - 1.0
    print(" Mismatch : %.2f%%"%(strain*100.))
    Lnew = (Lgra + Lfecl3) * 0.5
    _apply_hydro_strain(fecl3, (Lnew - Lfecl3)/Lfecl3)
    _apply_hydro_strain(graphene, (Lnew - Lgra)/Lgra)
    
    ## merge the structures
    gic = _get_superposed_structure(
            graphene, fecl3, distance=distance, tgra=tgra)
    
    if rectangular:
        print("AAAAAAAAAAAAAAAAAAAAA rectangular")
        exit()
    else:
        return gic

def _get_superposed_structure(
        gra, mol, distance=3.0, tgra=3.35, iax=2, tol=1e-5):
    
    for j in range(3):
        if j == iax:
            continue
        if np.linalg.norm(gra.cell[j] - mol.cell[j]) > tol:
            print("Error")
            exit()
    
    ## get positions
    zgra_max = np.max(gra.positions[:,iax])
    zmol_min= np.min(mol.positions[:,iax])
    zmol_max= np.max(mol.positions[:,iax])
    tmol = zmol_max - zmol_min
    
    ## translate vector for the molecule
    disp_mol = np.array([0., 0., zgra_max - zmol_min + distance])
    
    ## adjust cell size
    zheight = gra.cell[iax,iax] - tgra + 2.*distance + tmol
    gra.cell[iax,iax] = zheight
    
    ## superpose
    gic = ase.Atoms(cell=gra.cell)
    for atom in gra:
        gic.append(atom)
    for atom in mol:
        gic.append(ase.Atom(
            symbol=atom.symbol,
            position=atom.position + disp_mol
            ))
    return gic

def _apply_hydro_strain(atoms, strain):
    scaled_positions = atoms.get_scaled_positions().copy()
    atoms.cell *= (1. + strain)
    atoms.set_scaled_positions(scaled_positions)

def _rotate_in_2D(atoms, a0=[1,0]):
    """ a translational vector will be matched to x-axis.
    """
    anow = atoms.cell[0,:2]
    c = np.inner(a0, anow) / np.linalg.norm(a0) / np.linalg.norm(anow)
    angle = np.arccos(np.clip(c, -1., 1.))
    atoms.rotate(angle * 180./np.pi, 'z', rotate_cell=True)
    
    ## check
    anow = atoms.cell[0,:2]
    c = np.inner(a0, anow) / np.linalg.norm(a0) / np.linalg.norm(anow)
    angle = np.arccos(np.clip(c, -1, 1))
    if abs(angle) > 1e-5:
        print("Error", angle)
        exit()
    atoms.wrap()

def get_FeCl3_primitive():    
    cell = np.array([
            [6.635793, -2.965941, 0.0], 
            [6.635793, 2.965941, 0.0], 
            [5.310133, 0.0, 4.963169]
           ])
    positions = np.array([
            [6.19367716 ,  0.        ,  1.65432845],
            [12.38804184,  0.        ,  3.30884055],
            [ 7.87584922, -0.74996191,  0.33087462], 
            [ 7.04528877, -1.10827724,  3.44042905], 
            [11.2012265 , -1.85823915,  2.77771694], 
            [11.53643023,  1.10827724,  1.52273995], 
            [10.70586978,  0.74996191,  4.63229438], 
            [7.3804925  ,  1.85823915,  2.18545206]])
    symbols = ['Fe', 'Fe', 'Cl', 'Cl', 'Cl', 'Cl', 'Cl', 'Cl']
    
    atoms_prim = ase.Atoms(cell=cell)
    for ia, el in enumerate(symbols):
        atoms_prim.append(ase.Atom(el, positions[ia]))
    return atoms_prim

def get_FeCl3_structure(shape="standerdized", layer=True):
    """
    Parameters
    ------------
    shape : string
        primitive, standerdized, layer, rectangular
    layer : bool
        If True, only one layer will be returned.
    """
    atoms_prim = get_FeCl3_primitive()
    
    ## primitive cell
    if 'prim' in shape.lower():
        if layer:
            return _extract_layer(atoms_prim)
        else:
            return atoms_prim
    
    ## standerdized cell
    atoms_std = get_standerdized_cell(atoms_prim)
    if 'stand' in shape.lower():
        if layer:
            return _extract_layer(atoms_std)
        else:
            return atoms_std
    
    ## rectangular cell
    P = [[1,0,0], [1,2,0], [0,0,1]]
    atoms_rec = make_supercell(atoms_std, P)
    if 'rect' in shape.lower():
        if layer:
            return _extract_layer(atoms_rec)
        else:
            return atoms_rec

def _extract_layer(atoms, ilayer=1, tol=3.0):
    """
    Parameters
    -------------
    ilayer : integer
        index of layers (0, 1, ...)
    """
    zlayers, idx_layers = get_indices_at_layers(atoms, tol=3.0)
    layer = ase.Atoms(cell=atoms.cell)
    for ia in idx_layers[ilayer]:
        layer.append(atoms[ia])
    return layer


