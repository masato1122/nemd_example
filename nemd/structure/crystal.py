""" The functions in this script are general functions.
"""
import numpy as np
import spglib
import ase
from ase.build import make_supercell

def get_pubchem_structure(cid=24380):
    """
    Parameter
    -------------
    cid : integer
        PubChem ID

    Return
    ----------
    mol : Atoms object
    """
    c = pcp.Compound.from_cid(cid)
    mol = ase.Atoms()
    natoms = len(c.elements)
    mol.number_of_atoms = len(c.elements)
    for ia in range(natoms):
        el = c.elements[ia]
        coord = np.zeros(3)
        coord[0] = c.record['coords'][0]['conformers'][0]['x'][ia]
        coord[1] = c.record['coords'][0]['conformers'][0]['y'][ia]
        if c.coordinate_type.lower() == "2d":
            coord[2] = 0.
        else:
            coord[2] = c.record['coords'][0]['conformers'][0]['z'][ia]
        mol.append(
                ase.Atom(symbol=el, position=coord)
                )
    mol.wrap()
    return mol

def get_graphene_structure(acc=1.42, distance=3.35):
    """ Create ase.Atoms object for graphene
    """
    cell = np.array(
            [
                [acc*0.5*np.sqrt(3), -acc*1.5, 0.0],
                [acc*0.5*np.sqrt(3), acc*1.5, 0.0],
                [0.0, 0.0, distance]
                ]
            )
    atoms = ase.Atoms(cell=cell, pbc=[True,True,True])
    scaled_positions = np.array([[1./3., 2./3., 0.5], [2./3., 1./3., 0.5]])
    positions = np.dot(scaled_positions, cell)
    atoms.append(ase.Atom(symbol="C", position=positions[0]))
    atoms.append(ase.Atom(symbol="C", position=positions[1]))
    atoms.wrap()
    return atoms

def get_stacking_graphene(na=1, nb=1, nc=1, distance=3.35):
    """ Get ase.Atoms object for an AB stacking graphite
    """
    graphene = get_graphene_structure(distance=distance)

    ## make a supercell along x and y directions
    sc = make_supercell(graphene, [[na,0,0], [0,nb,0], [0,0,1]])
    natoms_layer = len(sc)

    ## make a supercell along z direction
    graphite = ase.Atoms(
            cell=[sc.cell[0], sc.cell[1], sc.cell[2]*nc],
            pbc=[True,True,True]
            )
    elements = ['C', 'C']
    for iz in range(nc):

        ## displacment to form AB stacking
        disp = np.zeros(3)
        if iz %2 == 1:
            disp += (graphene.cell[0,:] + graphene.cell[1,:] * 2.) * 2./3.
        disp += sc.cell[2,:] * iz

        ## add atoms at the layer
        for ia in range(len(sc)):
            graphite.append(ase.Atom(
                symbol=elements[iz%2],
                position=sc.positions[ia,:] + disp,
                tag=iz%2+1
                ))
    ##
    graphite.wrap()
    return graphite

def get_standerdized_cell(atoms_orig):
    cell, scaled_positions, numbers = spglib.refine_cell(atoms_orig)
    return _make_new_atoms(cell, scaled_positions, numbers)

def _make_new_atoms(cell, scaled_positions, numbers, pbc=True, center=True):
    """
    cell : shape=(3,3)
    scaled_positions : shape=(natoms,3)
        scaled positions
    numbers : shape=(natoms)
        IDs of chemical symbol
    """
    # --- chemical symbols in the primitive cell
    symbols = []
    for ia, num in enumerate(numbers):
        symbols.append(ase.data.chemical_symbols[num])

    # --- make the primitive cell
    atoms_new = ase.Atoms(
            symbols, 
            np.dot(scaled_positions, cell), 
            cell=cell,
            pbc=[True, True, True]
            )
    if pbc:
        atoms_new.pbc = pbc
    if center:
        atoms_new.center()
    return atoms_new

