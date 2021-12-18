"""
How to Use
-------------
>>> 
>>> ncells = [1, 1, 3]
>>> 
>>> atoms = get_FeCl3_intercalated_graphite(
>>>        nglayers=2, distance=3.0, ncells=ncells)
>>> 
>>> #ase.io.write("POSCAR", atoms, vasp5=True, direct=True)
>>> 
>>> atoms = get_ordered_structure(atoms, iax=2)
>>> idx4nemd = _get_atom_index_range_for_nemd(atoms, nlayers_thermo=2)
>>> 
>>> ## output LAMMPS scripts
>>> write_lammps_data('data.lammps', images=atoms, charge=True)
>>> write_nemd_inputs(atoms, idx4nemd, charge=True)
>>> 

"""
from nemd.structure.fecl3 import (
        get_FeCl3_intercalated_graphite,
        get_FeCl3_structure
        )
from nemd.structure.gic import (
        get_indices_at_layers,
        get_ordered_structure,
        get_layer_information
        )

