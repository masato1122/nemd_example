""" Make FeCl3-intercalated graphite or graphite structure.

How To Use
-------------
>>> 
>>> ### Stractural Parameters
>>> stage=2
>>> direction="out"; la=5; lb=5; lc=20; fn=data_out.lammps
>>> 
>>> ### NEMD parameters
>>> thermo=langevin
>>> tnpt=1000; tnemd=20000; n=4; m=3; ncenter=5
>>> 
>>> python mknemd_fecl3-gic.py \
>>>     --datafile $fn \
>>>     --stage $stage \
>>>     --direction $direction \
>>>     --la $la --lb $lb --lc $lc \
>>>     --thot 310 --tcold 290 \
>>>     --time_npt $tnpt \
>>>     --time_nemd $tnemd \
>>>     --nloop 3 \
>>>     --tdamp 0.5 \
>>>     --pdamp 1.0 \
>>>     --thermostat ${thermo}
>>> 

"""
import math
import numpy as np
from optparse import OptionParser

from nemd.structure import (
        get_FeCl3_intercalated_graphite,
        get_FeCl3_structure,
        get_indices_at_layers,
        get_ordered_structure
        )
from nemd.build import build_graphite

from nemd import write_lammps_data, write_nemd_inputs
import ase.io
from ase.build import make_supercell

def output_nemd_structure(outfile, nemd, lthermo=0.1, iaxis=2):
    pzs = nemd.get_positions()[:,iaxis]
    lz = nemd.cell[iaxis,iaxis]
    lmiddle = 0.5 - lthermo
    atoms = nemd.copy()
    for ia in range(len(nemd)):
        zi = nemd.get_positions()[ia,iaxis]
        if 0. <= zi <= lz*lmiddle:
            sym = "C"
        elif lz*lmiddle < zi <= lz*0.5:
            sym = "O"
        elif lz*0.5 < zi <= lz*(0.5+lmiddle):
            sym = "C"
        else:
            sym = "Si"
        atoms[ia].symbol = sym
    ##
    ase.io.write(outfile, atoms, vasp5=True, direct=True, sort=True)
    print(" Output", outfile) 

def switch_axes(atoms, axis0, axis1):
    """ switch axis0 and axis1
    """
    iax = {'x':0, 'y':1, 'z':2}
    iax0 = iax[axis0]
    iax1 = iax[axis1]
    
    ## cell size and positions
    natoms = len(atoms)
    cell_new = np.zeros((3,3))
    positions = np.zeros((natoms,3))
    for i0 in range(3):
        if i0 == iax0:
            cell_new[i0,i0] = atoms.cell[iax1,iax1]
            positions[:,i0] = atoms.get_positions()[:,iax1]
        elif i0 == iax1:
            cell_new[i0,i0] = atoms.cell[iax0,iax0]
            positions[:,i0] = atoms.get_positions()[:,iax0]
        else:
            cell_new[i0,i0] = atoms.cell[i0,i0]
            positions[:,i0] = atoms.get_positions()[:,i0]
    
    ##
    new = atoms.copy()
    new.set_cell(cell_new)
    new.set_positions(positions)
    return new

def reget_tags(atoms, iax=2, ntypes_C=4):
    
    ## get layer information
    zlayers, idx_layers = get_indices_at_layers(atoms, iax_out=iax, gap=1.5)
    
    if "Fe" in list(atoms.get_chemical_symbols()):
        tag_list = {"Fe":1, "Cl":2, "C":3}
        ## set tags
        for iat, sym in enumerate(list(atoms.get_chemical_symbols())):
            atoms[iat].tag = tag_list[sym]
        
        ## get indices of carbon layers
        nlayers = len(idx_layers)
        il_carbon = []
        for il in range(nlayers):
            iat = idx_layers[il][0]
            if atoms[iat].symbol == 'C':
                il_carbon.append(il)
        
        ## set tag for the carbon layers
        for ii, il in enumerate(il_carbon):
            for iat in idx_layers[il]:
                atoms[iat].tag = ii%4 + 3
        return atoms
    else:
        ## for graphite
        tag_new = np.zeros(len(atoms), dtype=int)
        for il, indices in enumerate(idx_layers):
            tag_new[indices] = il%ntypes_C + 1
        atoms.set_tags(tag_new)
        return atoms

    return None

def print_box_size(atoms, na, nb, nc):
    print(" === box size (nm) ===")
    for i1 in range(3):
        for j in range(3):
            print("%7.3f"%(atoms.cell[i1,j] * 0.1), end=" ")
        print("")
    print(" (na,nb,nc) = (%d, %d, %d)"%(na, nb, nc))
    print("")

def main(options):
    
    ## create FeCl3-GIC or graphite structure
    ## prim is a rectangular unit
    if options.stage >= 1:
        prim = get_FeCl3_intercalated_graphite(
                nglayers=options.stage, 
                distance=options.distance,
                ncells=[1,1,1])
        charge = True
    elif options.stage == 0:
        prim = build_graphite(n=1, m=1, distance=3.35)
        charge = False
    else:
        prim = get_FeCl3_structure(shape='rect', layer=False)
        init_charges = np.zeros(len(prim))
        for ia in range(len(prim)):
            if prim[ia].symbol == 'Fe':
                init_charges[ia] = 3.
            else:
                init_charges[ia] = -1.
        prim.set_initial_charges(init_charges)
        charge = True
        
    ## make a supercell
    na = int(math.ceil(10*options.la/prim.cell[0,0]))
    nb = int(math.ceil(10*options.lb/prim.cell[1,1]))
    nc = int(math.ceil(10*options.lc/prim.cell[2,2]))
    P = [[na, 0, 0], [0, nb, 0], [0, 0, nc]]
    nemd = make_supercell(prim, P)
    nemd.wrap()
    
    ## print box size
    print_box_size(nemd, na, nb, nc)
    
    ## switch axes so that the heat flows along z-direction.
    if options.direction == 'in':
        nemd = switch_axes(nemd, 'x', 'z')
        IAXIS_DIR = 0
    else:
        IAXIS_DIR = 2
    
    ### get atom indices 
    nemd = get_ordered_structure(nemd, iax=IAXIS_DIR)
    nemd = reget_tags(nemd, iax=IAXIS_DIR)
    
    ## output LAMMPS scripts
    if options.stage == 0:
        write_lammps_data(options.datafile, images=nemd, charge=False)
    else:
        write_lammps_data(options.datafile, images=nemd, charge=True)
    
    ## output POSCAR file
    outfile = "POSCAR.nemd"
    output_nemd_structure(outfile, nemd, lthermo=options.lthermo)
    
    ## areas for thermostat
    if options.direction == 'out':
        ## thermostat is applied to some layers
        from nemd.group import get_nemd_groupids_MP
        lthermo = None
        index_nemd = get_nemd_groupids_MP(
                nemd, thermo_ratio=options.lthermo, iaxis=IAXIS_DIR)
    else:
        ## thermostat is applied to regions
        lthermo = options.lthermo
        index_nemd = None
    
    ## write NEMD input file
    pdamps = [options.pdamp_xy, options.pdamp_xy, options.pdamp_z]
    write_nemd_inputs(
            nemd, 
            lthermo=lthermo, index_nemd=index_nemd,
            charge=charge,
            datafile=options.datafile,
            time_npt=options.time_npt,
            time_increase=options.time_increase,
            time_nemd=options.time_nemd,
            tcold=options.tcold, thot=options.thot,
            nloop=options.nloop,
            damping_time=options.tdamp,
            pressure_damping=pdamps,
            thermostat=options.thermostat
            )

if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("--datafile", dest="datafile", type="string",
            default="data.lammps", help="name of LAMMPS data file")
    
    parser.add_option("--stage", dest="stage", type="int",
            default=0, 
            help="stage: (>0:FeCl3-GIC, =0:graphite, <0:FeCl3")
    
    parser.add_option("--direction", dest="direction", type="string",
            default='in', help="Direction of the heat flow (in or out)")
    
    parser.add_option("--distance", dest="distance", type="float",
            default=3.0, help=
            "distance between a graphene layer and top or bottom of "\
                    "FeCl3 layer")
    
    ## number of units for NEMD
    ## The unit denotes the minimum rectangular cell of FeCl3-GIC
    parser.add_option("--la", dest="la", type="float",
            default=5, 
            help="length of a in-plane direction (nm). a-direction will be" \
                    "the direction of the heat flow for in-plane analysis.")
    
    parser.add_option("--lb", dest="lb", type="float",
            default=5, help="length of a in-plane direction (nm)")
    
    parser.add_option("--lc", dest="lc", type="float",
            default=5, help="length of the out-of-plane direction (nm)")
    
    parser.add_option("--lthermo", dest="lthermo", type="float",
            default=0.1, help="scaled length of the thermostat")
    
    # --- number of the NEMD loop
    parser.add_option("--nloop", dest="nloop", type="int",
            default=2,
            help="number of loops for NEMD calculations (default: 2)")
    
    ## temperatures for NEMD
    parser.add_option("--thot", dest="thot", type="float",
            default=310., help="hot temperature (default: 300 K)")
    
    parser.add_option("--tcold", dest="tcold", type="float",
            default=290., help="cold temperature (default: 290 K)")
    
    ## damping parameters for temperature and pressure
    parser.add_option("--tdamp", dest="tdamp", type="float",
            default=0.5,
            help="damping time for Nose-Hoover (ps) (default: 0.5)")
    
    parser.add_option("--pdamp_xy", dest="pdamp_xy", type="float",
            default=1.0,
            help="damping time for NPT along xy (ps) (default: 1.0)")
    
    parser.add_option("--pdamp_z", dest="pdamp_z", type="float",
            default=10.0,
            help="damping time for NPT along z (ps) (default: 10.0)")
    
    ## simulation times for NEMD
    parser.add_option("--time_npt", dest="time_npt", type="float",
            default=1000., help="NPT time (default: 1000 ps)")
    
    parser.add_option("--time_increase", dest="time_increase", type="float",
            default=None,
            help="time for incrasing temperature (default: 500 ps)")
    
    parser.add_option("--time_nemd", dest="time_nemd", type="float",
            default=2000., help="NEMD time (default: 2000 ps)")
    
    ## type of thermostat (default: Langevin)
    parser.add_option("--thermostat", dest="thermostat", type="string",
            default="lan",
            help="Nose-Hoover (nose) or Langevin (lan) (Default: Langevin)")
    
    (options, args) = parser.parse_args()
    
    main(options)

