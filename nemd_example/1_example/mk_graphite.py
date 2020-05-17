# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
import pymatgen as mg
#from ase.build import graphene_nanoribbon as build_ribbon
import ase, ase.io, ase.build

def build_graphite(n=3, m=4, type='armchair', distance=3.35, C_C=1.42):
    ''' Make a graphite using graphene_nanoribbon in ASE library
    Parameters
    -------------
    n, m : chiral numbers
    type : string
        direction (default: armchair)
    distance : float
        layer distance (default: 3.35 A)
    C_C : float
        C-C bond length (default: 1.42 A)

    Return
    -----------
    graphite : Atoms object
        graphite structure

    Example
    ------------
    > bukld_graphite(n=3, m=4)
    '''
    ## make the first layer
    layer1 = ase.build.graphene_nanoribbon(
            n, m, type=type, C_C=C_C,
            saturated=False, sheet=True)
    
    ## set the layer distance
    layer1.cell[1,1] = 2. * distance
    
    ## make the second layer
    layer2 = layer1.copy()
    disp = np.zeros(3)
    disp[0] = 0.
    disp[1] = distance
    disp[2] = C_C
    layer2.translate(np.tile(disp, (layer1.get_number_of_atoms(), 1)))
    
    ## merge two layers
    graphite = layer1
    for ia in range(layer2.get_number_of_atoms()):
        graphite.append(layer2[ia])

    print(graphite.cell)
    
    ## The initial out-of-plane direction is y-axis.
    ## It is chnaged to z-axis.
    graphite.rotate(90.0, 'x', rotate_cell=True)
    graphite.rotate('y', 'z', rotate_cell=True)
    
    print(graphite.cell)
    return graphite

    def main(options):
    
    graphite = build_graphite(n=options.n, m=options.m, 
            distance=options.distance)
    
    ase.io.write(
            'POSCAR', graphite, format='vasp', 
            direct=True, vasp5=True, sort=False)
    print(" Output POSCAR")


    
if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("--n", dest="n", type="int",
            default=4, help="chiral number 1")
    parser.add_option("--m", dest="m", type="int",
            default=4, help="chiral number 2")
    
    parser.add_option("--distance", dest="distance", type="float",
            default=3.35, help="layer distance [A] (default: 3.35)")

    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

