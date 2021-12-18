# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
#import pymatgen as mg
import ase, ase.build

def build_graphite(n=3, m=4, type='armchair', distance=3.35, C_C=1.42,
        thd_layer=1.0):
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
    thd_layer : float
        minimum distance between layers

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
    layer1.positions[:,1] = 0.25 * distance
    
    ## make the second layer
    layer2 = layer1.copy()
    disp = np.zeros(3)
    disp[0] = 0.
    disp[1] = distance
    disp[2] = C_C
    layer2.translate(np.tile(disp, (len(layer1), 1)))
    
    ## merge two layers
    graphite = layer1
    for ia in range(len(layer2)):
        graphite.append(layer2[ia])

    ## The out-of-plane direction is chnaged from y-axis to z-axis.
    pos_tmp = np.copy(graphite.positions)
    cell_tmp = np.copy(graphite.cell)
    graphite.positions[:,1] = pos_tmp[:,2]
    graphite.positions[:,2] = pos_tmp[:,1]
    graphite.cell[1] = np.array([0., cell_tmp[2,2], 0.])
    graphite.cell[2] = np.array([0., 0., cell_tmp[1,1]])
    
    ## label layers
    tag_cur = 1
    zmax = graphite.positions[0,2]
    for ia in range(len(graphite)):
        if abs(graphite.positions[ia,2] - zmax) < thd_layer:
            pass
        else:
            tag_cur += 1
        
        graphite[ia].tag = tag_cur
        zmax = max(graphite.positions[ia,2], zmax)

    return graphite

