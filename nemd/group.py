# -*- coding: utf-8 -*-
import numpy as np
from nemd.structure import get_indices_at_layers

def get_nemd_groupids_MP(atoms, thermo_ratio=0.1, iaxis=2):
    
    Lsystem = atoms.cell[iaxis,iaxis]
    Lthermo = Lsystem * 0.5 * thermo_ratio
    Lmid = Lsystem * 0.5 * (1. - thermo_ratio)
    if np.min(np.diff(atoms.get_positions()[:,iaxis])) < 0.0:
        print("Error: atoms should be ordered properly.")
        exit()

    zlayers, idx_layers = get_indices_at_layers(atoms, iax_out=iaxis, gap=1.5)
    ##
    nlayers = len(zlayers)
    idx_regions = {}
    for il in range(nlayers):
        if 0. <= zlayers[il] <= Lmid:
            region = 'mid1'
        elif Lmid < zlayers[il] <= Lmid + Lthermo:
            region = 'hot'
        elif Lmid+Lthermo < zlayers[il] <= Lmid+Lthermo+Lmid:
            region = 'mid2'
        else:
            region = 'cold'
        ##
        if region not in idx_regions.keys():
            idx_regions[region] = []
        idx_regions[region].extend(idx_layers[il]) 
    
    ##
    idx_nemd = {}
    for region in idx_regions.keys():
        idx = np.asarray(idx_regions[region], dtype=int)
        idx_nemd[region] = [np.min(idx), np.max(idx)]
    
    return idx_nemd

