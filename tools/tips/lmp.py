# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import ase, ase.io

def get_style_of_lammps_data(filename):

    ifs = open(filename, 'r')
    lines = ifs.readlines()
    ifs.close()

    iatoms = None
    for il, line in enumerate(lines):
        data = line.split()
        if len(data) != 0:
            if data[0] == 'Atoms':
                iatoms = il
                break
    if iatoms is None:
        print("")
        print(' Error1 in', filename)
        print(" cannot find Atoms")
        print("")
        exit()
    ##
    ndat = len(lines[iatoms+5].split())
    style = None
    if ndat == 5:
        style = 'atomic'
    elif ndat == 6:
        style = 'charge'
    else:
        print(" Error in", filename)
        print(ndat, lines[iatoms+5])
        exit()
    return style

def read_log_file(filename, ilog=0):
    ifs = open(filename, 'r')
    lines = ifs.readlines()
    ifs.close()

    nlines = len(lines)
    count = 0
    count2 = 0
    il0 = None; il1 = None
    for il in range(nlines):
        if "Step" in lines[il]:
            count += 1
            if count == ilog + 1:
                il0 = il + 1
        if "Loop" in lines[il]:
            count2 += 1
            if count2 == ilog + 1:
                il1 = il - 1
        if "Time step" in lines[il]:
            data = lines[il].split()
            tstep = float(data[-1])
    ## get lines
    if il0 is None:
        return None
    if il1 is None:
        il1 = nlines - 1

    ##
    ##Step Temp E_pair E_mol TotEng Press Volume
    labels = lines[il0-1].split()
    nlabs = len(labels)

    ## get data
    dall = []
    for il in range(il0, il1+1):
        each = lines[il].split()
        dall.append([])
        for dd in each:
            dall[-1].append(float(dd))
    dall = np.asarray(dall)

    ## DataFrame
    df = pd.DataFrame()
    istep = [i for i in range(nlabs) if labels[i] == "Step"][0]
    df['Time'] = dall[:,istep] * tstep * 1e-3    # ns
    for lab in labels:
        ilab = [i for i in range(nlabs) if labels[i] == lab][0]
        df[lab] = dall[:,ilab]
    return df

def read_lammps_data(filename, style='reax', sort_by_id=False, units="metal"):
    from ase.io.lammpsdata import read_lammps_data as read_lammps_data_ase
    atoms = read_lammps_data_ase(filename, style=style,
            sort_by_id=True, units=units)
    set_symbols(atoms)
    return atoms

def read_lammps_dump(filename, nskip=9, wrap=True):
    """
    Parameters
    --------------
    filename : string
        dump file name

    Return
    --------
    trajectory : list of ase.Atoms object
    dfs : list of DataFrame object
    """
    ## get positions
    trajectory = ase.io.read(filename, format='lammps-dump-text', index=':')
    nstructures = len(trajectory)
    for i in range(nstructures):
        trajectory[i].set_pbc(True)
        if wrap:
            trajectory[i].wrap()

    ## get the data
    ifs = open(filename, "r")
    lines = ifs.readlines()
    ifs.close()

    ## get labels
    line = lines[nskip-1]
    data = line.split()
    labels = []
    for i in range(2,len(data)):
        labels.append(data[i])
    nlabels = len(labels)

    ## read temperatures
    dfs = []
    natoms = len(trajectory[0])
    temperatures = np.zeros((nstructures, natoms))
    for ist in range(nstructures):
        istart = ist * (natoms + nskip) + nskip

        ## read temperature of each atom
        tmp = np.zeros((natoms, nlabels))
        for ia in range(natoms):
            data = lines[istart+ia].split()
            for ii in range(nlabels):
                tmp[ia,ii] = float(data[ii])

        ## get data as DataFrame
        dfs.append(pd.DataFrame())
        for ii in range(nlabels):
            dfs[-1][labels[ii]] = tmp[:,ii]
    return trajectory, dfs

