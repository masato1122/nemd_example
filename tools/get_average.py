# -*- coding: utf-8 -*-
import os.path
import numpy as np
import pandas as pd
from optparse import OptionParser
from tips.lmp import read_log_file

def read_file(filename, zmod=10.):
    ifs = open(filename, "r")
    lines = ifs.readlines()
    ifs.close()
    #
    tmp = {}
    #
    for line in lines:
        data = line.split()
        if line[0] == "#":
            if len(data) >= 6:
                if ("number of structures" in line or 
                        "numbre of structures" in line):
                    nstruct = int(data[-1])
            if line[1] == "#":
                label = data[1]
                tmp[label] = []
            continue
        ##
        if len(data) < 3:
            continue
        iat = int(data[0])
        z = float(data[1]) * zmod
        temp = float(data[2])
        tmp[label].append([iat, z, temp])
    ##
    for key in tmp.keys():
        tmp[key] = np.asarray(tmp[key])
    return nstruct, tmp

def get_all_data(directory, prefix, irange):
    if len(irange.split(":")) != 2:
        print(" Error in the range of index", irange)
        exit()
    i0 = int(irange.split(":")[0])
    i1 = int(irange.split(":")[1]) + 1
    ns = []
    dall = []
    for i in range(i0,i1):
        filename = "%s/%s%d.txt"%(directory, prefix, i)
        nst, each = read_file(filename)
        ns.append(nst)
        dall.append(each)
    return ns, dall

def get_average_values(ns, dall):
    nall = np.sum(ns)
    keys = dall[0].keys()
    ave = {}
    for i, nn in enumerate(ns):
        for key in keys:
            if i == 0:
                ave[key] = dall[i][key] * nn / nall
            else:
                ave[key] += dall[i][key] * nn / nall
    return nall, ave

def write_average(outfile, nall, ave):
    ofs = open(outfile, 'w')
    ofs.write("# number of structures : %d\n"%(nall))
    ofs.write("# natoms positions[nm] temperature[K]\n")
    for key in ave.keys():
        ofs.write("## %s\n"%(key))
        for data in ave[key]:
            ofs.write("%3d  %15.8f  %10.5f\n"%(int(data[0]), data[1], data[2]))
        ofs.write("\n")
    ofs.close()
    print(" Output", outfile)

def get_log_data(directory='.', prefix='log', irange=None, ilog=0):
    if len(irange.split(":")) != 2:
        print("Error in irange", irange)
        exit()
    i0 = int(irange.split(":")[0])
    i1 = int(irange.split(":")[1])
    df_log = []
    dump = {'Time':[], 'Step':[], 'f_t1':[], 'f_t2':[]}
    for ii in range(i0, i1+1):
        fn = "%s/%s%d.txt"%(directory, prefix, ii)
        df_each = read_log_file(fn, ilog=ilog)
        ##
        for key in dump.keys():
            if ii == i0:
                dump[key].extend(df_each[key].values)
            else:
                v0 = dump[key][-1]
                vi = (df_each[key].values[1:] -
                        df_each[key].values[0] + v0)
                dump[key].extend(vi)
    
    ## convert to dict to DataFrame
    df_all = pd.DataFrame()
    for key in dump.keys():
        df_all[key] = dump[key]
    return df_all

def get_cross_sectional_area(dump='npt.dump', iaxis=2):
    """
    iaxis :
        direction of the heat transport
    """
    import ase.io
    if os.path.exists(dump) == False:
        return None
    atoms = ase.io.read(dump, format='lammps-dump-text', index=-1)
    Across = 1.
    for j in range(3):
        if j != iaxis:
            Across *= atoms.cell[j,j]
    return Across 

def main(options):
    
    ns, dall = get_all_data(
            options.directory, options.prefix, options.irange)
    
    df_log = get_log_data(
            directory=options.directory, 
            prefix='log',
            irange=options.irange,
            )
    nptdump = "%s/%s"%(options.directory, options.nptdump)
    Across = get_cross_sectional_area(dump=nptdump)
    
    nall, ave = get_average_values(ns, dall)
    
    write_average(options.outfile, nall, ave)
    
    from nemd.plot import plot_temperature_profile
    plot_temperature_profile(
            ave, log=df_log, 
            figname=options.figname, 
            plot_layer=False,
            cross_section=Across
            )

if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-d", "--directory", dest="directory", type="string",
                default=".", help="directory")
    
    parser.add_option("-p", "--prefix", dest="prefix", type="string",
                default="temp_atom", help="input file name")
    
    parser.add_option("--nptdump", dest="nptdump", type="string",
                default="npt.dump", 
                help="npt.dump is given to get the cross-sectional area.")
    
    parser.add_option("-i", "--irange", dest="irange", type="string",
            default="1:5", help="range of indices")
    
    parser.add_option("-o", "--outfile", dest="outfile", type="string",
            default="average.txt", help="output file")
    
    parser.add_option("--figname", dest="figname", type="string",
            default="fig_ave.png", help="figure name")
    
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

