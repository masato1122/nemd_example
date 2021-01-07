# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser

def read_file(filename):
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
        z = float(data[1])
        temp = float(data[2])
        tmp[label].append([iat, z, temp])
    ##
    for key in tmp.keys():
        tmp[key] = np.asarray(tmp[key])
    return nstruct, tmp

def get_all_data(prefix, irange):
    i0 = int(irange.split(":")[0])
    i1 = int(irange.split(":")[1]) + 1
    ns = []
    dall = []
    for i in range(i0,i1):
        filename = "%s%d.txt"%(prefix, i)
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

def main(options):
    
    ns, dall = get_all_data(options.prefix, options.irange)
    
    nall, ave = get_average_values(ns, dall)
    
    write_average(options.outfile, nall, ave)
    
    from nemd.plot import plot_temperature_profile
    plot_temperature_profile(ave, figname=options.figname, plot_layer=False)

if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-p", "--prefix", dest="prefix", type="string",
                default="temp_atom", help="input file name")
    
    parser.add_option("-i", "--irange", dest="irange", type="string",
            default="1:5", help="range of indices")
    
    parser.add_option("-o", "--outfile", dest="outfile", type="string",
            default="average.txt", help="output file")
    
    parser.add_option("--figname", dest="figname", type="string",
            default="fig_ave.png", help="figure name")
    
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

