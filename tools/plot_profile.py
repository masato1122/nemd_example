# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser

from nemd.plot import plot_temperature_profile_atom

def main(options):
    
    ## read original structure
    from ase.io.lammpsdata import read_lammps_data
    atoms = read_lammps_data(
            options.lmpdata, style='atomic', 
            sort_by_id=True, units='metal')
    
    ## plot temperature profile
    plot_temperature_profile_atom(
            options.lmpdump, atoms=atoms,
            lmpinput=options.lmpinput,
            figname=options.figname,
            outfile=options.outfile
            )
    
if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("--figname", dest="figname", type="string",
            default="fig_profile.png", help="figure name")
    parser.add_option("--outfile", dest="outfile", type="string",
            default="temperatures_atom.txt", help="output file name")
    
    parser.add_option("--lmpinput", dest="lmpinput", type="string",
            default='nemd.in', help="LAMMPS input script")
    parser.add_option("--lmpdata", dest="lmpdata", type="string",
            default='data.lammps', help="LAMMPS data file")
    parser.add_option("--lmpdump", dest="lmpdump", type="string",
            default='nemd.dump', help="dump file name")
    
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

