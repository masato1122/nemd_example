# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
from nemd.file import get_averaged_temperatures_layer
from tips.lmp import get_style_of_lammps_data

def main(options):
    
    style = get_style_of_lammps_data(options.lmpdata)

    ## read original structure
    from ase.io.lammpsdata import read_lammps_data
    atoms = read_lammps_data(
            options.lmpdata, style=style,
            sort_by_id=True, units='metal')
    
    ##
    zt_layers = get_averaged_temperatures_layer(
            options.lmpdump, lmpinput=options.lmpinput, atoms=atoms,
            outfile=options.outfile, gap=options.gap)
    
if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("--lmpinput", dest="lmpinput", type="string",
            default='nemd0.in', help="LAMMPS input script")
    
    parser.add_option("--lmpdata", dest="lmpdata", type="string",
            default='data.lammps', help="LAMMPS data file")
    
    parser.add_option("--lmpdump", dest="lmpdump", type="string",
            default='nemd1.dump', help="dump file name")
    
    parser.add_option("--outfile", dest="outfile", type="string",
            default=None, help="output file name")
    
    parser.add_option("--gap", dest="gap", type="float",
            default=1.5, help="distance to judge a gap between layers")
    
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

