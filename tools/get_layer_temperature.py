# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser

from nemd.file import get_averaged_temperatures_layer

def main(options):
    
    ## read original structure
    from ase.io.lammpsdata import read_lammps_data
    atoms = read_lammps_data(
            options.lmpdata, style='atomic',
            sort_by_id=True, units='metal')
    
    ##
    zt_layers = get_averaged_temperatures_layer(
            options.lmpdump, lmpinput=options.lmpinput, atoms=atoms,
            outfile=options.outfile)
    
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
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

