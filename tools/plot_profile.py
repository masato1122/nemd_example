# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
from nemd.plot import plot_temperature_profile_atom
from tips.lmp import get_style_of_lammps_data, read_log_file

def main(options):
    
    ## read log file
    energy_labels = ['f_t1', 'f_t2']
    df_log = None
    if options.log is not None:
        df_log = read_log_file(options.log, ilog=0)
        columns = list(df_log.columns)
        if energy_labels[0] in columns and energy_labels[1] in columns:
            plot_log = True
        else:
            df_log = None

    ## judge style of data.lammps
    style = get_style_of_lammps_data(options.lmpdata)
    
    ## read original structure
    from ase.io.lammpsdata import read_lammps_data
    atoms = read_lammps_data(
            options.lmpdata, style=style, 
            sort_by_id=True, units='metal')
    
    ## plot temperature profile
    plot_temperature_profile_atom(
            options.lmpdump, atoms=atoms,
            iaxis=options.iaxis,
            lmpinput=options.lmpinput,
            figname=options.figname,
            outfile=options.outfile,
            log=df_log
            )
    
if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("--figname", dest="figname", type="string",
            default="fig_profile.png", help="figure name")
    parser.add_option("--outfile", dest="outfile", type="string",
            default="temperatures_atom.txt", help="output file name")
    
    parser.add_option("--lmpinput", dest="lmpinput", type="string",
            default='nemd.in', help="LAMMPS input script")
    parser.add_option("--log", dest="log", type="string",
            default='log.txt', help="LAMMPS log file")
    parser.add_option("--lmpdata", dest="lmpdata", type="string",
            default='data.lammps', help="LAMMPS data file")
    parser.add_option("--lmpdump", dest="lmpdump", type="string",
            default='nemd.dump', help="dump file name")
    
    parser.add_option("--iaxis", dest="iaxis", type="int",
            default=2, help="axial index (x:0, y:1, z:2)")
    
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

