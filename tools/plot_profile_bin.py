"""
Hot To Use
---------------
>>> 
>>> python plot_profile_bin.py \
>>>         -d nemd1.dump \
>>>         -i nemd1.in \
>>>         --lbin 5.0 \
>>>         --direction in
>>> 
"""
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import math
from optparse import OptionParser
import ase.io
from tips.lmp import read_lammps_dump, read_log_file

def read_lmpinput(filename):
    ifs = open(filename, 'r')
    lines = ifs.readlines()
    ifs.close()

    nlines = len(lines)
    regions = {}
    groups = {}
    for il in range(nlines):
        line = lines[il]
        data = line.split()
        if len(data) == 0:
            continue
        if line[0] == "#":
            continue
        if data[0] == 'region':
            z0 = data[7]
            z1 = data[8]
            regions[data[1]] = np.zeros(2)
            if z0 == 'INF':
                regions[data[1]][0] = -1e5
            else:
                regions[data[1]][0] = float(z0)
            ##
            if z1 == 'INF':
                regions[data[1]][1] = 1e5
            else:
                regions[data[1]][1] = float(z1)
        if data[0] == 'group':
            if data[2] == 'region':
                groups[data[3]] = data[1]
    ##
    ranges = {}
    for key in regions.keys():
        ranges[groups[key]] = regions[key]
    return ranges

def get_temperatures_in_each_region(
        coords, temperatures, zrange=[1e-5, 1e5], lbin=0.5, nstructures=None):
    """
    Parameters
    ---------------
    coords, temperatures : np.shape(ndata)
        dumped data at different times
    zrange : 
        range of coordinate
    lbin : unit[A]
        length of each bin
    """
    
    nbins = int(math.ceil((zrange[1]-zrange[0])/lbin))
    leach = (zrange[1] - zrange[0]) / nbins
    
    ##
    df = pd.DataFrame()

    ## set bins
    bins = np.zeros(nbins+1)
    df['position'] = np.zeros(nbins)
    for ib in range(nbins+1):
        bins[ib] = zrange[0] + leach * ib 
        if ib != 0:
            df['position'][ib-1] = (bins[ib-1] + bins[ib]) * 0.5
    
    ## prepare DataFrame
    idx_bins = np.digitize(coords, bins)
    df['temperature'] = np.zeros(nbins)
    df['ndata'] = np.zeros(nbins, dtype=int)
    if nstructures is not None:
        df['natoms'] = np.zeros(nbins)
    
    ## get averaged temperatures for each bin
    for ib in range(nbins):
        idx_atoms = [iat for iat, ii in enumerate(idx_bins) if ii == ib+1]
        df['temperature'].values[ib] = np.average(temperatures[idx_atoms])
        df['ndata'].values[ib] = len(idx_atoms)
        if nstructures is not None:
            df['natoms'].values[ib] = len(idx_atoms) / nstructures
    return df

def get_temperatures_bin(dumpfile, lbin=1.0, regions=None, iax_heat=2,
        label_temp='f_tempave'):
    """ 
    Note that the first data (0th) data in dumpfile should be omitted.
    """
    ## read .dump file
    trajectory, dfs_dump = read_lammps_dump(dumpfile, wrap=True)
    nstructures = len(trajectory) - 1
    
    ## system length
    Lsystem = trajectory[0].cell[iax_heat,iax_heat]
    
    ## prepare dict
    temp_dump = {}
    for name in regions.keys():
        temp_dump[name] = []
    
    ## put all data in an array
    for ist in range(1,nstructures):
        positions = trajectory[ist].get_positions()
        for name in regions.keys():
            z0 = regions[name][0]
            z1 = regions[name][1]
            idx_name = [iat for iat, zz in 
                    enumerate(positions[:,iax_heat]) 
                    if z0 < zz <= z1]
            for ia in idx_name:
                temp_dump[name].append([
                    positions[ia,iax_heat],
                    dfs_dump[ist][label_temp].values[ia]
                        ])
    
    ## conver to ndarray
    for name in regions.keys():
        temp_dump[name] = np.asarray(temp_dump[name])
    
    ## get temperatures in each region
    df_ave = {}
    for name in regions.keys():
        zrange = regions[name].copy()
        zrange[0] = max(zrange[0], np.min(temp_dump[name][:,0]))
        zrange[1] = min(zrange[1], np.max(temp_dump[name][:,0]))
        df_ave[name] = \
                get_temperatures_in_each_region(
                        temp_dump[name][:,0],
                        temp_dump[name][:,1],
                        zrange=zrange,
                        lbin=lbin,
                        nstructures=nstructures
                        )
    return df_ave

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tips.initialize import (set_matplot, set_axis, set_legend)
    
def make_figure(fontsize=None, fig_width=None, aspect=None, nfigs=None):
    
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    plt.subplots_adjust(left=0.20, bottom=0.17, right=0.98,
        top=0.98, wspace=0.10, hspace=0.25)
    axes = []
    for ifig in range(nfigs):
        axes.append(plt.subplot(nfigs,1,ifig+1))
    return fig, axes
    
def plot_temperature_profile(ax, df_ave, lw=0.8, ms=2.0):
    
    ax.set_xlabel('Position (nm)')
    ax.set_ylabel('Temperature (K)')
    
    markers = ['o', 'v', 's', '^', 'd']
    for ii, name in enumerate(df_ave.keys()):
        xdat = df_ave[name]['position'].values * 0.1
        ydat = df_ave[name]['temperature'].values
        if 'hot' in name:
            col = 'red'
        elif 'col' in name:
            col = 'blue'
        elif 'mid' in name:
            col = 'purple'
        ax.plot(xdat, ydat, 
                linestyle='None', 
                lw=lw, marker=markers[ii%len(markers)], 
                markersize=ms, c=col,
                mfc='none', mew=lw,
                label=name)
    
    set_axis(ax)
    set_legend(ax, fs=5, alpha=0.5)

def main(options):
    
    ## read log file
    plot_log = False
    energy_labels = ['f_t1', 'f_t2']
    if options.log is not None:
        df_log = read_log_file(options.log, ilog=0)
        columns = list(df_log.columns)
        if energy_labels[0] in columns and energy_labels[1] in columns:
            plot_log = True
    
    if plot_log:
        nfigs = 2
        fig_width = 2.5
        aspect = 1.3
    else:
        nfigs = 1
        fig_width = 2.3
        aspect = 0.9
    
    ## read lammps input script
    regions = read_lmpinput(options.lmpinput)
    for key in regions.keys():
        print(" %5s : %7.3f to %7.3f"%(key, regions[key][0], regions[key][1]))
    
    ##
    #if options.direction == 'in':
    df_ave = get_temperatures_bin(
            options.lmpdump, lbin=options.lbin, regions=regions)
    #else:
    #read_temperatures_out(
    #        options.lmpdump, regions=regions)
     
    ## output results
    os.makedirs(options.outdir, exist_ok=True)
    for name in df_ave.keys():
        outfile = "%s/%s.csv"%(options.outdir, name)
        df_ave[name].to_csv(outfile, index=False)
    
    ## plot data
    fig, axes = make_figure(
            fontsize=6, fig_width=fig_width, aspect=aspect, nfigs=nfigs)
    plot_temperature_profile(axes[0], df_ave)
    if plot_log:
        
        ## get cross section
        atoms = ase.io.read(
                options.lmpdump, format='lammps-dump-text', index=0)
        cross = atoms.cell[0,0] * atoms.cell[1,1]      # A^2
        
        from tips.plot import plot_energy
        heat = plot_energy(
                axes[1], df_log, 
                labels=energy_labels,
                cross_section=cross)
        
    fig.savefig(options.figname, dpi=300, bbox_inches='tight')
    print(" Output", options.figname)
    return fig

if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("--figname", dest="figname", type="string",
            default="fig_profile.png", help="figure name")
    
    parser.add_option("-d", "--lmpdump", dest="lmpdump", type="string",
            default="nemd1.dump", help=".dump file")
    
    parser.add_option("-i", "--lmpinput", dest="lmpinput", type="string",
            default="nemd1.in", help=".in file name")
    
    parser.add_option("--log", dest="log", type="string",
            default=None, help="log file")
    
    parser.add_option("--lbin", dest="lbin", type="float",
            default=5.0, 
            help="length of each bin [A]. This is used only for in-plane.")
    
    parser.add_option("--direction", dest="direction", type="string",
            default="out", help="in or out")
    
    parser.add_option("--outdir", dest="outdir", type="string",
            default="out_ave", help="output directory")
    
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)


