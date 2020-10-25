# -*- coding: utf-8 -*-
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nemd.mpl.initialize import (set_matplot, set_axis, set_legend)

def _set_temperature_profile(fontsize=None, fig_width=None, aspect=None,
        xunit=None, xlabel="Position", ylable="Temperature (K)"):
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    plt.subplots_adjust(left=0.20, bottom=0.17, right=0.98,
        top=0.98, wspace=0, hspace=0)
    
    ax = plt.subplot()
    if xunit == 'nm':
        _xunit = 'nm'
    else:
        _xunit = "$\\AA$"
    ax.set_xlabel('Position (%s)'%(_xunit))
    ax.set_ylabel('Temperature (K)')
    return fig, ax

def plot_temperature_profile_atom(
        dumpfile, atoms=None, lmpinput='nemd.dump', iaxis=2,
        figname=None, outfile=None,
        plot_layer=True, thd_layer=1.0,
        groups4plot=['hot', 'middle', 'cold'], xunit='nm',
        dpi=300, fontsize=7, fig_width=3.3, aspect=0.9, lw=0.5, ms=7.0):
    
    """ Plot temperature for each atom or segment. Use original atomic positions
    (atoms) and temperatures in "dumpfile". Note that LAMMPS use angstrom as the
    unit of length.
    """
    ## Read temperatures
    from nemd.file import get_averaged_temperatures_atom
    temp_ave, nstructures = \
            get_averaged_temperatures_atom(dumpfile, natoms=len(atoms))
     
    ## Read group IDs
    from nemd.file import read_group_ids
    ids = read_group_ids(lmpinput)
    
    ## Setting for the figure
    fig, ax = _set_temperature_profile(
            fontsize=fontsize, fig_width=fig_width, aspect=aspect, xunit=xunit)
    
    ## unit
    if xunit == 'nm':
        xmod = 0.1
    else:
        xmod = 1.0
    
    ## plot
    colors = ['red', 'green', 'blue']
    for i, group in enumerate(groups4plot):
        i0 = ids[group][0] - 1
        i1 = ids[group][1]
        ax.scatter(
                atoms.positions[i0:i1,iaxis] * xmod, 
                temp_ave[i0:i1],
                s=ms, marker='o', linewidth=lw,
                edgecolor=colors[i], facecolor='None', label=group)
    
    set_axis(ax)
    set_legend(ax, alpha=0.5)
    
    ## save a figure
    if figname is not None:
        fig.savefig(figname, dpi=dpi, bbox_inches='tight')
        print(" Output", figname)
    
    ## output text file
    if outfile is not None:
        from nemd.file import write_nemd_temperatures
        
        ## get data
        temperatures = {}
        for i, group in enumerate(ids):
            i0 = ids[group][0] - 1
            i1 = ids[group][1]
            temperatures[group] = {}
            temperatures[group]['natoms'] = np.ones(i1-i0)
            temperatures[group]['positions'] = atoms.positions[i0:i1,iaxis]*xmod
            temperatures[group]['temperatures'] = temp_ave[i0:i1]
        
        ## output a file
        write_nemd_temperatures(outfile, temperatures, unit_length=xunit,
                nstructures=nstructures)
    return fig

