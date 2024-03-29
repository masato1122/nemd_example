# -*- coding: utf-8 -*-
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nemd.mpl.initialize import (set_matplot, set_axis, set_legend)
    
from nemd.structure.gic import get_layer_information

def _set_temperature_profile(fontsize=None, fig_width=None, aspect=None,
        xunit=None, xlabel="Position", ylable="Temperature (K)", nfigs=1):
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    plt.subplots_adjust(left=0.20, bottom=0.17, right=0.98,
        top=0.98, wspace=0, hspace=0.25)
    
    axes = []
    for ifig in range(nfigs):
        axes.append(plt.subplot(nfigs,1,ifig+1))
    
    if xunit == 'nm':
        _xunit = 'nm'
    else:
        _xunit = "$\\AA$"
    axes[0].set_xlabel('Position (%s)'%(_xunit))
    axes[0].set_ylabel('Temperature (K)')
    return fig, axes

def plot_temperature_profile_atom(
        dumpfile, atoms=None, lmpinput='nemd.dump', iaxis=2,
        figname=None, outfile=None,
        xunit='nm', log=None,
        dpi=300, fontsize=7, fig_width=2.8, aspect=0.9, lw=0.3, ms=5.0
        ):
    ##groups4plot=None, 
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

    ##
    groups4plot = []
    for name in ids.keys():
        if name is not 'fix':
            groups4plot.append(name)

    ##
    temps = {}
    for group in groups4plot:
        idx = np.arange(ids[group][0]-1, ids[group][1])
        n = len(idx)
        temps[group] = np.zeros((n, 3))
        temps[group][:,0] = np.ones(n)
        temps[group][:,1] = atoms.get_positions()[idx,iaxis]
        temps[group][:,2] = temp_ave[idx]

    if outfile is not None:
        from nemd.file import write_nemd_temperatures
        write_nemd_temperatures(outfile, temps,
                unit_length=xunit,
                nstructures=nstructures)
    
    ## cross-sectional area
    cross = 1.
    for j in range(3):
        if j != iaxis:
            cross *= atoms.cell[j,j]      # A^2
    ##
    plot_temperature_profile(temps,
            xunit=xunit,
            figname=figname,
            log=log,
            cross_section=cross
            )

def plot_temperature_profile(temps, log=None,
        cross_section=None,
        xunit='nm', 
        figname=None, plot_layer=True,
        dpi=300, fontsize=6, lw=0.6, ms=2.0
        ):
    """
    Parameters
    ---------------
    temps : dict of ndarray
        temps['hot'][n,3] : (number of atoms, position[A], temperature[K])
    log : DataFrame object
    cross_section : unit=[A^2]
    """
    cmap = plt.get_cmap("tab10")
    ##
    if log is not None:
        nfigs = 2
        fig_width = 2.5
        aspect = 1.3
    else:
        nfigs = 1
        fig_width = 2.3
        aspect = 0.9
    
    ## Setting for the figure
    fig, axes = _set_temperature_profile(
            fontsize=fontsize, fig_width=fig_width, aspect=aspect, xunit=xunit,
            nfigs=nfigs)
    
    ## plot energy
    flux = None
    if log is not None:
        from tips.plot import plot_energy
        if cross_section is None:
            print("")
            print(" Error: cross sectional area is not given.")
            print("")
        else:
            ## flux : heat flux [W/m^2]
            flux = plot_energy(axes[1], log, cross_section=cross_section)   # Watt
    
    ## unit
    if xunit == 'nm':
        xmod = 0.1
    else:
        xmod = 1.0
    
    ## plot
    for i, group in enumerate(temps.keys()):
        if group == 'hot':
            col = 'red'
            marker = '^'
        elif group == 'cold':
            col = 'blue'
            marker = 'v'
        elif 'mid' in group:
            col = 'purple'
            marker = 'o'
        else:
            continue
        
        ##
        if 'mid' in group:
            out, xfit, yfit = get_gradient(
                    temps[group][:,1], temps[group][:,2])
            gradient = out[0]
        else:
            gradient = None

        ##
        if plot_layer:
            zpos, idx_layers = get_layer_information(temps[group][:,1], gap=1.5)
            temp_ave = np.zeros(len(zpos))
            for il in range(len(zpos)):
                temp_ave[il] = np.average(temps[group][idx_layers[il],2])
        else:
            zpos = temps[group][:,1]
            temp_ave = temps[group][:,2]

        axes[0].scatter(
                zpos * xmod,
                temp_ave,
                s=ms, marker=marker, linewidth=lw,
                edgecolor=col,
                facecolor='None', label=group)
        
        ##
        if gradient is not None:
            axes[0].plot(
                    xfit * xmod, 
                    yfit, 
                    linestyle='-', lw=lw*3, c=cmap(1), alpha=0.7,
                    marker='None', markersize=ms
                    )
            nfit = len(xfit)
            line = "%.3e K/nm\n"%(gradient / xmod)
            if flux is not None:
                kappa = flux / abs(gradient) / 1e10   # W/m-K
                line += "%.3f W/m-K"%(kappa)
                ##
                print(" Kappa: %.3f W/m-K"%(kappa))
            ##
            axes[0].text(
                    xfit[int(nfit*0.5)] * xmod, 
                    np.min(yfit), 
                    line,
                    fontsize=5, 
                    transform=axes[0].transData,
                    horizontalalignment="center", 
                    verticalalignment="top"
                    )
        ##
        if group == 'hot' or group == 'cold':
            axes[0].text(
                    zpos[0] * xmod,
                    np.average(temp_ave),
                    "%.2f"%(np.average(temp_ave)), 
                    fontsize=5, 
                    transform=axes[0].transData,
                    horizontalalignment="right", 
                    verticalalignment="center"
                    )
                   


    set_axis(axes[0])
    set_legend(axes[0], alpha=0.5, fs=6)
    
    ## save a figure
    if figname is not None:
        fig.savefig(figname, dpi=dpi, bbox_inches='tight')
        print(" Output", figname)
    
    return fig

def get_gradient(xdat, temps, buff=0.1, nfit=11):
    """
    Parameters
    -------------
    xdat : unit=[A]
    temps : unit=[K]
    buff :
        both edges with the length of buffer*(total length) are excluded to get
        the average.
    Return
    -------------
    out : shape=(2)
        result of polyfit
    xfit, yfit : shape=(nfit)
    """
    xmin = np.min(xdat)
    xmax = np.max(xdat)
    length = xmax - xmin
    x0 = xmin + length * buff
    x1 = xmax - length * buff
    n = len(xdat)
    idx_ave = [i for i in range(n) if x0 <= xdat[i] <= x1]
    out = np.polyfit(xdat[idx_ave], temps[idx_ave], 1)
    xfit = np.linspace(x0, x1, nfit)
    yfit = np.polyval(out, xfit)
    return out, xfit, yfit

