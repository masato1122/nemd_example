import numpy as np
from tips.initialize import (set_matplot, set_axis, set_legend)

EvToJ = 1.60218e-19     # J/eV

def plot_energy(
        ax, df_log, lw=0.3, ms=0.6, 
        labels = ['f_t1', 'f_t2'],
        cross_section=None
        ):
    """ Return heat flux [W]
    """
    alpha_y = 1e-3
    
    markers = ['+', '_']
    heats = []
    for i, lab in enumerate(labels):
        xdat = df_log['Time'].values        # ns
        ydat = df_log[lab].values           # eV
        ax.plot(xdat, 
                ydat * alpha_y,
                linestyle='None', lw=lw, marker=markers[i],
                c='black', markersize=ms)

        ##
        n = len(xdat)
        params = np.polyfit(xdat, ydat, 1)
        heat = params[0] * EvToJ * 1e9      # J/s = W
        text = "%.3e W"%(heat)
        ax.text(xdat[int(n*0.5)], 
                ydat[int(n*0.4)] * alpha_y, 
                text,
                fontsize=5, transform=ax.transData,
                horizontalalignment="left", verticalalignment="center")
        heats.append(heat)
    
    ## analyze and print heat flux
    heat = abs((heats[0] - heats[1]) * 0.5)   # Watt
    flux = heat * 1e20 / cross_section        # W/m^2
    print(" HEAT FLUX : %.3e W/m^2"%(flux))
    text = "q = %.3e W/m$^2$"%(flux)
    ax.text(0.5, 0.5, text,
            fontsize=6, transform=ax.transAxes,
            horizontalalignment="center", verticalalignment="center")
    
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Energy (keV)")
    set_axis(ax)
    return flux


