# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mytool.mpl.initialize import (set_matplot, set_axis, set_legend)

def plot_temperature(filename='dump.md', figname=None,
        dpi=300, fontsize=7, fig_width=3.3, aspect=1.0, lw=1.0, ms=2.0):
    
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    plt.subplots_adjust(left=0.20, bottom=0.17, right=0.98,
        top=0.98, wspace=0, hspace=0)
    
    ax = plt.subplot()
    ax.set_xlabel('')
    ax.set_ylabel('')
    
    ax.plot(xdat, ydat, linestyle='None', lw=lw, marker='o', markersize=ms)
    #ax.plot(xdat, ydat, c='#187FC4', marker="o", c='orange',
    #   mew=lw, ms=5, lw=lw)
    #ax.scatter(xdat, ydat, c=zdat, cmap=cm.rainbow, marker="o", s=10)
    
    set_axis(ax)
    
    
    if figname is not None:
        fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    return fig


