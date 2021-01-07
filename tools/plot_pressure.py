# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from optparse import OptionParser

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tips.initialize import (set_matplot, set_axis, set_legend)

def get_data(filename, inpt=1):
    ifs = open(filename, 'r')
    lines = ifs.readlines()
    ifs.close()

    nlines = len(lines)
    count = 0
    count2 = 0
    il0 = None; il1 = None
    for il in range(nlines):
        if "Step" in lines[il]:
            count += 1
            if count == inpt + 1:
                il0 = il + 1
        if "Loop" in lines[il]:
            count2 += 1
            if count2 == inpt + 1:
                il1 = il - 1
        if "Time step" in lines[il]:
            data = lines[il].split()
            tstep = float(data[-1])
    ## get lines
    if il0 is None:
        return None
    if il1 is None:
        il1 = nlines - 1
    
    ##
    ##Step Temp E_pair E_mol TotEng Press Volume
    labels = lines[il0-1].split()
    nlabs = len(labels)

    ## get data
    dall = []
    for il in range(il0, il1+1):
        each = lines[il].split()
        dall.append([])
        for dd in each:
            dall[-1].append(float(dd))
    dall = np.asarray(dall)

    ## DataFrame
    df = pd.DataFrame()
    istep = [i for i in range(nlabs) if labels[i] == "Step"][0]
    df['Time'] = dall[:,istep] * tstep * 1e-3    # ns
    for lab in labels:
        ilab = [i for i in range(nlabs) if labels[i] == lab][0]
        df[lab] = dall[:,ilab]
    return df

def plot_data(df, figname='fig_npt.png', 
        dpi=300, fontsize=7, fig_width=2.8, aspect=0.6, lw=0.5, ms=2.0):
    
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    plt.subplots_adjust(left=0.20, bottom=0.17, right=0.98,
        top=0.98, wspace=0, hspace=0)
    
    ax = plt.subplot()
    ax2 = ax.twinx()
    
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('P (bar)')
    ax2.set_ylabel('T (K)')
    
    cmap = plt.get_cmap("tab10")
    columns = list(df.columns)
    lns = []
    count = 0
    markers = ['o', '^', 's', 'v', 'd']
    for lab_target in ['Press', 'Pxx', 'Pyy', 'Pzz']:
        if lab_target in columns:
            ll = ax.plot(df['Time'], df[lab_target], 
                    linestyle='-', lw=lw, marker=markers[count%len(markers)],
                    markersize=ms,
                    mfc='none', mew=lw, label=lab_target)
            lns.append(ll)
            count += 1
    ##
    if 'Temp' in columns:
        ll = ax2.plot(df['Time'], df['Temp'], c='black',
                linestyle='-', lw=lw, marker='x', markersize=ms,
                mfc='none', mew=lw, mec='black', label='T')
        lns.append(ll)
    
    ##
    xmax = max(1.05, np.max(df['Time'].values))
    ax.set_xlim(xmax=xmax)
    
    ax.set_ylim([-10000, 10000])
    ax2.set_ylim([280, 320])
    
    ## set legend
    lns_all = lns[0]
    for i in range(1,len(lns)):
        lns_all += lns[i]
    labs = [l.get_label() for l in lns_all]
    set_axis(ax)
    set_axis(ax2)
    ax.legend(lns_all, labs, loc=0, fontsize=5)
    
    ##
    ax.tick_params(labelright=False, right=False)
    ax2.tick_params(labelleft=False, left=False)

    fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    print(' Output', figname)
    return fig

def main(options):
    df = get_data(options.filename, inpt=options.inpt)
    plot_data(df)
    
if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-f", "--filename", dest="filename", type="string",
            default="log0.txt", help="input file name")
    
    parser.add_option("-i", "--inpt", dest="inpt", type="int",
            default=1, help="inpt >= 1")
    
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)


