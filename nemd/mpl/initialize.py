# -*- coding: utf-8 -*-
#import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.ticker import *
#from matplotlib import cm
import matplotlib.gridspec as gridspec

#---------- Start matplot -------------#
def set_matplot(fontsize=9):
    lw_bor = 0.5
    plt.rcParams["font.size"] = fontsize
    #plt.rcParams['font.family'] = 'Myriad Pro'
    #plt.rcParams["mathtext.fontset"] = 'stix'
    plt.rcParams["mathtext.fontset"] = 'dejavusans'
    plt.rcParams['axes.linewidth'] = lw_bor
    plt.rcParams['xtick.major.width'] = lw_bor
    plt.rcParams['xtick.minor.width'] = lw_bor
    plt.rcParams['ytick.major.width'] = lw_bor
    plt.rcParams['ytick.minor.width'] = lw_bor 

def set_spaces(left=0.14, bottom=0.14, right=0.98, top=0.98, ratio=1.0):
    FIG_WIDTH = 3.3
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH*ratio))
    plt.subplots_adjust(
            left=left, bottom=bottom,
            right=right, top=top, wspace=0, hspace=0)
    return fig, plt

def set_axis(ax, xformat="linear", yformat="linear", 
        xticks=None, mxticks=None, yticks=None, myticks=None,
        labelbottom=None, length=2.4, width=0.5):
    ax.tick_params(axis='both', which='major', 
            direction='in', length=length, width=width)
    ax.tick_params(axis='both', which='minor',
            direction='in', length=length*0.6, width=width)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    
    #--- for linear scale
    if xticks is not None:
        ax.xaxis.set_major_locator(tick.MultipleLocator(xticks))
    if mxticks is not None:
        interval = float(xticks) / float(mxticks)
        ax.xaxis.set_minor_locator(tick.MultipleLocator(interval))
    if yticks is not None:
        ax.yaxis.set_major_locator(tick.MultipleLocator(yticks))
    if myticks is not None:
        interval = float(yticks) / float(myticks)
        ax.yaxis.set_minor_locator(tick.MultipleLocator(interval))
    #--- for logscale
    if xformat.lower() == "log":
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(tick.LogLocator(base=10.0, numticks=15))
    if yformat.lower() == "log":
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(tick.LogLocator(base=10.0, numticks=15))
    return ax

def get_both_axis(erange, ylabel, ksym, klabels, x2label):
    gs = gridspec.GridSpec(1,3)
    ax1 = plt.subplot(gs[0,0:2])
    ax2 = plt.subplot(gs[0,-1])
    set_axis(ax1)
    set_axis(ax2)
    if ksym is not None:
        ax1.set_xticks(ksym)
    if klabels is not None:
        ax1.set_xticklabels(klabels)

    ax1.set_ylim(erange)
    ax2.set_ylim(erange)
    ax2.yaxis.set_major_formatter(NullFormatter())
    ax1.set_ylabel(ylabel)
    ax2.set_xlabel(x2label)
    return ax1, ax2

def set_legend(plt, ncol=1, fs=7, loc="upper right", loc2=None,
        alpha=1.0, lw=0.2, length=1.0, labelspacing=0.3, borderpad=None):
    leg = plt.legend(
            loc=loc, ncol=ncol, fontsize=fs, fancybox=False, 
            facecolor="white", edgecolor="black", handletextpad=0.4,
            handlelength=length, labelspacing=labelspacing,
            borderpad=borderpad)
    if loc2 is not None:
        leg.set_bbox_to_anchor([loc2[0], loc2[1]])
    leg.get_frame().set_alpha(alpha)
    leg.get_frame().set_linewidth(lw)	
    return leg

def set4bandos():
    FIG_WIDTH = 3.3
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH*0.9))
    plt.subplots_adjust(
            left=0.14, bottom=0.14,
            right=0.98, top=0.98, wspace=0, hspace=0)
    return fig, plt

def set_axis_lim(ax, data, axis='x', alpha=0.05, scale='linear'):
    if scale == 'linear':
        vmin = np.min(data)
        vmax = np.max(data)
        x0 = vmin - alpha*(vmax - vmin)
        x1 = vmax + alpha*(vmax - vmin)
    elif scale == 'log':
        cmin = np.log10(np.min(data))
        cmax = np.log10(np.max(data))
        c0 = cmin - alpha*(cmax - cmin)
        c1 = cmax + alpha*(cmax - cmin)
        x0 = np.power(10, c0)
        x1 = np.power(10, c1)
    else:
        return None
    if axis == 'x':
        ax.set_xlim([x0, x1])
    else:
        ax.set_ylim([x0, x1])

