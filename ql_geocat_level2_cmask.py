import os, sys
from os import path, uname, mkdir
from glob import glob
import string, logging, traceback
from time import time

import numpy as np
from  numpy import ma as ma
import scipy as scipy

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import Colormap,normalize,LinearSegmentedColormap,ListedColormap,LogNorm

from matplotlib.figure import Figure

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl

from pyhdf.SD import SD

import viirs_edr_data


CMD = viirs_edr_data.CloudMaskData
cmap = ListedColormap(CMD.ViirsCMfillColours[0][1])

L2_file='./geocatL2.GOES-15.2014086.210000.hdf'
truth_L2_file='CSPP_GEO_DEMO_RESULTS/geocatL2.GOES-15.2014086.210000.hdf'

def get_L2(L2_file):
    L2_file_obj = SD(L2_file)
    cmask_obj = L2_file_obj.select('baseline_cmask_goes_nop_cloud_mask')
    cmask = cmask_obj.get()
    cmask_obj.endaccess()
    L2_file_obj.end()

    return cmask


def plot_cmask(cmask,plotTitle,pngFile):

    cmap = ListedColormap(CMD.ViirsCMfillColours[0][1])

    numCats = np.array(CMD.ViirsCMfillColours[0][1]).size
    numBounds = numCats + 1

    tickPos = np.arange(float(numBounds))/float(numCats)
    tickPos = tickPos[0 :-1] + tickPos[1]/2.

    # Create figure with default size, and create canvas to draw on
    scale=1.5
    fig = Figure(figsize=(scale*5,scale*5))
    canvas = FigureCanvas(fig)

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
    ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

    # Granule axis title
    ax_title = ppl.setp(ax,title=plotTitle)
    ppl.setp(ax_title,fontsize=12)
    ppl.setp(ax_title,family="sans-serif")

    vmin,vmax = 0,3

    im = ax.imshow(cmask,interpolation='nearest',vmin=vmin,vmax=vmax,cmap=cmap)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)
    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)

    # add a colorbar axis
    cax_rect = [0.05 , 0.05, 0.9 , 0.10 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    ppl.setp(cax.get_xticklabels(),fontsize=9)
    ppl.setp(cax.get_xticklines(),visible=False)

    # Set the colourbar tick locations and ticklabels
    #ppl.setp(cb.ax,xticks=CMD.ViirsCMTickPos) # In colorbar axis coords (0..1)
    tickpos_data_coords = vmax*tickPos
    cb.set_ticks(tickpos_data_coords) # In data coords (0..3)
    ppl.setp(cb.ax,xticklabels=CMD.ViirsCMtickNames[0][1])

    # Colourbar title
    cax_title = ppl.setp(cax,title='Cloud Mask')
    ppl.setp(cax_title,fontsize=10)

    # Redraw the figure
    canvas.draw()

    canvas.print_figure(pngFile,dpi=200)
    print "Writing to {}...".format(pngFile)


def plot_cmask_diff(cmask,truth_cmask,plotTitle,pngFile):

    # Create figure with default size, and create canvas to draw on
    scale=1.5
    fig = Figure(figsize=(scale*5,scale*5))
    canvas = FigureCanvas(fig)

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
    ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

    # Granule axis title
    ax_title = ppl.setp(ax,title=plotTitle)
    ppl.setp(ax_title,fontsize=12)
    ppl.setp(ax_title,family="sans-serif")

    vmin,vmax = -3.5,3.5
    numCats = 7

    cmap = cm.get_cmap('RdBu', numCats)

    numBounds = 2*numCats + 1

    tickPos = (np.linspace(0.,numBounds,numBounds)/numBounds)[1:-1]

    im = ax.imshow(truth_cmask-cmask,interpolation='nearest',vmin=vmin,vmax=vmax,cmap=cmap)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)
    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)

    # add a colorbar axis
    cax_rect = [0.05 , 0.05, 0.9 , 0.10 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    ppl.setp(cax.get_xticklabels(),fontsize=9)
    ppl.setp(cax.get_xticklines(),visible=True)

    # Set the colourbar tick locations and ticklabels
    tickpos_data_coords = tickPos * numCats + vmin
    cb.set_ticks(tickpos_data_coords[::2]) # In data coords (0..3)
    ticklabels = np.arange(numCats) + vmin + 0.5
    ppl.setp(cb.ax,xticklabels=ticklabels)

    # Colourbar title
    cax_title = ppl.setp(cax,title='Cloud Mask difference')
    ppl.setp(cax_title,fontsize=10)

    # Redraw the figure
    canvas.draw()

    canvas.print_figure(pngFile,dpi=200)
    print "Writing to {}...".format(pngFile)

