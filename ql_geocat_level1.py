#!/usr/bin/env python
# encoding: utf-8
"""
ql_geocat_level1.py

Purpose: Plot a dataset from a geocat level-1 HDF4 file.

Preconditions:
    * pyhdf HDF4 python module

Optional:
    * 

Minimum commandline:

    python ql_geocat_level1.py  INPUTFILE DATASET

where...

    INPUTFILE: The fully qualified path to the geocat level-1 input files.

    DATASET: One of .


Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2015-03-04.
Copyright (c) 2015 University of Wisconsin Regents. All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os, sys, logging, traceback
from os import path,uname,environ
import string
import re
import uuid
from shutil import rmtree,copyfile
from glob import glob
from time import time
from datetime import datetime,timedelta

import numpy as np
from numpy import ma
import copy

from scipy import vectorize

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl

from mpl_toolkits.basemap import Basemap

from pyhdf.SD import SD

# every module should have a LOG object
LOG = logging.getLogger(__file__)


class GOES_L1():

    def __init__(self,l1_file):

        self.l1_file = l1_file
        self.file_obj = SD(l1_file)
        self.data_dict = self.file_obj.datasets()
        self.datanames = self.data_dict.keys()
        self.datanames.sort()


    class Dataset():

        def __init__(selfd,L2_obj,dataname):

            selfd.dataname = dataname
            selfd.dset_obj = L2_obj.file_obj.select(dataname)
            selfd.attrs = selfd.dset_obj.attributes()
            selfd.dset = ma.masked_equal(selfd.dset_obj.get(),selfd.attrs['_FillValue'])
            selfd.dset = selfd.dset * selfd.attrs['scale_factor'] + selfd.attrs['add_offset']


    def close(self):
        self.file_obj.end()


    def plot_L1(self,data,png_file,**plot_options):

        # Copy the plot options to local variables
        title         = plot_options['title']
        cbar_title    = plot_options['cbar_title']
        units         = plot_options['units']
        plotMin       = plot_options['plotMin']
        plotMax       = plot_options['plotMax']
        cmap          = plot_options['cmap']
        dpi           = plot_options['dpi']

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*5,scale*5))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

        # Granule axis title
        ax_title = ppl.setp(ax,title=title)
        ppl.setp(ax_title,fontsize=12)
        ppl.setp(ax_title,family="sans-serif")

        vmin,vmax = None,None

        im = ax.imshow(data,interpolation='nearest',vmin=vmin,vmax=vmax,cmap=cmap)
        ppl.setp(ax.get_xticklabels(), visible=False)
        ppl.setp(ax.get_yticklabels(), visible=False)
        ppl.setp(ax.get_xticklines(),visible=False)
        ppl.setp(ax.get_yticklines(),visible=False)

        # add a colorbar axis
        cax_rect = [0.05 , 0.05, 0.9 , 0.08 ] # [left,bottom,width,height]
        cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

        # Plot the colorbar.
        cb = fig.colorbar(im, cax=cax, orientation='horizontal')
        ppl.setp(cax.get_xticklabels(),fontsize=9)
        ppl.setp(cax.get_xticklines(),visible=False)

        # Colourbar title
        cax_title = ppl.setp(cax,title=cbar_title)
        ppl.setp(cax_title,fontsize=10)

        # Redraw the figure
        canvas.draw()

        canvas.print_figure(png_file,dpi=200)
        print "Writing to {}...".format(png_file)


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    prodChoices=[
                 'channel_14_brightness_temperature',
                 'channel_16_brightness_temperature',
                 'channel_2_reflectance',
                 'channel_7_brightness_temperature',
                 'channel_7_emissivity',
                 'channel_7_reflectance',
                 'channel_9_brightness_temperature',
                 'pixel_ecosystem_type',
                 'pixel_latitude',
                 'pixel_longitude',
                 'pixel_relative_azimuth_angle',
                 'pixel_satellite_zenith_angle',
                 'pixel_solar_zenith_angle',
                 'pixel_surface_type'
                ]

    map_res_choice = ['c','l','i']

    defaults = {
                'input_file':None,
                'dataset':'channel_2_reflectance',
                'stride':1,
                'plotMin'  : None,
                'plotMax'  : None,
                'output_file':None,
                'outputFilePrefix' : None,
                'dpi':200
                }

    description = '''Create a plot of a level-1 dataset from a geocat HDF4 file.'''

    usage = "usage: %prog [mandatory args] [options]"
    version = 0.1

    parser = argparse.ArgumentParser(
                                     description=description
                                     )

    # Mandatory/positional arguments
    
    parser.add_argument(
                      action='store',
                      dest='input_file',
                      type=str,
                      help='''The fully qualified path to a single geocat level-1 HDF4 input file.'''
                      )

    parser.add_argument(
                      action="store",
                      dest="dataset",
                      default=defaults["dataset"],
                      type=str,
                      choices=prodChoices,
                      help='''The geocat level-1 dataset to plot.
                              Possible values are...
                              {}.
                              [default: {}]
                           '''.format(prodChoices.__str__()[1:-1],
                               defaults["dataset"])
                      )

    # Optional arguments 

    parser.add_argument('-S','--stride',
                      action="store",
                      dest="stride",
                      default=defaults["stride"],
                      type=int,
                      help='''Sample every STRIDE rows and columns in the data. 
                      [default: {}]'''.format(defaults["stride"])
                      )

    parser.add_argument('--plotMin',
                      action="store",
                      dest="plotMin",
                      default=defaults["plotMin"],
                      type=float,
                      help="Minimum value to plot.".format(defaults["plotMin"])
                      )

    parser.add_argument('--plotMax',
                      action="store",
                      dest="plotMax",
                      default=defaults["plotMax"],
                      type=float,
                      help="Maximum value to plot.".format(defaults["plotMax"])
                      )

    parser.add_argument('-d','--dpi',
                      action="store",
                      dest="dpi",
                      default='200.',
                      type=float,
                      help='''The resolution in dots per inch of the output 
                      png file. 
                      [default: {}]'''.format(defaults["dpi"])
                      )

    parser.add_argument('-o','--output_file',
                      action="store",
                      dest="output_file",
                      default=defaults["output_file"],
                      type=str,
                      help='''The filename of the output png file. 
                      '''
                      )

    parser.add_argument('-O','--output_file_prefix',
                      action="store",
                      dest="outputFilePrefix",
                      default=defaults["outputFilePrefix"],
                      type=str,
                      help="""String to prepend to the automatically generated 
                      png names. [default: {}]""".format(defaults["outputFilePrefix"])
                      )

    parser.add_argument("-v", "--verbose",
                      dest='verbosity',
                      action="count", 
                      default=0,
                      help='''each occurrence increases verbosity 1 level from 
                      ERROR: -v=WARNING -vv=INFO -vvv=DEBUG'''
                      )


    args = parser.parse_args()

    # Set up the logging
    console_logFormat = '%(asctime)s : (%(levelname)s):%(filename)s:%(funcName)s:%(lineno)d:  %(message)s'
    dateFormat = '%Y-%m-%d %H:%M:%S'
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[args.verbosity], 
            format = console_logFormat, 
            datefmt = dateFormat)

    return args


def main():
    '''
    The main method.
    '''

    # Read in the options
    options = _argparse()

    input_file = options.input_file
    dataset = options.dataset
    stride = options.stride
    plotMin = options.plotMin
    plotMax = options.plotMax
    output_file  = options.output_file
    outputFilePrefix  = options.outputFilePrefix
    dpi = options.dpi

    # Create and populate the GOES-L1 object

    goes_l1_obj = GOES_L1(input_file)

    lats = goes_l1_obj.Dataset(goes_l1_obj,'pixel_latitude').dset
    lons = goes_l1_obj.Dataset(goes_l1_obj,'pixel_longitude').dset
    data_obj = goes_l1_obj.Dataset(goes_l1_obj,dataset)


    data = data_obj.dset
    plot_title = "{}".format(input_file)
    cbar_title = "{} ({})".format(data_obj.dataname,data_obj.attrs['units'])

    input_file = path.basename(input_file)

    goes_l1_obj.close()
    

    # Determine the filename
    file_suffix = "{}".format(dataset)

    if output_file==None and outputFilePrefix==None :
        output_file = "{}.{}.png".format(input_file,file_suffix)
    if output_file!=None and outputFilePrefix==None :
        pass
    if output_file==None and outputFilePrefix!=None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)
    if output_file!=None and outputFilePrefix!=None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)

    # Define a colormap for each dataset
    cmap_dict={
                 'channel_14_brightness_temperature': cm.Spectral_r,
                 'channel_16_brightness_temperature': cm.Spectral_r,
                 'channel_2_reflectance': cm.gray,
                 'channel_7_brightness_temperature': cm.Spectral_r,
                 'channel_7_emissivity': cm.Spectral_r,
                 'channel_7_reflectance': cm.gray,
                 'channel_9_brightness_temperature': cm.Spectral_r,
                 'pixel_ecosystem_type': cm.Spectral_r,
                 'pixel_latitude': cm.Spectral_r,
                 'pixel_longitude': cm.Spectral_r,
                 'pixel_relative_azimuth_angle': cm.Spectral_r,
                 'pixel_satellite_zenith_angle': cm.Spectral_r,
                 'pixel_solar_zenith_angle': cm.Spectral_r,
                 'pixel_surface_type': cm.Spectral_r

    }
    # Default plot labels
    plot_options = {}
    plot_options['title'] = plot_title
    plot_options['cbar_title'] = cbar_title
    plot_options['units'] = cbar_title
    plot_options['plotMin'] = plotMin
    plot_options['plotMax'] = plotMax
    plot_options['cmap'] = cmap_dict[dataset]
    plot_options['dpi'] = dpi

    # Create the plot
    goes_l1_obj.plot_L1(data,output_file,**plot_options)

    return 0


if __name__=='__main__':
    sys.exit(main())  
