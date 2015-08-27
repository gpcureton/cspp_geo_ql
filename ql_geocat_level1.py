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
from netCDF4 import Dataset
from netCDF4 import num2date

import geocat_l1_data
from ql_geocat_common import GOES_NetCDF
from ql_geocat_common import set_plot_navigation_bm as set_plot_navigation
from ql_geocat_common import set_plot_styles

# every module should have a LOG object
LOG = logging.getLogger(__file__)


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

    goes_choice = ['goes_e','goes_w']
    goes_region_choice = ['FD','CONUS','MESO']

    defaults = {
                'input_file':None,
                'dataset':'channel_2_reflectance',
                'stride':1,
                'lon_0':None,
                'llcrnrx':None,
                'llcrnry':None,
                'urcrnrx':None,
                'urcrnry':None,
                'plotMin'  : None,
                'plotMax'  : None,
                'region' : None,
                'map_axis' : [0.10, 0.15, 0.80, 0.8],
                'cbar_axis' : [0.10 , 0.05, 0.8 , 0.05],
                'image_size' : [7.5, 7.5],
                'scatter_plot':False,
                'unnavigated':False,
                'list_datasets':False,
                'pointSize':1,
                'map_res':'c',
                'cmap':None,
                'output_file':None,
                'outputFilePrefix' : None,
                'dpi':200
                }

    description = '''Create a plot of a level-1 dataset from a geocat netCDF4 file.'''

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
                      help='''The geocat level-1 dataset to plot.
                              [default: {}]
                           '''.format(defaults["dataset"])
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

    parser.add_argument('--lon_0',
                      action="store",
                      dest="lon_0",
                      type=float,
                      help="Center the plot over this longitude."
                      )

    parser.add_argument('--viewport',
                      action="store",
                      dest="viewport",
                      type=float,
                      nargs=4,
                      metavar=('LLCRNRX', 'LLCRNRY', 'URCRNRX', 'URCRNRY'),
                      help="""Lower-left and upper-right coordinates 
                      [*llcrnrx*, *llcrnry*, *urcrnrx*, *urcrnry*] of the projection 
                      viewport, in the range [-0.5,+0.5] (for navigated plots only)"""
                      )

    parser.add_argument('--image_size',
                      action="store",
                      dest="image_size",
                      default=defaults["image_size"],
                      type=float,
                      nargs=2,
                      metavar=('WIDTH', 'HEIGHT'),
                      help="""The size of the output image [*width*, *height*]
                      in inches. [default: '{}']""".format(defaults["image_size"])
                      )

    parser.add_argument('-m','--map_res',
                      action="store",
                      dest="map_res",
                      default=defaults["map_res"],
                      type=str,
                      choices=map_res_choice,
                      help="""The map coastline resolution. Possible values are 
                      'c' (coarse),'l' (low) and 'i' (intermediate). 
                      [default: '{}']""".format(defaults["map_res"])
                      )

    parser.add_argument('--region',
                      action="store",
                      dest="region",
                      default=defaults["region"],
                      type=str,
                      choices=goes_region_choice,
                      help="""The GOES region. 
                      [default: '{}']""".format(defaults["region"])
                      )

    parser.add_argument('--map_axis',
                      action="store",
                      dest="map_axis",
                      default=defaults["map_axis"],
                      type=float,
                      nargs=4,
                      metavar=('LEFT', 'BOTTOM', 'WIDTH', 'HEIGHT'),
                      help="""Set the map axes at position [*left*, *bottom*, *width*, *height*] 
                      where all quantities are in fractions of figure width and height. 
                      [default: '{}']""".format(defaults["map_axis"])
                      )

    parser.add_argument('--cbar_axis',
                      action="store",
                      dest="cbar_axis",
                      default=defaults["cbar_axis"],
                      type=float,
                      nargs=4,
                      metavar=('LEFT', 'BOTTOM', 'WIDTH', 'HEIGHT'),
                      help="""Set the colorbar axes at position [*left*, *bottom*, *width*, *height*] 
                      where all quantities are in fractions of figure width and height. 
                      [default: '{}']""".format(defaults["cbar_axis"])
                      )

    parser.add_argument('--satellite',
                      action="store",
                      dest="satellite",
                      type=str,
                      choices=goes_choice,
                      help="""The GOES satellite."""
                      )

    parser.add_argument('--scatter_plot',
                      action="store_true",
                      dest="doScatterPlot",
                      default=defaults["scatter_plot"],
                      help="Generate the plot using a scatterplot approach."
                      )

    parser.add_argument('--unnavigated',
                      action="store_true",
                      dest="unnavigated",
                      default=defaults["unnavigated"],
                      help="Do not navigate the data, just display the image."
                      )

    parser.add_argument('--list_datasets',
                      action="store_true",
                      dest="list_datasets",
                      default=defaults["list_datasets"],
                      help="""List the available datasets, and exit. Specify the
                      required dataset as 'None'."""
                      )

    parser.add_argument('-P','--pointSize',
                      action="store",
                      dest="pointSize",
                      default=defaults["pointSize"],
                      type=float,
                      help='''Size of the plot point used to represent each pixel. 
                      [default: {}]'''.format(defaults["pointSize"])
                      )

    parser.add_argument('--cmap',
                      action="store",
                      dest="cmap",
                      default=defaults["cmap"],
                      type=str,
                      help="""The matplotlib colormap to use. 
                      [default: '{}']""".format(defaults["cmap"])
                      )

    parser.add_argument('--plot_title',
                      action="store",
                      dest="plot_title",
                      type=str,
                      help='''The plot title. 
                      '''
                      )

    parser.add_argument('--cbar_title',
                      action="store",
                      dest="cbar_title",
                      type=str,
                      help='''The colourbar title. 
                      '''
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

    output_file  = options.output_file
    outputFilePrefix  = options.outputFilePrefix

    # Create and populate the GOES object
    #goes_l1_obj = GOES_HDF4(options.input_file)
    goes_l1_obj = GOES_NetCDF(options.input_file)

    lats = goes_l1_obj.Dataset(goes_l1_obj,'pixel_latitude').dset
    lons = goes_l1_obj.Dataset(goes_l1_obj,'pixel_longitude').dset
    sat_zenith_angle = goes_l1_obj.Dataset(goes_l1_obj,'pixel_satellite_zenith_angle').dset

    # If we want to list the datasets, do that here and exit
    if options.list_datasets:
        LOG.info('Datasets in {}:'.format(options.input_file))
        for dsets in goes_l1_obj.datanames:
            print "\t{}".format(dsets)
        goes_l1_obj.close()
        sys.exit(1)

    # Read in the desired dataset
    try:

        LOG.debug('options.dataset name: {}'.format(options.dataset))

        data_obj = goes_l1_obj.Dataset(goes_l1_obj,options.dataset)
    except :
        LOG.error('"{}" is not a valid options.dataset in {}, aborting.'.format(options.dataset,options.input_file))
        goes_l1_obj.close()
        return 1

    # Use the solar zenith angle to mask off-disk pixels...
    data = ma.masked_array(data_obj.dset,mask=sat_zenith_angle.mask)

    if ma.is_masked(data):
        if data.mask.shape == ():
            data_mask = np.ones(data.shape,dtype='bool')
        else:
            data_mask = data.mask
    else: 
        data_mask = np.zeros(data.shape,dtype='bool')

    goes_l1_obj.close()
    
    # Determine the filename
    file_suffix = "{}".format(options.dataset)

    if output_file==None and outputFilePrefix==None :
        output_file = "{}.{}.png".format(path.basename(options.input_file),file_suffix)
    if output_file!=None and outputFilePrefix==None :
        pass
    if output_file==None and outputFilePrefix!=None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)
    if output_file!=None and outputFilePrefix!=None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)

    # Determine the correct key for the channel name...
    chan_convention = goes_l1_obj.attrs['Channel_Number_Convention']
    if 'instrument-native' in chan_convention:
        spacecraft = goes_l1_obj.attrs['Spacecraft_Name']
        dataset_prefix = "{}_".format(string.replace(spacecraft.lower(),'-','_'))
        dataset = string.replace(options.dataset,dataset_prefix,"")

    # Get the dataset options
    dataset_options = geocat_l1_data.Dataset_Options.data[dataset]

    # Set the navigation 
    plot_nav_options = set_plot_navigation(lats,lons,goes_l1_obj,options)

    # Set the plot styles 
    plot_style_options = set_plot_styles(goes_l1_obj,data_obj,
            dataset_options,options,plot_nav_options)

    # Get pointers to the desired plotting routines
    plot_image = plot_style_options['plot_image']
    plot_map = plot_style_options['plot_map']

    # Create the plot
    if options.unnavigated :

        plot_image(data, data_mask, output_file, dataset_options, plot_style_options)

    else :

        plot_map(lats,lons, data, data_mask, output_file,
                dataset_options, plot_nav_options, plot_style_options)

    return 0


if __name__=='__main__':
    sys.exit(main())  
