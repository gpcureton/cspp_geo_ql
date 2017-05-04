#!/usr/bin/env python
# encoding: utf-8
"""
ql_geocat_level2.py

Purpose: Plot a dataset from a geocat level-2 output file.

Preconditions:

Optional:

Minimum commandline:

    python ql_geocat_level2.py  INPUTFILE DATASET

where...

    INPUTFILE: The fully qualified path to the geocat level-2 input files.

    DATASET: Name of a dataset in the HDF4 or NetCDF file.


Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2015-03-05.
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

import os
import sys
import logging
import traceback
from os import path, uname, environ
import string
import re
import uuid
import argparse
from shutil import rmtree, copyfile
from glob import glob
from time import time
from datetime import datetime, timedelta

import numpy as np
from numpy import ma
import copy

from scipy import vectorize

from pyhdf.SD import SD
from netCDF4 import Dataset
from netCDF4 import num2date

import geocat_l2_data
from ql_geocat_common import Satellite_NetCDF
from ql_geocat_common import set_plot_navigation_bm as set_plot_navigation
from ql_geocat_common import set_plot_styles
from ql_geocat_common import list_l2_datasets as list_datasets

# every module should have a LOG object
LOG = logging.getLogger(__file__)


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    # Do we need the expert help messages...
    is_expert = False
    if '--expert' in sys.argv:
        expert_index = sys.argv.index('--expert')
        sys.argv[expert_index] = '--help'
        is_expert = True
    elif '-x' in sys.argv  :
        expert_index = sys.argv.index('-x')
        sys.argv[expert_index] = '--help'
        is_expert = True
    else:
        pass

    map_res_choice = ['c','l','i']

    defaults = {
                'input_file':None,
                'stride':1,
                'viewport_radius': 1.,
                'llcrnrx':None,
                'llcrnry':None,
                'urcrnrx':None,
                'urcrnry':None,
                'plotMin'  : None,
                'plotMax'  : None,
                'map_axis' : [0.10, 0.15, 0.80, 0.8],
                'cbar_axis' : [0.10 , 0.05, 0.8 , 0.05],
                'image_size' : [7.5, 7.5],
                'scatter_plot':False,
                'unnavigated':False,
                'list_datasets':False,
                'pointSize':1,
                'font_scale':1.,
                'map_res':'c',
                'cmap':None,
                'output_file':None,
                'outputFilePrefix' : None,
                'dpi':200
                }

    description = '''Create a plot of a level-2 dataset from a geocat netCDF4 file.'''

    usage = "usage: %prog [mandatory args] [options]"
    version = 'cspp-geo-geocat-1.0a3'

    parser = argparse.ArgumentParser(
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     add_help=False,
                                     )

    # Mandatory/positional arguments

    parser.add_argument(
                      action='store',
                      dest='input_file',
                      type=str,
                      help='''The fully qualified path to a single geocat level-2 NetCDF4''' \
                              ''' input file.'''
                      )

    parser.add_argument(
                      action="store",
                      dest="dataset",
                      type=str,
                      help='''The geocat level-2 dataset to plot. See the --list_datasets''' \
                              ''' option for available datasets.'''
                      )

    # Optional arguments

    parser.add_argument('--cbar-axis',
                      action="store",
                      dest="cbar_axis",
                      default=defaults["cbar_axis"],
                      type=float,
                      nargs=4,
                      metavar=('LEFT', 'BOTTOM', 'WIDTH', 'HEIGHT'),
                      help='''Set the colorbar axes within the figure at position''' \
                              ''' [*left*, *bottom*, *width*,\n*height*] where all quantities''' \
                              ''' are in the range [0..1]. [default: '{}']'''.format(
                                  defaults["cbar_axis"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--cbar-title',
                      action="store",
                      dest="cbar_title",
                      type=str,
                      help='''The colourbar title. Must be placed in double quotes.
                      ''' if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--cmap',
                      action="store",
                      dest="cmap",
                      default=defaults["cmap"],
                      type=str,
                      help='''The matplotlib colormap to use. See the --list_datasets option''' \
                              ''' for details and default values.
                              ''' if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('-d','--dpi',
                      action="store",
                      dest="dpi",
                      default='200.',
                      type=float,
                      help='''The resolution in dots per inch of the output png file.''' \
                              ''' [default: {}]'''.format(defaults["dpi"]
                                  ) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--font-scale',
                      action="store",
                      dest="font_scale",
                      default=defaults["font_scale"],
                      type=float,
                      help='''The scale factor to apply to the default font size for the plot''' \
                              ''' labels. [default: {}]'''.format(defaults["font_scale"]
                          ) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--image-size',
                      action="store",
                      dest="image_size",
                      default=defaults["image_size"],
                      type=float,
                      nargs=2,
                      metavar=('WIDTH', 'HEIGHT'),
                      help='''The size of the output image [*width*, *height*] in inches.''' \
                              ''' [default: '{}']'''.format(defaults["image_size"]
                          ) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--list-datasets',
                      action="store_true",
                      dest="list_datasets",
                      default=defaults["list_datasets"],
                      help='''List the available datasets, the default colormap and whether a''' \
                              ''' log plot is created\nby default, then exit. The required''' \
                              ''' dataset must be given as 'None'.''')

    parser.add_argument('--logscale',
                      action="store_true",
                      dest="logscale",
                      help='''Plot the dataset using a logarithmic scale.
                      ''' if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--map-axis',
                      action="store",
                      dest="map_axis",
                      default=defaults["map_axis"],
                      type=float,
                      nargs=4,
                      metavar=('LEFT', 'BOTTOM', 'WIDTH', 'HEIGHT'),
                      help='''Set the map axes within the figure at position [*left*, *bottom*,''' \
                              ''' *width*, *height*]\nwhere all quantities are in the range''' \
                              ''' [0..1]. [default: '{}']'''.format(defaults["map_axis"]
                                  ) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('-m','--map-res',
                      action="store",
                      dest="map_res",
                      default=defaults["map_res"],
                      type=str,
                      choices=map_res_choice,
                      help='''The map coastline resolution. Possible values are 'c' (coarse),''' \
                              ''' 'l' (low) and\n'i' (intermediate). [default: '{}']'''.format(
                                  defaults["map_res"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--no-logscale',
                      action="store_true",
                      dest="no_logscale",
                      help='''Plot the dataset using a linear scale.
                      ''' if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('-o','--output-file',
                      action="store",
                      dest="output_file",
                      default=defaults["output_file"],
                      type=str,
                      help='''The filename of the output png file.'''
                      )

    parser.add_argument('-O','--output-file-prefix',
                      action="store",
                      dest="outputFilePrefix",
                      default=defaults["outputFilePrefix"],
                      type=str,
                      help='''String to prepend to the automatically generated png names.''' \
                              ''' [default: {}]'''.format(defaults["outputFilePrefix"])
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

    parser.add_argument('--plot-title',
                      action="store",
                      dest="plot_title",
                      type=str,
                      help='''The plot title. Must be placed in double quotes.''' \
                              if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('-P','--pointSize',
                      action="store",
                      dest="pointSize",
                      default=defaults["pointSize"],
                      type=float,
                      help='''Size of the plot point used to represent each pixel.''' \
                              ''' [default: {}]'''.format(defaults["pointSize"]
                          ) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--scatter-plot',
                      action="store_true",
                      dest="doScatterPlot",
                      default=defaults["scatter_plot"],
                      help="Generate the plot using a scatterplot approach."
                      )

    parser.add_argument('-S','--stride',
                      action="store",
                      dest="stride",
                      default=defaults["stride"],
                      type=int,
                      help='''Sample every STRIDE rows and columns in the data.''' \
                              ''' [default: {}]'''.format(defaults["stride"])
                      )

    parser.add_argument('--unnavigated',
                      action="store_true",
                      dest="unnavigated",
                      default=defaults["unnavigated"],
                      help="Do not navigate the data, just display the image."
                      )

    parser.add_argument('--subset',
                      action="store",
                      dest="viewport",
                      type=float,
                      nargs=4,
                      metavar=('LLCRNRX', 'LLCRNRY', 'URCRNRX', 'URCRNRY'),
                      help='''Lower-left and upper-right coordinates [*llcrnrx*, *llcrnry*,''' \
                              ''' *urcrnrx*, *urcrnry*]\nof the projection viewport, where the''' \
                              ''' default is [-0.5,-0.5,+0.5,+0.5] for a full\ndisk (for''' \
                              ''' navigated plots only)''' if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--subset-lat0',
                      action="store",
                      dest="lat_0",
                      type=float,
                      help='''Center the plot over this latitude in the range [-90..90]''' \
                              ''' degrees. Must be used\nwith the option --subset-lon0.
                              ''' if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--subset-lon0',
                      action="store",
                      dest="lon_0",
                      type=float,
                      help='''Center the plot over this longitude in the range [-180..180]''' \
                              ''' degrees. Must be used\nwith the option --subset-lat0.
                              ''' if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--subset-radius',
                      action="store",
                      dest="viewport_radius",
                      type=float,
                      default=defaults["viewport_radius"],
                      help='''The radius in degrees of the subset region centered on the''' \
                              ''' coordinate selected\nby --subset-lat0 and --subset-lon.'''
                              ''' [default: {} degrees]'''.format(defaults["viewport_radius"]
                                  ) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument("-v", "--verbosity",
                      dest='verbosity',
                      action="count",
                      default=2,
                      help='''each occurrence increases verbosity 1 level from ERROR:''' \
                              ''' -v=WARNING -vv=INFO -vvv=DEBUG'''
                      )

    parser.add_argument('-V', '--version',
                      action='version',
                      version=version,
                      help='''Print the CSPP Geo package version'''
                      )

    parser.add_argument("-h", "--help",
                      action='help',
                      help='''Show this help message and exit.'''
                      )

    parser.add_argument('-x','--expert',
                      action="store_true",
                      dest="is_expert",
                      default=False,
                      help="Display all help options, including the expert ones."
                      )

    args = parser.parse_args()

    # Set up the logging
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    level = levels[args.verbosity if args.verbosity < 4 else 3]

    if level == logging.DEBUG :
        console_logFormat = '%(asctime)s.%(msecs)03d (%(levelname)s) : %(filename)s : %(funcName)s : %(lineno)d:%(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
    else:
        console_logFormat = '%(asctime)s.%(msecs)03d (%(levelname)s) : %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

    logging.basicConfig(stream = sys.stdout,
                        level = level,
                        format = console_logFormat,
                        datefmt = date_format)

    #
    # Enforce any mutual exclusivity or other relationships between various options.
    #

    # We can't give --subset and --subset-lat0 or --subset-lon0 simultaneously...
    if args.viewport is not None and (args.lat_0 is not None or args.lon_0 is not None):
        parser.error('''Cannot give --subset and (--subset-lat0/--subset-lon0) simultaneously.''')

    # Both subset  latitude and longitude must be used together...
    if (args.lat_0 is not None and args.lon_0 is None) or \
       (args.lat_0 is None and args.lon_0 is not None):
        parser.error('''--subset-lat0 and --subset-lon0 must be used together.''')

    return args,version


def main():
    '''
    The main method.
    '''

    # Read in the options
    options,cspp_geo_version = _argparse()

    output_file  = options.output_file
    outputFilePrefix  = options.outputFilePrefix

    # Create and populate the satellite object
    sat_l2_obj = Satellite_NetCDF(options.input_file)

    subsatellite_lon = sat_l2_obj.attrs['Subsatellite_Longitude']
    LOG.info("File subsatellite_lon = {}".format(subsatellite_lon))

    # If we want to list the datasets, do that here and exit
    if options.list_datasets:
        list_datasets(options,sat_l2_obj,geocat_l2_data)
        return 0

    lats = sat_l2_obj.Dataset(sat_l2_obj,'pixel_latitude').dset
    lons = sat_l2_obj.Dataset(sat_l2_obj,'pixel_longitude').dset
    sat_zenith_angle = sat_l2_obj.Dataset(sat_l2_obj,'pixel_satellite_zenith_angle').dset

    # Read in the desired dataset
    try:
        #chan_convention = sat_l2_obj.attrs['Channel_Number_Convention']
        #LOG.info('chan_convention: {}'.format(chan_convention))
        if 'Channel_Number_Convention' in sat_l2_obj.attrs.keys():
            LOG.warn('Channel_Number_Convention attribute in \n\t{}, is this a level-1 file? Aborting.\n'
                    .format(options.input_file))
            return 1

        data_obj = sat_l2_obj.Dataset(sat_l2_obj,options.dataset)
    except Exception:
        LOG.debug(traceback.format_exc())
        LOG.error('"{}" is not a valid options.dataset in {}, aborting.'.format(options.dataset,options.input_file))
        sat_l2_obj.close_file()
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

    sat_l2_obj.close_file()

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

    dataset = options.dataset

    if 'goes' in  sat_l2_obj.attrs['Sensor_Name']:
        sat_obj = geocat_l2_data.Satellite.factory('GOES_NOP')
    elif 'himawari' in  sat_l2_obj.attrs['Sensor_Name']:
        sat_obj = geocat_l2_data.Satellite.factory('Himawari')
    else:
        LOG.error("Unsupported satellite {}, aborting...".format(sat_l2_obj.attrs['Sensor_Name']))

    sat_obj.set_subsatellite_lon(subsatellite_lon)

    # Get the dataset options
    try:
        dataset_options = sat_obj.data[dataset]
    except KeyError:
        dataset_options = sat_obj.data['unknown']
        dataset_options['name'] = dataset

    # Set the navigation
    if options.unnavigated :
        plot_nav_options = {}
    else:
        plot_nav_options = set_plot_navigation(lats, lons, sat_l2_obj, options)

    # Set the plot styles
    plot_style_options = set_plot_styles(sat_l2_obj,data_obj,
            dataset_options,options,plot_nav_options)

    plot_style_options['version'] = cspp_geo_version

    # Get pointers to the desired plotting routines
    plot_image = plot_style_options['plot_image']
    plot_map = plot_style_options['plot_map']

    # Create the plot
    if options.unnavigated :

        plot_image(data, data_mask, output_file, dataset_options, plot_style_options)

    else :

        plot_map(lats, lons, data, data_mask, output_file,
                dataset_options, plot_nav_options, plot_style_options)

    return 0


if __name__=='__main__':
    sys.exit(main())
