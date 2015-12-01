#!/usr/bin/env python
# encoding: utf-8
"""
ql_geocat_level2.py

Purpose: Plot a dataset from a geocat level-2 HDF4 file.

Preconditions:
    * pyhdf HDF4 python module

Optional:
    * 

Minimum commandline:

    python ql_geocat_level2.py  INPUTFILE DATASET

where...

    INPUTFILE: The fully qualified path to the geocat level-2 input files.

    DATASET: One of .


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

from pyhdf.SD import SD
from netCDF4 import Dataset
from netCDF4 import num2date

import geocat_l2_data
from ql_geocat_common import GOES_NetCDF
from ql_geocat_common import set_plot_navigation_bm as set_plot_navigation
from ql_geocat_common import set_plot_styles
from ql_geocat_common import list_l2_datasets as list_datasets

# every module should have a LOG object
LOG = logging.getLogger(__file__)


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    prodChoices=[
        'pixel_latitude',
        'pixel_longitude',
        'pixel_solar_zenith_angle',
        'pixel_satellite_zenith_angle',
        'pixel_relative_azimuth_angle',
        'pixel_surface_type',
        'pixel_ecosystem_type',
        'nwp_x_index',
        'nwp_y_index',
        'baseline_cmask_goes_nop_cloud_mask',
        'goesnp_ctype_cloud_type',
        'goesnp_ctype_cloud_phase',
        'ACHA_mode_6_cloud_top_temperature',
        'ACHA_mode_6_cloud_top_pressure',
        'ACHA_mode_6_cloud_top_height',
        'ACHA_mode_6_cloud_emissivity',
        'ACHA_mode_7_goes_cloud_optical_depth_vis',
        'ACHA_mode_7_goes_cloud_top_temperature',
        'ACHA_mode_7_goes_cloud_top_pressure',
        'ACHA_mode_7_goes_cloud_top_height',
        'ACHA_mode_7_goes_cloud_emissivity',
        'ACHA_mode_7_goes_cloud_particle_effective_radius',
        'DCOMP_mode_3_cloud_optical_depth_vis',
        'DCOMP_mode_3_cloud_particle_effective_radius',
        'DCOMP_mode_3_cloud_liquid_water_path',
        'DCOMP_mode_3_cloud_ice_water_path',
        'DCOMP_mode_3_cloud_albedo',
        'goesr_fog_fog_mask',
        'goesr_fog_MVFR_fog_probability',
        'goesr_fog_LIFR_fog_probability',
        'goesr_fog_IFR_fog_probability',
        'goesr_fog_IFR_RHonly_Fog_Probability',
        'goesr_fog_fog_depth',
        'goesr_fog_ems7_atmospherically_corrected',
        'goesr_fog_surface_temperature_bias',
        'goesr_fog_Surface_Temperature_Bias_Global',
        'goesr_fog_ref2_stddev',
        'goesr_fog_Sfc_Emiss_Chn7',
        'goesr_fog_Ems7_Composite',
        'goesr_fog_Ref7_Composite',
        'goesr_fog_Ref2_Stddev_Composite',
        'goesr_fog_bt14_stddev',
        'goesr_fog_Max_RH_500ft_Layer_AGL',
        'goesr_fog_Max_RH_1000ft_Layer_AGL',
        'goesr_fog_Max_RH_3000ft_Layer_AGL',
        'goesr_fog_Surface_RH',
        'goesr_fog_Refl_Chn2_StdDev_Lrc',
        ]

    map_res_choice = ['c','l','i']

    #goes_choice = ['goes_e','goes_w']
    #goes_region_choice = ['FD','CONUS','MESO']
    goes_region_choice = ['FD']

    defaults = {
                'input_file':None,
                'dataset':'baseline_cmask_goes_nop_cloud_mask',
                'stride':1,
                #'lon_0':None,
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
                'font_scale':1.,
                'map_res':'c',
                'cmap':None,
                'output_file':None,
                'outputFilePrefix' : None,
                'dpi':200
                }

    description = '''Create a plot of a level-2 dataset from a geocat netCDF4 file.'''

    usage = "usage: %prog [mandatory args] [options]"
    version = 'CSPP Geo GEOCAT v1.0beta'

    parser = argparse.ArgumentParser(
                                     description=description
                                     )

    # Mandatory/positional arguments
    
    parser.add_argument(
                      action='store',
                      dest='input_file',
                      type=str,
                      help='''The fully qualified path to a single geocat level-2 
                      NetCDF4 input file.'''
                      )

    parser.add_argument(
                      action="store",
                      dest="dataset",
                      type=str,
                      help='''The geocat level-2 dataset to plot. See the 
                      --list_datasets option for available datasets.'''
                      )

    # Optional arguments 

    parser.add_argument('--cbar_axis',
                      action="store",
                      dest="cbar_axis",
                      default=defaults["cbar_axis"],
                      type=float,
                      nargs=4,
                      metavar=('LEFT', 'BOTTOM', 'WIDTH', 'HEIGHT'),
                      help="""Set the colorbar axes within the figure at position 
                      [*left*, *bottom*, *width*, *height*] where all quantities 
                      are in the range [0..1].[default: '{}']
                      """.format(defaults["cbar_axis"])
                      )

    parser.add_argument('--cbar_title',
                      action="store",
                      dest="cbar_title",
                      type=str,
                      help='''The colourbar title. Must be placed in double quotes.
                      '''
                      )

    parser.add_argument('--cmap',
                      action="store",
                      dest="cmap",
                      default=defaults["cmap"],
                      type=str,
                      help="""The matplotlib colormap to use. See the --list_datasets
                      option for details and default values."""
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

    parser.add_argument('--font_scale',
                      action="store",
                      dest="font_scale",
                      default=defaults["font_scale"],
                      type=float,
                      help='''The scale factor to apply to the default font size
                      for the plot labels. [default: {}]'''.format(defaults["font_scale"])
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

    parser.add_argument('--list_datasets',
                      action="store_true",
                      dest="list_datasets",
                      default=defaults["list_datasets"],
                      help="""List the available datasets, the default colormap
                      and whether a log plot is created by default, then exit. 
                      The required dataset may be given as 'None'."""
                      )

    parser.add_argument('--logscale',
                      action="store_true",
                      dest="logscale",
                      help="""Plot the dataset using a logarithmic scale."""
                      )

    #parser.add_argument('--lon_0',
                      #action="store",
                      #dest="lon_0",
                      #type=float,
                      #help="Center the plot over this longitude."
                      #)

    parser.add_argument('--map_axis',
                      action="store",
                      dest="map_axis",
                      default=defaults["map_axis"],
                      type=float,
                      nargs=4,
                      metavar=('LEFT', 'BOTTOM', 'WIDTH', 'HEIGHT'),
                      help="""Set the map axes within the figure at position 
                      [*left*, *bottom*, *width*, *height*] where all quantities 
                      are in the range [0..1].[default: '{}']
                      """.format(defaults["map_axis"])
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

    parser.add_argument('--no_logscale',
                      action="store_true",
                      dest="no_logscale",
                      help="""Plot the dataset using a linear scale."""
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

    parser.add_argument('--plot_title',
                      action="store",
                      dest="plot_title",
                      type=str,
                      help='''The plot title. Must be placed in double quotes.
                      '''
                      )

    parser.add_argument('-P','--pointSize',
                      action="store",
                      dest="pointSize",
                      default=defaults["pointSize"],
                      type=float,
                      help='''Size of the plot point used to represent each pixel. 
                      [default: {}]'''.format(defaults["pointSize"])
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

    #parser.add_argument('--satellite',
                      #action="store",
                      #dest="satellite",
                      #type=str,
                      #choices=goes_choice,
                      #help="""The GOES satellite."""
                      #)

    parser.add_argument('--scatter_plot',
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
                      help='''Sample every STRIDE rows and columns in the data. 
                      [default: {}]'''.format(defaults["stride"])
                      )

    parser.add_argument('--unnavigated',
                      action="store_true",
                      dest="unnavigated",
                      default=defaults["unnavigated"],
                      help="Do not navigate the data, just display the image."
                      )

    parser.add_argument('--viewport',
                      action="store",
                      dest="viewport",
                      type=float,
                      nargs=4,
                      metavar=('LLCRNRX', 'LLCRNRY', 'URCRNRX', 'URCRNRY'),
                      help="""Lower-left and upper-right coordinates 
                      [*llcrnrx*, *llcrnry*, *urcrnrx*, *urcrnry*] of the projection 
                      viewport, where the default is [-0.5,-0.5,+0.5,+0.5] for a 
                      full disk (for navigated plots only)"""
                      )

    parser.add_argument("-v", "--verbosity",
                      dest='verbosity',
                      action="count", 
                      default=2,
                      help='''each occurrence increases verbosity 1 level from 
                      ERROR: -v=WARNING -vv=INFO -vvv=DEBUG'''
                      )

    parser.add_argument('-V', '--version', 
                      action='version',
                      version=version,
                      help='''Print the CSPP Geo package version'''
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

    logging.basicConfig(level=level, format=console_logFormat, datefmt=date_format)

    return args,version


def main():
    '''
    The main method.
    '''

    # Read in the options
    options,cspp_geo_version = _argparse()

    output_file  = options.output_file
    outputFilePrefix  = options.outputFilePrefix

    # Create and populate the GOES object
    #goes_l2_obj = GOES_HDF4(options.input_file)
    goes_l2_obj = GOES_NetCDF(options.input_file)

    lats = goes_l2_obj.Dataset(goes_l2_obj,'pixel_latitude').dset
    lons = goes_l2_obj.Dataset(goes_l2_obj,'pixel_longitude').dset
    sat_zenith_angle = goes_l2_obj.Dataset(goes_l2_obj,'pixel_satellite_zenith_angle').dset

    # If we want to list the datasets, do that here and exit
    if options.list_datasets:
        list_datasets(options,goes_l2_obj,geocat_l2_data)
        return 0

    # Read in the desired dataset
    try:

        LOG.debug('options.dataset name: {}'.format(options.dataset))
        if 'Channel_Number_Convention' in goes_l2_obj.attrs.keys():
            LOG.warn('Channel_Number_Convention attribute in \n\t{}, is this a level-1 file? Aborting.\n'
                    .format(options.input_file))
            return 1
        
        data_obj = goes_l2_obj.Dataset(goes_l2_obj,options.dataset)
    except Exception:
        LOG.debug(traceback.format_exc())
        LOG.error('"{}" is not a valid options.dataset in {}, aborting.'.format(options.dataset,options.input_file))
        goes_l2_obj.close_netcdf_file()
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

    goes_l2_obj.close_netcdf_file()

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

    # Get the dataset options
    try:
        dataset_options = geocat_l2_data.Dataset_Options.data[dataset]
    except KeyError:
        dataset_options = geocat_l2_data.Dataset_Options.data['unknown']
        dataset_options['name'] = dataset

    # Set the navigation 
    if options.unnavigated :
        plot_nav_options = {}
    else:
        plot_nav_options = set_plot_navigation(lats,lons,goes_l2_obj,options)

    # Set the plot styles 
    plot_style_options = set_plot_styles(goes_l2_obj,data_obj,
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
