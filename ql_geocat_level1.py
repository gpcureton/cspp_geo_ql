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
        cax_rect = [0.05 , 0.10, 0.9 , 0.05 ] # [left,bottom,width,height]
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

        canvas.print_figure(png_file,dpi=dpi)
        print "Writing to {}...".format(png_file)


    def plot_L1_Map(self,lat,lon,data,data_mask,pngName,**plot_options):
            
        # Copy the plot options to local variables
        title         = plot_options['title']
        cbar_title    = plot_options['cbar_title']
        units         = plot_options['units']
        stride        = plot_options['stride']
        lat_0         = plot_options['lat_0']
        lon_0         = plot_options['lon_0']
        latMin        = plot_options['latMin']
        lonMin        = plot_options['lonMin']
        latMax        = plot_options['latMax']
        lonMax        = plot_options['lonMax']
        plotMin       = plot_options['plotMin']
        plotMax       = plot_options['plotMax']
        map_res       = plot_options['map_res']
        cmap          = plot_options['cmap']
        doScatterPlot = plot_options['scatterPlot']
        pointSize     = plot_options['pointSize']
        dpi           = plot_options['dpi']

        '''
        Plot the input dataset in mapped to particular projection
        '''

        # If our data is all missing, return
        if (np.sum(data_mask) == data.size):
            LOG.warn("Entire {} dataset is missing, aborting".\
                    format(cbar_title))
            return -1

        # Compute the central lat and lon if they are not specified
        if (lat_0==None) and (lon_0==None):
            geo_shape = lat.shape
            nrows, ncols = geo_shape[0],geo_shape[1]
            LOG.debug("nrows,ncols= ({},{})".format(nrows,ncols))
            # Non lat/lon pair given, use the central unmasked values.
            row_idx = int(nrows/2.)
            col_idx = int(ncols/2.)
            lat_0 = lat[row_idx,col_idx]
            lon_0 = lon[row_idx,col_idx]
            LOG.info("No lat/lon pair given, using ({:4.2f},{:4.2f})".
                    format(lat_0,lon_0))

        if (latMin==None) and (latMax==None):
            LOG.info("Calculating lat extent...")
            latMin = np.min(lat)
            latMax = np.max(lat)
            LOG.info("Latitude extent: ({:4.2f},{:4.2f})".format(latMin,latMax))
        if (lonMin==None) and (lonMax==None):
            LOG.info("Calculating lon extent...")
            lonMin = np.min(lon)
            lonMax = np.max(lon)
            LOG.info("Longitude extent: ({:4.2f},{:4.2f})".format(lonMin,lonMax))

        # General Setup
        figWidth,figHeight = 5.,5.
        ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]

        fig = Figure(figsize=(figWidth,figHeight))
        canvas = FigureCanvas(fig)

        ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

        m = Basemap(projection='geos',lon_0=lon_0,ax=ax,fix_aspect=True,resolution=map_res)


        x,y=m(lon[::stride,::stride],lat[::stride,::stride])

        m.drawcoastlines(ax=ax,color='white')
        m.drawcountries(ax=ax,color='white')
        m.fillcontinents(color='0.85',zorder=0)
        m.drawparallels(np.arange( -90, 91,30), color = '0.25', 
                linewidth = 0.5)
        m.drawmeridians(np.arange(-180,180,30), color = '0.25', 
                linewidth = 0.5)

        data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

        plotMin = np.min(data) if plotMin==None else plotMin
        plotMax = np.max(data) if plotMax==None else plotMax
        LOG.debug("plotMin = {}".format(plotMin))
        LOG.debug("plotMax = {}".format(plotMax))

        if doScatterPlot:
            cs = m.scatter(x,y,s=pointSize,c=data,axes=ax,edgecolors='none',
                    vmin=plotMin,vmax=plotMax,cmap=cmap)
        else:
            cs = m.pcolor(x,y,data,axes=ax,edgecolors='none',antialiased=False,
                    vmin=plotMin,vmax=plotMax,cmap=cmap)

        txt = ax.set_title(title,fontsize=11)

        ppl.setp(ax.get_xticklines(),visible=False)
        ppl.setp(ax.get_yticklines(),visible=False)
        ppl.setp(ax.get_xticklabels(), visible=False)
        ppl.setp(ax.get_yticklabels(), visible=False)

        #ax.set_aspect('equal')

        cax_rect = [0.05 , 0.05, 0.9 , 0.05 ] # [left,bottom,width,height]
        cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes
        cb = fig.colorbar(cs, cax=cax, orientation='horizontal')

        txt = cax.set_title(cbar_title)

        #
        # Add a small globe with the swath indicated on it #
        #

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        glax_rect = [0.81, 0.75, 0.18, 0.20 ] # [left,bottom,width,height]
        glax = fig.add_axes(glax_rect)

        m_globe = Basemap(lat_0=0.,lon_0=0.,\
            ax=glax,resolution='c',area_thresh=10000.,projection='robin')

        # If we previously had a zero size data array, increase the pointSize
        # so the data points are visible on the global plot
        if (np.shape(lon[::stride,::stride])[0]==2) :
            pointSize = 5.

        x,y=m_globe(lon[::stride,::stride],lat[::stride,::stride])
        swath = np.zeros(np.shape(x),dtype=int)

        m_globe.drawmapboundary(linewidth=0.1)
        m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
        m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

        p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none',zorder=2)

        # Redraw the figure
        canvas.draw()
        LOG.info("Writing image file {}".format(pngName))
        canvas.print_figure(pngName,dpi=dpi)


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
                'lat_0':None,
                'lat_0':None,
                'plotMin'  : None,
                'plotMax'  : None,
                'scatter_plot':False,
                'pointSize':1,
                'map_res':'c',
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

    parser.add_argument('--lat_0',
                      action="store",
                      dest="lat_0",
                      type=float,
                      help="Center latitude of plot."
                      )

    parser.add_argument('--lon_0',
                      action="store",
                      dest="lon_0",
                      type=float,
                      help="Center longitude of plot."
                      )

    parser.add_argument('--latMin',
                      action="store",
                      dest="latMin",
                      type=float,
                      help="Minimum latitude to plot."
                      )

    parser.add_argument('--latMax',
                      action="store",
                      dest="latMax",
                      type=float,
                      help="Maximum latitude to plot."
                      )

    parser.add_argument('--lonMin',
                      action="store",
                      dest="lonMin",
                      type=float,
                      help="Minimum longitude to plot."
                      )

    parser.add_argument('--lonMax',
                      action="store",
                      dest="lonMax",
                      type=float,
                      help="Maximum longitude to plot."
                      )

    parser.add_argument('--scatter_plot',
                      action="store_true",
                      dest="doScatterPlot",
                      default=defaults["scatter_plot"],
                      help="Generate the plot using a scatterplot approach."
                      )

    parser.add_argument('-P','--pointSize',
                      action="store",
                      dest="pointSize",
                      default=defaults["pointSize"],
                      type=float,
                      help='''Size of the plot point used to represent each pixel. 
                      [default: {}]'''.format(defaults["pointSize"])
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
    lat_0  = options.lat_0
    lon_0  = options.lon_0
    latMin = options.latMin
    lonMin = options.lonMin
    latMax = options.latMax
    lonMax = options.lonMax
    lonMax = options.lonMax
    plotMin = options.plotMin
    plotMax = options.plotMax
    doScatterPlot = options.doScatterPlot
    pointSize = options.pointSize
    map_res = options.map_res
    output_file  = options.output_file
    outputFilePrefix  = options.outputFilePrefix
    dpi = options.dpi

    # Create and populate the GOES-L1 object

    goes_l1_obj = GOES_L1(input_file)

    lats = goes_l1_obj.Dataset(goes_l1_obj,'pixel_latitude').dset
    lons = goes_l1_obj.Dataset(goes_l1_obj,'pixel_longitude').dset
    data_obj = goes_l1_obj.Dataset(goes_l1_obj,dataset)

    LOG.info('Subsatellite_Longitude = {}'.format(goes_l1_obj.file_obj.attrs['units']))

    sys.exit(0)

    data = data_obj.dset
    data_mask = data.mask
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
    plot_options['stride'] = stride
    plot_options['lat_0'] = lat_0
    plot_options['lon_0'] = lon_0
    plot_options['latMin'] = latMin
    plot_options['lonMin'] = lonMin
    plot_options['latMax'] = latMax
    plot_options['lonMax'] = lonMax
    plot_options['plotMin'] = plotMin
    plot_options['plotMax'] = plotMax
    plot_options['map_res'] = map_res
    plot_options['cmap'] = cmap_dict[dataset]
    plot_options['scatterPlot'] = doScatterPlot
    plot_options['pointSize'] = pointSize
    plot_options['dpi'] = dpi

    # Create the plot
    #goes_l1_obj.plot_L1(data,output_file,**plot_options)
    #goes_l1_obj.plot_L1_Map(lats,lons,data,data_mask,output_file,**plot_options)

    return 0


if __name__=='__main__':
    sys.exit(main())  
