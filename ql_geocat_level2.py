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

import geocat_l2_data

# every module should have a LOG object
LOG = logging.getLogger(__file__)


class GOES_L2():

    def __init__(self,l2_file):

        self.l2_file = l2_file
        self.file_obj = SD(l2_file)
        self.attrs = self.file_obj.attributes()
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
        LOG.info('Closing the file {}'.format(self.l2_file))
        self.file_obj.end()


    def plot_L2(self,data,data_mask,png_file,dataset,**plot_options):

        l2_plot_options = geocat_l2_data.Plot_Options.data[dataset]

        # Copy the plot options to local variables
        title         = plot_options['title']
        cbar_title    = plot_options['cbar_title']
        units         = plot_options['units']
        stride        = plot_options['stride']
        plotMin       = plot_options['plotMin']
        plotMax       = plot_options['plotMax']
        plotLims      = plot_options['plotLims']
        cmap          = plot_options['cmap']
        dpi           = plot_options['dpi']

        '''
        Plot the input dataset in in native data coordinates
        '''

        # If our data is all missing, return
        if (np.sum(data_mask) == data.size):
            LOG.warn("Entire {} dataset is missing, aborting".\
                    format(cbar_title))
            return -1

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*5,scale*5))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.10, 0.15, 0.80, 0.8  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

        # Granule axis title
        ax_title = ppl.setp(ax,title=title)
        ppl.setp(ax_title,fontsize=12)
        ppl.setp(ax_title,family="sans-serif")

        LOG.debug('data.shape = {}'.format(data.shape))
        LOG.debug('data_mask.shape = {}'.format(data_mask.shape))
        data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

        vmin,vmax = plotLims[0],plotLims[-1]

        im = ax.imshow(data,interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax,cmap=cmap)

        #txt = ax.set_title(title,fontsize=11)

        ppl.setp(ax.get_xticklabels(), visible=False)
        ppl.setp(ax.get_yticklabels(), visible=False)
        ppl.setp(ax.get_xticklines(),visible=False)
        ppl.setp(ax.get_yticklines(),visible=False)

        # add a colorbar axis
        cax_rect = [0.10 , 0.05, 0.8 , 0.05 ] # [left,bottom,width,height]
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

        # Write the figure to file
        canvas.print_figure(png_file,dpi=dpi)
        LOG.info("Writing to {}...".format(png_file))


    def plot_L2_discrete(self,data,data_mask,png_file,dataset,**plot_options):

        l2_plot_options = geocat_l2_data.Plot_Options.data[dataset]

        # Copy the plot options to local variables
        title         = plot_options['title']
        cbar_title    = plot_options['cbar_title']
        units         = plot_options['units']
        stride        = plot_options['stride']
        plotMin       = plot_options['plotMin']
        plotMax       = plot_options['plotMax']
        plotLims      = plot_options['plotLims']
        cmap          = plot_options['cmap']
        dpi           = plot_options['dpi']

        '''
        Plot the input dataset in in native data coordinates
        '''

        # If our data is all missing, return
        if (np.sum(data_mask) == data.size):
            LOG.warn("Entire {} dataset is missing, aborting".\
                    format(dataset))
            return -1

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*5,scale*5))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.10, 0.15, 0.80, 0.8  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

        # Granule axis title
        ax_title = ppl.setp(ax,title=title)
        ppl.setp(ax_title,fontsize=12)
        ppl.setp(ax_title,family="sans-serif")

        fill_colours = l2_plot_options['fill_colours']
        cmap = ListedColormap(fill_colours)

        numCats = np.array(fill_colours).size
        numBounds = numCats + 1

        tickPos = np.arange(float(numBounds))/float(numCats)
        tickPos = tickPos[0 :-1] + tickPos[1]/2.

        LOG.debug('data.shape = {}'.format(data.shape))
        LOG.debug('data_mask.shape = {}'.format(data_mask.shape))
        data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

        data_values = l2_plot_options['values']
        vmin,vmax = data_values[0],data_values[-1]

        im = ax.imshow(data,interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax,cmap=cmap)

        #txt = ax.set_title(title,fontsize=11)

        ppl.setp(ax.get_xticklabels(), visible=False)
        ppl.setp(ax.get_yticklabels(), visible=False)
        ppl.setp(ax.get_xticklines(),visible=False)
        ppl.setp(ax.get_yticklines(),visible=False)

        # add a colorbar axis
        cax_rect = [0.10 , 0.05, 0.8 , 0.05 ] # [left,bottom,width,height]
        cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

        # Plot the colorbar.
        cb = fig.colorbar(im, cax=cax, orientation='horizontal')
        ppl.setp(cax.get_xticklabels(),fontsize=9)
        ppl.setp(cax.get_xticklines(),visible=False)

        # Set the colourbar tick locations and ticklabels
        #ppl.setp(cb.ax,xticks=CMD.ViirsCMTickPos) # In colorbar axis coords (0..1)
        tickpos_data_coords = vmax*tickPos
        cb.set_ticks(tickpos_data_coords) # In data coords (0..3)
        tick_names = l2_plot_options['tick_names']
        ppl.setp(cb.ax,xticklabels=tick_names)

        # Colourbar title
        cbar_title = l2_plot_options['name']
        cax_title = ppl.setp(cax,title=cbar_title)
        ppl.setp(cax_title,fontsize=10)

        # Redraw the figure
        canvas.draw()

        # Write the figure to file
        canvas.print_figure(png_file,dpi=dpi)
        LOG.info("Writing to {}...".format(png_file))


    def plot_L2_Map(self,lat,lon,data,data_mask,pngName,dataset,**plot_options):

        l2_plot_options = geocat_l2_data.Plot_Options.data[dataset]

        # Copy the plot options to local variables
        title         = plot_options['title']
        cbar_title    = plot_options['cbar_title']
        units         = plot_options['units']
        stride        = plot_options['stride']
        lat_0         = plot_options['lat_0']
        lon_0         = plot_options['lon_0']
        plotMin       = plot_options['plotMin']
        plotMax       = plot_options['plotMax']
        plotLims      = plot_options['plotLims']
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
                    format(dataset))
            return -1

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*5,scale*5))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.05, 0.15, 0.9, 0.8  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

        # Granule axis title
        ax_title = ppl.setp(ax,title=title)
        ppl.setp(ax_title,fontsize=12)
        ppl.setp(ax_title,family="sans-serif")

        # Setup the map
        m = Basemap(projection='geos',lon_0=lon_0,ax=ax,fix_aspect=True,resolution=map_res)

        x,y=m(lon[::stride,::stride],lat[::stride,::stride])

        coastline_color = l2_plot_options['coastline_color']
        country_color = l2_plot_options['country_color']

        m.drawcoastlines(ax=ax,color=coastline_color,linewidth = 0.5)
        m.drawcountries(ax=ax,color=country_color)
        m.fillcontinents(color='0.85',zorder=0)
        m.drawparallels(np.arange( -90, 91,30), color = '0.25', 
                linewidth = 0.5)
        m.drawmeridians(np.arange(-180,180,30), color = '0.25', 
                linewidth = 0.5)

        LOG.debug('data.shape = {}'.format(data.shape))
        LOG.debug('data_mask.shape = {}'.format(data_mask.shape))
        data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

        vmin,vmax = plotLims[0],plotLims[-1]

        if doScatterPlot:
            cs = m.scatter(x,y,s=pointSize,c=data,axes=ax,edgecolors='none',
                    vmin=vmin,vmax=vmax,cmap=cmap)
        else:
            cs = m.pcolor(x,y,data,axes=ax,edgecolors='none',antialiased=False,
                    vmin=vmin,vmax=vmax,cmap=cmap)

        #txt = ax.set_title(title,fontsize=11)

        ppl.setp(ax.get_xticklines(),visible=False)
        ppl.setp(ax.get_yticklines(),visible=False)
        ppl.setp(ax.get_xticklabels(), visible=False)
        ppl.setp(ax.get_yticklabels(), visible=False)

        # add a colorbar axis
        cax_rect = [0.05 , 0.05, 0.9 , 0.05 ] # [left,bottom,width,height]
        cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

        # Plot the colorbar.
        cb = fig.colorbar(cs, cax=cax, orientation='horizontal')
        ppl.setp(cax.get_xticklabels(),fontsize=9)
        ppl.setp(cax.get_xticklines(),visible=False)

        # Colourbar title
        cax_title = ppl.setp(cax,title=cbar_title)
        ppl.setp(cax_title,fontsize=10)

        # Redraw the figure
        canvas.draw()

        # Write the figure to file
        canvas.print_figure(pngName,dpi=dpi)
        LOG.info("Writing image file {}".format(pngName))


    def plot_L2_Map_discrete(self,lat,lon,data,data_mask,pngName,dataset,**plot_options):
            
        l2_plot_options = geocat_l2_data.Plot_Options.data[dataset]

        # Copy the plot options to local variables
        title         = plot_options['title']
        cbar_title    = plot_options['cbar_title']
        units         = plot_options['units']
        stride        = plot_options['stride']
        lat_0         = plot_options['lat_0']
        lon_0         = plot_options['lon_0']
        plotMin       = plot_options['plotMin']
        plotMax       = plot_options['plotMax']
        plotLims      = plot_options['plotLims']
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
                    format(dataset))
            return -1

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*5,scale*5))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.05, 0.15, 0.9, 0.8  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

        # Granule axis title
        ax_title = ppl.setp(ax,title=title)
        ppl.setp(ax_title,fontsize=12)
        ppl.setp(ax_title,family="sans-serif")

        fill_colours = l2_plot_options['fill_colours']
        cmap = ListedColormap(fill_colours)

        numCats = np.array(fill_colours).size
        numBounds = numCats + 1

        tickPos = np.arange(float(numBounds))/float(numCats)
        tickPos = tickPos[0 :-1] + tickPos[1]/2.

        # Setup the map
        m = Basemap(projection='geos',lon_0=lon_0,ax=ax,fix_aspect=True,resolution=map_res)

        x,y=m(lon[::stride,::stride],lat[::stride,::stride])

        coastline_color = l2_plot_options['coastline_color']
        country_color = l2_plot_options['country_color']

        m.drawcoastlines(ax=ax,color=coastline_color,linewidth = 0.5)
        m.drawcountries(ax=ax,color=country_color)
        m.fillcontinents(color='0.85',zorder=0)
        m.drawparallels(np.arange( -90, 91,30), color = '0.25', 
                linewidth = 0.5)
        m.drawmeridians(np.arange(-180,180,30), color = '0.25', 
                linewidth = 0.5)

        LOG.debug('data.shape = {}'.format(data.shape))
        LOG.debug('data_mask.shape = {}'.format(data_mask.shape))
        data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

        data_values = l2_plot_options['values']
        vmin,vmax = data_values[0],data_values[-1]

        if doScatterPlot:
            cs = m.scatter(x,y,s=pointSize,c=data,axes=ax,edgecolors='none',
                    vmin=vmin,vmax=vmax,cmap=cmap)
        else:
            cs = m.pcolor(x,y,data,axes=ax,edgecolors='none',antialiased=False,
                    vmin=vmin,vmax=vmax,cmap=cmap)

        #txt = ax.set_title(title,fontsize=11)

        ppl.setp(ax.get_xticklines(),visible=False)
        ppl.setp(ax.get_yticklines(),visible=False)
        ppl.setp(ax.get_xticklabels(), visible=False)
        ppl.setp(ax.get_yticklabels(), visible=False)

        # add a colorbar axis
        cax_rect = [0.05 , 0.05, 0.9 , 0.05 ] # [left,bottom,width,height]
        cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

        # Plot the colorbar.
        cb = fig.colorbar(cs, cax=cax, orientation='horizontal')
        ppl.setp(cax.get_xticklabels(),fontsize=9)
        ppl.setp(cax.get_xticklines(),visible=False)

        # Set the colourbar tick locations and ticklabels
        #ppl.setp(cb.ax,xticks=l2_plot_options.ViirsCMTickPos) # In colorbar axis coords (0..1)
        tickpos_data_coords = vmax*tickPos
        cb.set_ticks(tickpos_data_coords) # In data coords (0..3)
        tick_names = l2_plot_options['tick_names']
        ppl.setp(cb.ax,xticklabels=tick_names)

        # Colourbar title
        cbar_title = l2_plot_options['name']
        cax_title = ppl.setp(cax,title=cbar_title)
        ppl.setp(cax_title,fontsize=10)

        # Redraw the figure
        canvas.draw()

        # Write the figure to file
        canvas.print_figure(pngName,dpi=dpi)
        LOG.info("Writing image file {}".format(pngName))


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

    defaults = {
                'input_file':None,
                'dataset':'baseline_cmask_goes_nop_cloud_mask',
                'stride':1,
                'lat_0':None,
                'lon_0':None,
                'plotMin'  : None,
                'plotMax'  : None,
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

    description = '''Create a plot of a level-2 dataset from a geocat HDF4 file.'''

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
                      help='''The fully qualified path to a single geocat level-2 HDF4 input file.'''
                      )

    parser.add_argument(
                      action="store",
                      dest="dataset",
                      default=defaults["dataset"],
                      type=str,
                      help='''The geocat level-2 dataset to plot.
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

    parser.add_argument('--cmap',
                      action="store",
                      dest="cmap",
                      default=defaults["cmap"],
                      type=str,
                      help="""The matplotlib colormap to use. 
                      [default: '{}']""".format(defaults["cmap"])
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

    input_file = options.input_file
    dataset = options.dataset
    stride = options.stride
    lat_0  = options.lat_0
    lon_0  = options.lon_0
    latMin = options.latMin
    lonMin = options.lonMin
    latMax = options.latMax
    lonMax = options.lonMax
    plotMin = options.plotMin
    plotMax = options.plotMax
    doScatterPlot = options.doScatterPlot
    unnavigated = options.unnavigated
    list_datasets = options.list_datasets
    pointSize = options.pointSize
    map_res = options.map_res
    cmap = options.cmap
    cbar_title = options.cbar_title
    output_file  = options.output_file
    outputFilePrefix  = options.outputFilePrefix
    dpi = options.dpi

    # Create and populate the GOES-L2 object

    goes_l2_obj = GOES_L2(input_file)

    lats = goes_l2_obj.Dataset(goes_l2_obj,'pixel_latitude').dset
    lons = goes_l2_obj.Dataset(goes_l2_obj,'pixel_longitude').dset
    sat_zenith_angle = goes_l2_obj.Dataset(goes_l2_obj,'pixel_satellite_zenith_angle').dset

    # If we want to list the datasets, do that here and exit
    if list_datasets:
        LOG.info('Datasets in {}:'.format(input_file))
        for dsets in goes_l2_obj.datanames:
            print "\t{}".format(dsets)
        goes_l2_obj.close()
        sys.exit(1)

    # Read in the desired dataset
    try:
        data_obj = goes_l2_obj.Dataset(goes_l2_obj,dataset)
    except :
        LOG.error('"{}" is not a valid dataset in {}, aborting.'.format(dataset,input_file))
        goes_l2_obj.close()
        sys.exit(1)

    LOG.info('Subsatellite_Longitude = {}'.format(goes_l2_obj.attrs['Subsatellite_Longitude']))
    lon_0 = goes_l2_obj.attrs['Subsatellite_Longitude'] if lon_0==None else lon_0

    # Use the solar zenith angle to mask off-disk pixels...
    data = ma.masked_array(data_obj.dset,mask=sat_zenith_angle.mask)

    if ma.is_masked(data):
        if data.mask.shape == ():
            data_mask = np.ones(data.shape,dtype='bool')
        else:
            data_mask = data.mask
    else: 
        data_mask = np.zeros(data.shape,dtype='bool')

    plot_title = "{}".format(path.basename(input_file))
    #if cbar_title==None:
        #cbar_title = "{} ({})".format(data_obj.dataname,data_obj.attrs['units']) 

    input_file = path.basename(input_file)

    goes_l2_obj.close()
    

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

    l2_plot_options = geocat_l2_data.Plot_Options.data[dataset]

    # Default plot options
    plot_options = {}
    plot_options['title'] = plot_title
    plot_options['cbar_title'] = l2_plot_options['name'] if cbar_title==None else cbar_title
    plot_options['units'] = cbar_title
    plot_options['stride'] = stride
    plot_options['lat_0'] = lat_0
    plot_options['lon_0'] = lon_0
    plot_options['latMin'] = latMin
    plot_options['lonMin'] = lonMin
    plot_options['latMax'] = latMax
    plot_options['lonMax'] = lonMax
    plot_options['plotMin'] = l2_plot_options['values'][0] if plotMin==None else plotMin
    plot_options['plotMax'] = l2_plot_options['values'][-1] if plotMax==None else plotMax
    plot_options['plotLims'] = [plot_options['plotMin'],plot_options['plotMax']]
    plot_options['map_res'] = map_res
    #plot_options['cmap'] = l2_plot_options['cmap'] if cmap==None else cmap
    plot_options['scatterPlot'] = doScatterPlot
    plot_options['pointSize'] = pointSize
    plot_options['dpi'] = dpi

    if cmap == None:
        plot_options['cmap'] = l2_plot_options['cmap']
    else :
        try:
            cmap = getattr(cm,options.cmap)
        except AttributeError:
            LOG.warning('Colormap {} does not exist, falling back to Spectral_r'.format(options.cmap))
            cmap = getattr(cm,'Spectral_r')

    is_l2_discrete = l2_plot_options['discrete']

    # Create the plot
    if unnavigated :
        if is_l2_discrete :
            goes_l2_obj.plot_L2_discrete(data,data_mask,output_file,dataset,**plot_options)
        else :
            goes_l2_obj.plot_L2(data,data_mask,output_file,dataset,**plot_options)
    else :
        if is_l2_discrete :
            goes_l2_obj.plot_L2_Map_discrete(lats,lons,data,data_mask,output_file,dataset,**plot_options)
        else :
            goes_l2_obj.plot_L2_Map(lats,lons,data,data_mask,output_file,dataset,**plot_options)

    return 0


if __name__=='__main__':
    sys.exit(main())  
