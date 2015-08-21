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

# every module should have a LOG object
LOG = logging.getLogger(__file__)


class GOES_L1():

    def __init__(self,l1_file):
        pass


    def plot_L1(self,data,data_mask,png_file,dataset,plot_style_options):

        l1_dataset_options = geocat_l1_data.Dataset_Options.data[dataset]

        # Copy the plot options to local variables
        title         = plot_style_options['title']
        cbar_title    = plot_style_options['cbar_title']
        units         = plot_style_options['units']
        stride        = plot_style_options['stride']
        plotMin       = plot_style_options['plotMin']
        plotMax       = plot_style_options['plotMax']
        plotLims      = plot_style_options['plotLims']
        cmap          = plot_style_options['cmap']
        dpi           = plot_style_options['dpi']

        '''
        Plot the input dataset in in native data coordinates
        '''

        # If our data is all missing, return
        if (np.sum(data_mask) == data.size):
            LOG.warn("Entire {} dataset is missing, aborting".\
                    format(cbar_title))
            return -1

        LOG.info("{} data array has shape {}".format(title,data.shape))
        data_aspect = float(data.shape[1])/float(data.shape[0])
        LOG.info("{} data array has aspect {}".format(title,data_aspect))

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*5*data_aspect,1.1*scale*5))
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

        LOG.debug('plotLims = {},{}'.format(plotLims[0],plotLims[1]))
        vmin,vmax = plotLims[0],plotLims[1]

        ax.set_aspect(data_aspect)

        #im = ax.imshow(data,interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax,cmap=cmap)
        im = ax.imshow(data,interpolation='nearest',vmin=vmin,vmax=vmax,cmap=cmap)

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


    def plot_L1_Map(self,lat,lon,data,data_mask,png_file,dataset,
            plot_nav_options,plot_style_options):

        l1_dataset_options = geocat_l1_data.Dataset_Options.data[dataset]

        # Copy the plot options to local variables
        llcrnrx        = plot_nav_options['llcrnrx']
        llcrnry        = plot_nav_options['llcrnry']
        urcrnrx        = plot_nav_options['urcrnrx']
        urcrnry        = plot_nav_options['urcrnry']

        lon_0         = plot_nav_options['lon_0']

        title         = plot_style_options['title']
        cbar_title    = plot_style_options['cbar_title']
        units         = plot_style_options['units']
        stride        = plot_style_options['stride']
        plotMin       = plot_style_options['plotMin']
        plotMax       = plot_style_options['plotMax']
        plotLims      = plot_style_options['plotLims']
        map_res       = plot_style_options['map_res']
        cmap          = plot_style_options['cmap']
        doScatterPlot = plot_style_options['scatterPlot']
        pointSize     = plot_style_options['pointSize']
        dpi           = plot_style_options['dpi']


        '''
        Plot the input dataset in mapped to particular projection
        '''

        # If our data is all missing, return
        if (np.sum(data_mask) == data.size):
            LOG.warn("Entire {} dataset is missing, aborting".\
                    format(cbar_title))
            return -1

        LOG.info("{} data array has shape {}".format(title,data.shape))
        data_aspect = float(data.shape[1])/float(data.shape[0])
        LOG.info("{} data array has aspect {}".format(title,data_aspect))

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*5,scale*5))
        #fig = Figure(figsize=(scale*5*data_aspect,1.1*scale*5))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.10, 0.15, 0.80, 0.8  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect,axis_bgcolor='black')

        # Granule axis title
        ax_title = ppl.setp(ax,title=title)
        ppl.setp(ax_title,fontsize=12)
        ppl.setp(ax_title,family="sans-serif")

        # Setup the map
        m = Basemap(projection='geos',lon_0=lon_0,ax=ax,fix_aspect=True,resolution=map_res,
                llcrnrx=llcrnrx,
                llcrnry=llcrnry,
                urcrnrx=urcrnrx,
                urcrnry=urcrnry
                )

        x,y=m(lon[::stride,::stride],lat[::stride,::stride])

        coastline_color = l1_dataset_options['coastline_color']
        country_color = l1_dataset_options['country_color']
        meridian_color = l1_dataset_options['meridian_color']

        m.drawcoastlines(ax=ax,color=coastline_color,linewidth = 0.3)
        m.drawcountries(ax=ax,color=country_color,linewidth = 0.2)
        m.drawstates(ax=ax,color=country_color,linewidth = 0.2)
        m.fillcontinents(color='0.',zorder=0)
        m.drawparallels(np.arange( -90, 91,30), color = meridian_color, 
                linewidth = 0.5)
        m.drawmeridians(np.arange(-180,180,30), color = meridian_color, 
                linewidth = 0.5)

        LOG.debug('data.shape = {}'.format(data.shape))
        LOG.debug('data_mask.shape = {}'.format(data_mask.shape))
        data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

        LOG.debug('plotLims = {},{}'.format(plotLims[0],plotLims[1]))
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
        cax_rect = [0.10 , 0.05, 0.8 , 0.05 ] # [left,bottom,width,height]
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
        canvas.print_figure(png_file,dpi=dpi)
        LOG.info("Writing image file {}".format(png_file))


class GOES_HDF4(GOES_L1):

    def __init__(self,l1_file):


        self.l1_file = l1_file

        LOG.debug('Opening {} with GOES_HDF4...'.format(self.l1_file))
        self.file_obj = SD(self.l1_file)
        
        # Dictionary of file object attributes
        self.attrs = self.file_obj.attributes()

        # Dictionary of dataset shape and type attributes
        self.data_dict = self.file_obj.datasets()

        # List of dataset names
        self.datanames = self.data_dict.keys()
        self.datanames.sort()


    class Dataset():

        def __init__(selfd,L1_obj,dataname):

            selfd.dataname = dataname

            selfd.dset_obj = L1_obj.file_obj.select(dataname)

            selfd.attrs = selfd.dset_obj.attributes()
            selfd.dset = ma.masked_equal(selfd.dset_obj.get(),selfd.attrs['_FillValue'])
            
            selfd.dset = selfd.dset * selfd.attrs['scale_factor'] + selfd.attrs['add_offset']

    def close(self):
        LOG.debug('Closing {}...'.format(self.l1_file))
        self.file_obj.end()


class GOES_NetCDF(GOES_L1):

    def __init__(self,l1_file):


        self.l1_file = l1_file

        LOG.debug('Opening {} with GOES_NetCDF...'.format(self.l1_file))
        self.file_obj = Dataset(self.l1_file)

        # Dictionary of file object attributes
        self.attrs = {}
        for attr_key in self.file_obj.ncattrs():
            self.attrs[attr_key] = getattr(self.file_obj,attr_key)
        
        # Ordered dictionary of dataset objects
        self.data_dict = self.file_obj.variables

        # List of dataset names
        self.datanames = self.data_dict.keys()
        self.datanames.sort()


    class Dataset():

        def __init__(selfd,L1_obj,dataname):

            selfd.dataname = dataname

            selfd.dset_obj = L1_obj.file_obj.variables[dataname]

            selfd.attrs = {}
            for attr_key in selfd.dset_obj.ncattrs():
                selfd.attrs[attr_key] = getattr(selfd.dset_obj,attr_key)
            
            selfd.dset = ma.masked_equal(selfd.dset_obj[:],selfd.attrs['_FillValue'])
            #selfd.dset = selfd.dset * selfd.attrs['scale_factor'] + selfd.attrs['add_offset']

    def close(self):
        LOG.debug('Closing {}...'.format(self.l1_file))
        self.file_obj.close()


def _tuple2args(parms):
    s = ' '.join( '+%s=%s' % (k,v) for (k,v) in parms )
    return s.encode('ascii')


def set_plot_navigation_proj(lats,lons,goes_l1_obj, options):
    """
    This method determines the appropriate plotting viewport values.
    """
    # Here are some attributes from the geocat Level-1 files...
    '''
    :Platform_Name = "GOES-13" ;
    :Subsatellite_Longitude = -74.9999542236328 ;
    :Latitude_Range = -80.8643646240234, 80.925651550293 ;
    :Longitude_Range = -156.056655883789, 6.18542814254761 ;
    :Earth-Sun_Distance = 1.0122344493866 ;
    :Satellite_Height = 42166112. ;
    :Sensor_ID = 180L ;
    '''
    from mpl_toolkits.basemap.pyproj import Proj

    nrows,ncols = lats.shape[0],lats.shape[1]
    nrows_div2,ncols_div2 = np.floor_divide(nrows,2),np.floor_divide(ncols,2)
    LOG.info('lats.shape = {}'.format(lats.shape))
    LOG.info('nrows,ncols = ({},{})'.format(nrows,ncols))
    LOG.info('nrows_div2,ncols_div2 = ({},{})'.format(nrows_div2,ncols_div2))

    corners = [(0,0), (0,-1), (-1,-1), (-1,0)]

    for crnr,crnr_idx in zip(range(4),corners):
        LOG.info("({}) (lon,lat):({},{})".format(crnr,lons[crnr_idx],lats[crnr_idx]))


    subsatellite_longitude = goes_l1_obj.attrs['Subsatellite_Longitude']
    LOG.info(' File subsatellite_Longitude = {}'.format(subsatellite_longitude))

    plot_nav_options = {}

    print "\n>>>>>>> Using Proj\n"
    p1 = Proj("+proj=geos +h=35774290 +a= 6378137 +b= 6378137 +lon_0=-75 +units=meters +no_defs")

    # Western edge
    x,y = p1(lons[:,0],lats[:,0])
    west_x = np.average(x.compressed())
    # Eastern edge
    x,y = p1(lons[:,-1],lats[:,-1])
    east_x = np.average(x.compressed())
    # Northern edge
    x,y = p1(lons[0,:],lats[0,:])
    north_y = np.average(y.compressed())
    # Southern edge
    x,y = p1(lons[-1,:],lats[-1,:])
    south_y = np.average(y.compressed())

    corner_x_fallback = [west_x,east_x,east_x,west_x]
    corner_y_fallback = [north_y,north_y,south_y,south_y]

    LOG.info("(west_x ,east_x ):({:10.1f},{:10.1f})".format(west_x,east_x))
    LOG.info("(north_y,south_y):({:10.1f},{:10.1f})".format(north_y,south_y))

    crnr_x_names_proj = ['ulcrnrx_proj','urcrnrx_proj','lrcrnrx_proj','llcrnrx_proj']
    crnr_y_names_proj = ['ulcrnry_proj','urcrnry_proj','lrcrnry_proj','llcrnry_proj']

    # The maximum extent of the full disk in the x and y directions
    plot_nav_options['extent_x'] = p1(subsatellite_longitude+81.,0)[0]*2.
    plot_nav_options['extent_y'] = p1(subsatellite_longitude,81.)[1]*2.
    LOG.info("(extent_x,extent_y):({},{})".format(plot_nav_options['extent_x'],plot_nav_options['extent_y']))

    # Generate the corner-origin coordinates from Basemap object...
    for crnr in range(4):
        crnr_idx = corners[crnr]
        lon,lat = lons[crnr_idx], lats[crnr_idx]
        if ma.is_masked(lon) and ma.is_masked(lat):
            x,y = corner_x_fallback[crnr],corner_y_fallback[crnr]
        else:
            x,y = p1(lon, lat)
        #if crnr_x_names_proj[crnr] is not None:
        plot_nav_options[crnr_x_names_proj[crnr]] = x
        #if crnr_y_names_proj[crnr] is not None:
        plot_nav_options[crnr_y_names_proj[crnr]] = y
        LOG.info(
                "({}) (lon,lat):({},{}), (x,y): ({:10.1f},{:10.1f})".format(crnr,lon,lat,x,y))


    print ""
    for x_name,y_name in zip(crnr_x_names_proj,crnr_y_names_proj):
        #if x_name is not None and y_name is not None:
        print "plot_nav_options['{}','{}'] = {:10.1f},{:10.1f}".format(x_name,y_name,plot_nav_options[x_name],plot_nav_options[y_name])

    # Default plot options
    plot_nav_options['llcrnrx'] = options.llcrnrx
    plot_nav_options['llcrnry'] = options.llcrnry
    plot_nav_options['urcrnrx'] = options.urcrnrx
    plot_nav_options['urcrnry'] = options.urcrnry


    return plot_nav_options


def set_plot_navigation_bm(lats,lons,goes_l1_obj, options):
    """
    This method determines the appropriate plotting viewport values.
    """

    nrows,ncols = lats.shape[0],lats.shape[1]
    nrows_div2,ncols_div2 = np.floor_divide(nrows,2),np.floor_divide(ncols,2)
    LOG.info('lats.shape = {}'.format(lats.shape))
    LOG.info('nrows,ncols = ({},{})'.format(nrows,ncols))
    LOG.info('nrows_div2,ncols_div2 = ({},{})'.format(nrows_div2,ncols_div2))

    corners = [(0,0), (0,-1), (-1,-1), (-1,0)]

    for crnr,crnr_idx in zip(range(4),corners):
        LOG.info("({}) (lon,lat):({},{})".format(crnr,lons[crnr_idx],lats[crnr_idx]))


    subsatellite_longitude = goes_l1_obj.attrs['Subsatellite_Longitude']
    LOG.info(' File subsatellite_Longitude = {}'.format(subsatellite_longitude))

    print "\n>>>>>>> Using Basemap\n"
    m1 = Basemap(projection='geos',lon_0=subsatellite_longitude,resolution=None)

    plot_nav_options = {}

    # The maximum extent of the full disk in the x and y directions
    plot_nav_options['extent_x'] = m1.urcrnrx
    plot_nav_options['extent_y'] = m1.urcrnry
    LOG.info("(extent_x,extent_y):({},{})".format(plot_nav_options['extent_x'],plot_nav_options['extent_y']))

    # Compute coordinates of sub-satellite point
    x_subsat,y_subsat = m1(subsatellite_longitude,0)

    is_full_disk = False
    if ma.is_masked(lats[corners[0]]) and ma.is_masked(lats[corners[1]]) and \
       ma.is_masked(lats[corners[2]]) and ma.is_masked(lats[corners[3]]):
        is_full_disk = True

    if options.region == "FD":
        is_full_disk = True

    if is_full_disk:

        LOG.info("This image is full disk")

        plot_nav_options['urcrnrx_map'] = plot_nav_options['extent_x'] - x_subsat
        plot_nav_options['urcrnry_map'] = plot_nav_options['extent_y'] - y_subsat
        plot_nav_options['llcrnrx_map'] = 0. - x_subsat
        plot_nav_options['llcrnry_map'] = 0  - y_subsat

    else:
        LOG.info("This image is NOT full disk")

        # Western edge
        x,y = m1(lons[:,0],lats[:,0])
        west_x = np.average(x.compressed())
        # Eastern edge
        x,y = m1(lons[:,-1],lats[:,-1])
        east_x = np.average(x.compressed())
        # Northern edge
        x,y = m1(lons[0,:],lats[0,:])
        north_y = np.average(y.compressed())
        # Southern edge
        x,y = m1(lons[-1,:],lats[-1,:])
        south_y = np.average(y.compressed())

        corner_x_fallback = [west_x,east_x,east_x,west_x]
        corner_y_fallback = [north_y,north_y,south_y,south_y]

        LOG.info("(west_x ,east_x ):({:10.1f},{:10.1f})".format(west_x,east_x))
        LOG.info("(north_y,south_y):({:10.1f},{:10.1f})".format(north_y,south_y))

        crnr_x_names = ['ulcrnrx','urcrnrx','lrcrnrx','llcrnrx']
        crnr_y_names = ['ulcrnry','urcrnry','lrcrnry','llcrnry']

        # Generate the center-origin coordinates from Proj object...
        for crnr in range(4):
            crnr_idx = corners[crnr]
            lon,lat = lons[crnr_idx], lats[crnr_idx]
            if ma.is_masked(lon) and ma.is_masked(lat):
                x,y = corner_x_fallback[crnr],corner_y_fallback[crnr]
            else:
                x,y = m1(lon, lat)
            plot_nav_options[crnr_x_names[crnr]] = x
            plot_nav_options[crnr_y_names[crnr]] = y
            LOG.info(
                    "({}) (lon,lat):({},{}), (x,y): ({:10.1f},{:10.1f})".format(crnr,lon,lat,x,y))


        print ""
        for x_name,y_name in zip(crnr_x_names,crnr_y_names):
            print "plot_nav_options['{}','{}'] = {:10.1f},{:10.1f}".format(x_name,y_name,plot_nav_options[x_name],plot_nav_options[y_name])

        plot_nav_options['ulcrnrx_map'] = plot_nav_options['ulcrnrx'] - x_subsat
        plot_nav_options['ulcrnry_map'] = plot_nav_options['ulcrnry'] - y_subsat
        plot_nav_options['urcrnrx_map'] = plot_nav_options['urcrnrx'] - x_subsat
        plot_nav_options['urcrnry_map'] = plot_nav_options['urcrnry'] - y_subsat
        plot_nav_options['lrcrnrx_map'] = plot_nav_options['lrcrnrx'] - x_subsat
        plot_nav_options['lrcrnry_map'] = plot_nav_options['lrcrnry'] - y_subsat
        plot_nav_options['llcrnrx_map'] = plot_nav_options['llcrnrx'] - x_subsat
        plot_nav_options['llcrnry_map'] = plot_nav_options['llcrnry'] - y_subsat

        print ""
        print "plot_nav_options['ulcrnrx_map','ulcrnry_map'] = {:10.1f},{:10.1f}".format(plot_nav_options['ulcrnrx_map'],plot_nav_options['ulcrnry_map'])
        print "plot_nav_options['urcrnrx_map','urcrnry_map'] = {:10.1f},{:10.1f}".format(plot_nav_options['urcrnrx_map'],plot_nav_options['urcrnry_map'])
        print "plot_nav_options['lrcrnrx_map','lrcrnry_map'] = {:10.1f},{:10.1f}".format(plot_nav_options['lrcrnrx_map'],plot_nav_options['lrcrnry_map'])
        print "plot_nav_options['llcrnrx_map','llcrnry_map'] = {:10.1f},{:10.1f}".format(plot_nav_options['llcrnrx_map'],plot_nav_options['llcrnry_map'])


    # Default plot options
    plot_nav_options['llcrnrx'] = plot_nav_options['llcrnrx_map'] if options.llcrnrx is None else options.llcrnrx * plot_nav_options['extent_x']
    plot_nav_options['llcrnry'] = plot_nav_options['llcrnry_map'] if options.llcrnry is None else options.llcrnry * plot_nav_options['extent_y'] 
    plot_nav_options['urcrnrx'] = plot_nav_options['urcrnrx_map'] if options.urcrnrx is None else options.urcrnrx * plot_nav_options['extent_x'] 
    plot_nav_options['urcrnry'] = plot_nav_options['urcrnry_map'] if options.urcrnry is None else options.urcrnry * plot_nav_options['extent_y'] 

    plot_nav_options['lon_0'] = subsatellite_longitude
    return plot_nav_options


def set_plot_styles(goes_l1_obj, dataset, options):

    l1_dataset_options = geocat_l1_data.Dataset_Options.data[dataset]

    plot_style_options = {}
    plot_style_options['title'] = " {}".format(path.basename(options.input_file))
    plot_style_options['cbar_title'] = l1_dataset_options['name'] if options.cbar_title==None else options.cbar_title
    plot_style_options['units'] = options.cbar_title
    plot_style_options['stride'] = options.stride
    plot_style_options['plotMin'] = l1_dataset_options['values'][0] if options.plotMin==None else options.plotMin
    plot_style_options['plotMax'] = l1_dataset_options['values'][-1] if options.plotMax==None else options.plotMax
    plot_style_options['plotLims'] = [plot_style_options['plotMin'],plot_style_options['plotMax']]
    plot_style_options['map_res'] = options.map_res
    plot_style_options['scatterPlot'] = options.doScatterPlot
    plot_style_options['pointSize'] = options.pointSize
    plot_style_options['dpi'] = options.dpi

    if options.cmap == None:
        plot_style_options['cmap'] = l1_dataset_options['cmap']
    else :
        try:
            plot_style_options['cmap'] = getattr(cm,options.cmap)
        except AttributeError:
            LOG.warning('Colormap {} does not exist, falling back to Spectral_r'.format(options.cmap))
            plot_style_options['cmap'] = getattr(cm,'Spectral_r')

    return plot_style_options


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
                #'lat_0':None,
                'lon_0':None,
                'extent_x':1.,
                'extent_y':1.,
                'offset_x':0.,
                'offset_y':0.,
                'llcrnrx':None,
                'llcrnry':None,
                'urcrnrx':None,
                'urcrnry':None,
                'plotMin'  : None,
                'plotMax'  : None,
                'region' : None,
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

    #parser.add_argument('--lat_0',
                      #action="store",
                      #dest="lat_0",
                      #type=float,
                      #help="Center latitude of plot."
                      #)

    parser.add_argument('--lon_0',
                      action="store",
                      dest="lon_0",
                      type=float,
                      help="Center the plot over this longitude."
                      )

    parser.add_argument('--extent_x',
                      action="store",
                      dest="extent_x",
                      default=defaults["extent_x"],
                      type=float,
                      help='''x-direction extent of the plot viewport in the range [0.,1] (for
                      navigated plots only)[default: {}]'''.format(defaults["extent_x"])
                      )

    parser.add_argument('--extent_y',
                      action="store",
                      dest="extent_y",
                      default=defaults["extent_y"],
                      type=float,
                      help='''y-direction extent of the plot viewport in the range [0.,1] (for
                      navigated plots only)[default: {}]'''.format(defaults["extent_y"])
                      )

    parser.add_argument('--offset_x',
                      action="store",
                      dest="offset_x",
                      default=defaults["offset_x"],
                      type=float,
                      help='''x-direction offset of the plot viewport in the range [-0.5,+0.5] (for
                      navigated plots only)[default: {}]'''.format(defaults["offset_x"])
                      )

    parser.add_argument('--offset_y',
                      action="store",
                      dest="offset_y",
                      default=defaults["offset_y"],
                      type=float,
                      help='''y-direction offset of the plot viewport in the range [-0.5,+0.5] (for
                      navigated plots only)[default: {}]'''.format(defaults["offset_y"])
                      )

    parser.add_argument('--llcrnrx',
                      action="store",
                      dest="llcrnrx",
                      default=defaults["llcrnrx"],
                      type=float,
                      help='''Lower left x-coordinate of the plot in the range [-0.5,+0.5] (for
                      navigated plots only)[default: {}]'''.format(defaults["llcrnrx"])
                      )

    parser.add_argument('--llcrnry',
                      action="store",
                      dest="llcrnry",
                      default=defaults["llcrnry"],
                      type=float,
                      help='''Lower left y-coordinate of the plot in the range [-0.5,+0.5] (for
                      navigated plots only)[default: {}]'''.format(defaults["llcrnry"])
                      )

    parser.add_argument('--urcrnrx',
                      action="store",
                      dest="urcrnrx",
                      default=defaults["urcrnrx"],
                      type=float,
                      help='''Lower left x-coordinate of the plot in the range [-0.5,+0.5] (for . 
                      navigated plots only)[default: {}]'''.format(defaults["urcrnrx"])
                      )

    parser.add_argument('--urcrnry',
                      action="store",
                      dest="urcrnry",
                      default=defaults["urcrnry"],
                      type=float,
                      help='''Lower left x-coordinate of the plot in the range [-0.5,+0.5] (for . 
                      navigated plots only)[default: {}]'''.format(defaults["urcrnry"])
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

    #input_file = options.input_file
    #dataset = options.dataset
    #stride = options.stride
    #lat_0  = options.lat_0
    #lon_0  = options.lon_0
    #llcrnrx = options.llcrnrx
    #llcrnry = options.llcrnry
    #urcrnrx = options.urcrnrx
    #urcrnry = options.urcrnry
    #plotMin = options.plotMin
    #plotMax = options.plotMax
    #doScatterPlot = options.doScatterPlot
    #unnavigated = options.unnavigated
    #list_datasets = options.list_datasets
    #pointSize = options.pointSize
    #map_res = options.map_res
    #cmap = options.cmap
    cbar_title = options.cbar_title
    output_file  = options.output_file
    outputFilePrefix  = options.outputFilePrefix
    dpi = options.dpi

    # Create and populate the GOES-L1 object
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
    LOG.info('options.dataset: {}'.format(options.dataset))
    try:

        LOG.info('options.dataset name: {}'.format(options.dataset))

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

    plot_title = "{}".format(path.basename(options.input_file))
    if options.cbar_title==None:
        cbar_title = "{} ({})".format(data_obj.dataname,data_obj.attrs['units']) 

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

    # Plot options relating to the plotting viewport...
    plot_style_options = set_plot_styles(goes_l1_obj,dataset,options)

    # Create the plot
    if options.unnavigated :
        goes_l1_obj.plot_L1(data,data_mask,output_file,dataset,
                plot_style_options)
    else :
        #plot_nav_options = set_plot_navigation_proj(lats,lons,goes_l1_obj,options)
        plot_nav_options = set_plot_navigation_bm(lats,lons,goes_l1_obj,options)
        goes_l1_obj.plot_L1_Map(lats,lons,data,data_mask,output_file,dataset,
                plot_nav_options,plot_style_options)

    return 0


if __name__=='__main__':
    sys.exit(main())  
