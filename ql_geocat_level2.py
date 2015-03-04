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

#from cspp_quicklooks.core import ql_common as qlc

'''
run -e geo_L1.py
geo_L1_obj = GOES_L1('geocatL1.GOES-13.2015001.123000.hdf')
LAT = geo_L1_obj.Dataset(geo_L1_obj,'pixel_latitude')
LON = geo_L1_obj.Dataset(geo_L1_obj,'pixel_longitude')

channel_14_brightness_temperature = geo_L1_obj.Dataset(geo_L1_obj,'channel_14_brightness_temperature')
channel_16_brightness_temperature = geo_L1_obj.Dataset(geo_L1_obj,'channel_16_brightness_temperature')
channel_2_reflectance = geo_L1_obj.Dataset(geo_L1_obj,'channel_2_reflectance')
channel_7_brightness_temperature = geo_L1_obj.Dataset(geo_L1_obj,'channel_7_brightness_temperature')
channel_7_emissivity = geo_L1_obj.Dataset(geo_L1_obj,'channel_7_emissivity')
channel_7_reflectance = geo_L1_obj.Dataset(geo_L1_obj,'channel_7_reflectance')
channel_9_brightness_temperature = geo_L1_obj.Dataset(geo_L1_obj,'channel_9_brightness_temperature')

calibration_solar_constant = geo_L1_obj.Dataset(geo_L1_obj,'calibration_solar_constant')

geo_L1_obj.plot_L1(channel_14_brightness_temperature,cmap=cm.Spectral_r)
geo_L1_obj.plot_L1(channel_16_brightness_temperature,cmap=cm.Spectral_r)
geo_L1_obj.plot_L1(channel_2_reflectance)
geo_L1_obj.plot_L1(channel_7_brightness_temperature,cmap=cm.Spectral_r)
geo_L1_obj.plot_L1(channel_7_emissivity)
geo_L1_obj.plot_L1(channel_7_reflectance)
geo_L1_obj.plot_L1(channel_9_brightness_temperature,cmap=cm.Spectral_r)


['bc1_planck',
 'bc2_planck',
 'calibration_offset',
 'calibration_slope',
 'calibration_slope_degrade',
 'calibration_solar_constant',
 'channel_14_brightness_temperature',
 'channel_16_brightness_temperature',
 'channel_2_reflectance',
 'channel_7_brightness_temperature',
 'channel_7_emissivity',
 'channel_7_reflectance',
 'channel_9_brightness_temperature',
 'channel_wavenumber',
 'fk1_planck',
 'fk2_planck',
 'pixel_ecosystem_type',
 'pixel_latitude',
 'pixel_longitude',
 'pixel_relative_azimuth_angle',
 'pixel_satellite_zenith_angle',
 'pixel_solar_zenith_angle',
 'pixel_surface_type',
 'scan_line_time']
'''


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


    def plot_L1(self,data_obj,png_file=None,cmap=cm.gray,plot_title=None,cbar_title=None):

        data = data_obj.dset
        plot_title = data_obj.dataname if plot_title==None else plot_title
        cbar_title = data_obj.attrs['units'] if cbar_title==None else cbar_title
        png_file = '{}-{}.png'.format(self.l1_file,data_obj.dataname) if png_file==None else png_file

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*5,scale*5))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

        # Granule axis title
        ax_title = ppl.setp(ax,title=plot_title)
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
