#!/usr/bin/env python
# encoding: utf-8
"""
geocat_l1_data.py

Purpose: Provide required data for GEOCAT Level-1 products.

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2015-04-07.
Copyright (c) 2012-2013 University of Wisconsin Regents. All rights reserved.

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

import numpy as np
import time

from matplotlib import cm as cm
import geocat_colormaps as g1_cmaps

class Plot_Options:
    """
    This class contains static data for the interpretation of the various  
    geocat discrete and continuous datasets.
    """

    #satellite_cmap_dict={
            #'goes_15_imager_channel_1_reflectance': cm.gray,
            #'goes_15_imager_channel_2_brightness_temperature': CPD.cmap_ice_water,
            #'goes_15_imager_channel_2_emissivity': CPD.cmap_ice_water,
            #'goes_15_imager_channel_2_reflectance': cm.gray,
            #'goes_15_imager_channel_3_brightness_temperature': CPD.cmap_ice_water,
            #'goes_15_imager_channel_4_brightness_temperature': CPD.cmap_ice_water,
            #'goes_15_imager_channel_6_brightness_temperature': CPD.cmap_ice_water,
            #'pixel_ecosystem_type': CPD.cmap_ice_water,
            #'pixel_latitude': CPD.cmap_ice_water,
            #'pixel_longitude': CPD.cmap_ice_water,
            #'pixel_relative_azimuth_angle': CPD.cmap_ice_water,
            #'pixel_satellite_zenith_angle': CPD.cmap_ice_water,
            #'pixel_solar_zenith_angle': CPD.cmap_ice_water,
            #'pixel_surface_type': CPD.cmap_ice_water
    #}

    #geocat_cmap_dict={
            #'channel_14_brightness_temperature': CPD.cmap_ice_water,
            #'channel_16_brightness_temperature': CPD.cmap_ice_water,
            #'channel_2_reflectance': cm.gray,
            #'channel_7_brightness_temperature': CPD.cmap_ice_water,
            #'channel_7_emissivity': CPD.cmap_ice_water,
            #'channel_7_reflectance': cm.gray,
            #'channel_9_brightness_temperature': CPD.cmap_ice_water,
            #'pixel_ecosystem_type': CPD.cmap_ice_water,
            #'pixel_latitude': CPD.cmap_ice_water,
            #'pixel_longitude': CPD.cmap_ice_water,
            #'pixel_relative_azimuth_angle': CPD.cmap_ice_water,
            #'pixel_satellite_zenith_angle': CPD.cmap_ice_water,
            #'pixel_solar_zenith_angle': CPD.cmap_ice_water,
            #'pixel_surface_type': CPD.cmap_ice_water
    #}

    data = {}

    data['channel_2_reflectance'] = {
                'name':'channel_2_reflectance',
                'discrete':False,
                'values':[0.,100.],
                'coastline_color':'white',
                'country_color':'white',
                'cmap':cm.gray,
                'n_levels':256
                }
    data['channel_7_reflectance'] = {
                'name':'channel_7_reflectance',
                'discrete':False,
                'values':[0.,100.],
                'coastline_color':'white',
                'country_color':'white',
                'cmap':cm.gray,
                'n_levels':256
                }
    data['channel_7_emissivity'] = {
                'name':'channel_7_emissivity',
                'discrete':False,
                'values':[0.,100.],
                'coastline_color':'black',
                'country_color':'black',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['channel_7_brightness_temperature'] = {
                'name':'channel_7_brightness_temperature',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'black',
                'country_color':'black',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['channel_9_brightness_temperature'] = {
                'name':'channel_9_brightness_temperature',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'black',
                'country_color':'black',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['channel_14_brightness_temperature'] = {
                'name':'channel_14_brightness_temperature',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'black',
                'country_color':'black',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['channel_16_brightness_temperature'] = {
                'name':'channel_16_brightness_temperature',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'black',
                'country_color':'black',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }

    data['imager_channel_1_reflectance'] = data['channel_2_reflectance']
    data['imager_channel_2_brightness_temperature'] = data['channel_7_brightness_temperature']
    data['imager_channel_2_emissivity'] = data['channel_7_emissivity']
    data['imager_channel_2_reflectance'] = data['channel_7_reflectance']
    data['imager_channel_3_brightness_temperature'] = data['channel_9_brightness_temperature']
    data['imager_channel_4_brightness_temperature'] = data['channel_14_brightness_temperature']
    data['imager_channel_6_brightness_temperature'] = data['channel_16_brightness_temperature']

    data['pixel_ecosystem_type'] = {
                'name':'pixel_ecosystem_type',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_latitude'] = {
                'name':'pixel_latitude',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_longitude'] = {
                'name':'pixel_longitude',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_relative_azimuth_angle'] = {
                'name':'pixel_relative_azimuth_angle',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_satellite_zenith_angle'] = {
                'name':'pixel_satellite_zenith_angle',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_solar_zenith_angle'] = {
                'name':'pixel_solar_zenith_angle',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_surface_type'] = {
                'name':'pixel_surface_type',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
