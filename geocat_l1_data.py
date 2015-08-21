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

class Dataset_Options:
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
                'name':'Channel 2 reflectance',
                'quantity':'reflectance',
                'discrete':False,
                'values':[0.,100.],
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.gray,
                'n_levels':256
                }
    data['channel_7_reflectance'] = {
                'name':'Channel 7 reflectance',
                'quantity':'reflectance',
                'discrete':False,
                'values':[0.,75.],
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.gray,
                'n_levels':256
                }
    data['channel_7_emissivity'] = {
                'name':'Channel 7 emissivity',
                'quantity':'emissivity',
                'discrete':False,
                'values':[0.,30.],
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.gray,
                'n_levels':256
                }
    data['channel_7_brightness_temperature'] = {
                'name':'Channel 7 brightness temperature',
                'quantity':'brightness temperature',
                'discrete':False,
                'values':[None,None],
                'units': 'K',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.gray_r,
                'n_levels':256
                }
    data['channel_9_brightness_temperature'] = {
                'name':'Channel 9 brightness temperature',
                'quantity':'brightness temperature',
                'discrete':False,
                'values':[None,None],
                'units': 'K',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.gray_r,
                'n_levels':256
                }
    data['channel_14_brightness_temperature'] = {
                'name':'Channel 14 brightness temperature',
                'quantity':'brightness temperature',
                'discrete':False,
                'values':[None,None],
                'units': 'K',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.gray_r,
                'n_levels':256
                }
    data['channel_16_brightness_temperature'] = {
                'name':'Channel 16 brightness temperature',
                'quantity':'brightness temperature',
                'discrete':False,
                'values':[None,None],
                'units': 'K',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.gray_r,
                'n_levels':256
                }

    data['imager_channel_1_reflectance'] = data['channel_2_reflectance']
    data['imager_channel_2_brightness_temperature'] = data['channel_7_brightness_temperature']
    data['imager_channel_2_emissivity'] = data['channel_7_emissivity']
    data['imager_channel_2_reflectance'] = data['channel_7_reflectance']
    data['imager_channel_3_brightness_temperature'] = data['channel_9_brightness_temperature']
    data['imager_channel_4_brightness_temperature'] = data['channel_14_brightness_temperature']
    data['imager_channel_6_brightness_temperature'] = data['channel_16_brightness_temperature']

    for dsets in ['imager_channel_1_reflectance','imager_channel_2_brightness_temperature',
            'imager_channel_2_emissivity','imager_channel_2_reflectance',
            'imager_channel_3_brightness_temperature','imager_channel_4_brightness_temperature',
            'imager_channel_6_brightness_temperature']:
        data[dsets]['name'] = " ".join(dsets.split('_')[1:])


    data['pixel_ecosystem_type'] = {
                'name':'pixel ecosystem type',
                'discrete':True,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'meridian_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_latitude'] = {
                'name':'pixel latitude',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'meridian_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_longitude'] = {
                'name':'pixel longitude',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'meridian_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_relative_azimuth_angle'] = {
                'name':'pixel relative azimuth angle',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'meridian_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_satellite_zenith_angle'] = {
                'name':'pixel satellite zenith angle',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'meridian_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_solar_zenith_angle'] = {
                'name':'pixel solar zenith angle',
                'discrete':False,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'meridian_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }
    data['pixel_surface_type'] = {
                'name':'pixel surface type',
                'discrete':True,
                'values':[None,None],
                'coastline_color':'white',
                'country_color':'white',
                'meridian_color':'white',
                'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'n_levels':256
                }

    navigation = {}

    navigation['goes_e'] = {
            'Subsatellite_Longitude': -75.,
            'FD':{
                'extent':{'ew':0.3036872 , 'ns':0.3036872},
                'offset':{'ew':0.        , 'ns':0.}
            },
            'CONUS':{
                'extent':{'ew':0.14        , 'ns':0.084},
                'offset':{'ew':-0.040248647, 'ns':0.084625052}
            },
            'MESO':{
                'extent':{'ew':0.028     , 'ns':0.028},
                'offset':{'ew':None      , 'ns':None}
            }

            #'extent':{
                #'FD'   : {'ew':0.3036872 , 'ns':0.3036872},
                #'CONUS': {'ew':0.14      , 'ns':0.084},
                #'MESO' : {'ew':0.028     , 'ns':0.028}
            #},
            #'img_center_offset':{
                #'FD'   : {'ew':0.          , 'ns':0.},
                #'CONUS': {'ew':-0.040248647, 'ns':0.084625052},
                #'MESO' : {'ew':None        , 'ns':None}
            #}
    }

    navigation['goes_w'] = {
            'Subsatellite_Longitude': -137.,
            'FD':{
                'extent':{'ew':0.3036872 , 'ns':0.3036872},
                'offset':{'ew':0.        , 'ns':0.}
            },
            'CONUS':{
                'extent':{'ew':0.14        , 'ns':0.084},
                'offset':{'ew':0.082900064, 'ns':0.083759424}
            },
            'MESO':{
                'extent':{'ew':0.028     , 'ns':0.028},
                'offset':{'ew':None      , 'ns':None}
            }


            #'extent':{
                #'FD'   : {'ew':0.3036872 , 'ns':0.3036872},
                #'CONUS': {'ew':0.14      , 'ns':0.084},
                #'MESO' : {'ew':0.028     , 'ns':0.028}
            #},
            #'img_center_offset':{
                #'FD'   : {'ew':0.         , 'ns':0.},
                #'CONUS': {'ew':0.082900064, 'ns':0.083759424},
                #'MESO' : {'ew':None       , 'ns':None}
            #}
    }
