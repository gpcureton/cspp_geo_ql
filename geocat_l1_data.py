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

import logging
from copy import copy
from matplotlib import cm as cm
import geocat_colormaps as g1_cmaps

# every module should have a LOG object
LOG = logging.getLogger(__file__)

class Satellite(object):

    def factory(type):
        if type == "GOES_NOP": 
            return GOES_NOP()
        if type == "Himawari": 
            return Himawari()
        assert 0, "Bad satellite creation: " + type
 
    factory = staticmethod(factory)

    refl_dict = {
                'quantity':'reflectance',
                'discrete':False,
                'values':[0.,100.],
                'mask_ranges':[],
                'logscale':False,
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.gray,
                'n_levels':256
                }
    emis_dict = {
                'quantity':'emissivity',
                'discrete':False,
                'values':[0.,30.],
                'mask_ranges':[],
                'logscale':False,
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.gray,
                'n_levels':256
                }
    bt_dict  =  {
                'quantity':'brightness temperature',
                'discrete':False,
                'values':[None,None],
                'mask_ranges':[],
                'logscale':False,
                'units': 'K',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.gray_r,
                'n_levels':256
                }

    common_data = {}
    common_data['pixel_ecosystem_type'] = {
                'name':'pixel ecosystem type',
                'quantity':'ecosystem type',
                'discrete':False,
                'values':[0,101],
                #'mask_values':[],
                'mask_ranges':[],
                'logscale':False,
                'units': None,
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                #'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'cmap':g1_cmaps.viridis,
                'n_levels':256
                }
    common_data['pixel_latitude'] = {
                'name':'pixel latitude',
                'quantity':'latitude',
                'discrete':False,
                'values':[None,None],
                'mask_ranges':[],
                'logscale':False,
                'units': None,
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g1_cmaps.viridis,
                'n_levels':256
                }
    common_data['pixel_longitude'] = {
                'name':'pixel longitude',
                'quantity':'longitude',
                'discrete':False,
                'values':[None,None],
                'mask_ranges':[],
                'logscale':False,
                'units': 'degrees',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g1_cmaps.viridis,
                'n_levels':256
                }
    common_data['pixel_relative_azimuth_angle'] = {
                'name':'pixel relative azimuth angle',
                'quantity':'relative azimuth angle',
                'discrete':False,
                'values':[None,None],
                'mask_ranges':[],
                'logscale':False,
                'units': 'degrees',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g1_cmaps.viridis,
                'n_levels':256
                }
    common_data['pixel_satellite_zenith_angle'] = {
                'name':'pixel satellite zenith angle',
                'quantity':'satellite zenith angle',
                'discrete':False,
                'values':[None,None],
                'mask_ranges':[],
                'logscale':False,
                'units': 'degrees',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g1_cmaps.viridis,
                'n_levels':256
                }
    common_data['pixel_solar_zenith_angle'] = {
                'name':'pixel solar zenith angle',
                'quantity':'solar zenith angle',
                'discrete':False,
                'values':[None,None],
                'mask_ranges':[],
                'logscale':False,
                'units': 'degrees',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g1_cmaps.viridis,
                'n_levels':256
                }
    common_data['pixel_surface_type'] = {
                'name':'pixel surface type',
                'quantity':'surface type',
                'discrete':False,
                'values':[0,13],
                #'mask_values':[],
                'mask_ranges':[],
                'logscale':False,
                'units': None,
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                #'cmap':g1_cmaps.geocat_colormaps.get_cmap_ice_water(),
                'cmap':g1_cmaps.viridis,
                'n_levels':256
                }
    common_data['unknown'] = {
                'name':None,
                'quantity':None,
                'discrete':False,
                'values':[None,None],
                'mask_ranges':[],
                'logscale':False,
                'units': None,
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g1_cmaps.viridis,
                'n_levels':256
                }


class GOES_NOP(Satellite):

    def __init__(self, *args, **kwargs):
        Satellite.__init__(self, *args, **kwargs)

        geocat_reflectance_channels = [2, 7]
        geocat_emissivity_channels = [7]
        geocat_btemperature_channels = [7, 9, 14, 16]

        data = {}

        for chan in geocat_reflectance_channels:
            data['channel_{}_reflectance'.format(chan)] = copy(Satellite.refl_dict)
            data['channel_{}_reflectance'.format(chan)]['name'] = 'Channel {} reflectance'.format(chan)
            data['channel_{}_reflectance'.format(chan)]['dset_name'] = 'channel_{}_reflectance'.format(chan)
        for chan in geocat_emissivity_channels:
            data['channel_{}_emissivity'.format(chan)] = copy(Satellite.emis_dict)
            data['channel_{}_emissivity'.format(chan)]['name'] = 'Channel {} emissivity'.format(chan)
            data['channel_{}_emissivity'.format(chan)]['dset_name'] = 'channel_{}_emissivity'.format(chan)
        for chan in geocat_btemperature_channels:
            data['channel_{}_brightness_temperature'.format(chan)] = copy(Satellite.bt_dict)
            data['channel_{}_brightness_temperature'.format(chan)]['name'] = 'Channel {} brightness temperature'.format(chan)
            data['channel_{}_brightness_temperature'.format(chan)]['dset_name'] = 'channel_{}_brightness_temperature'.format(chan)


        data.update(Satellite.common_data)

        for common_dset in Satellite.common_data.keys():
            data[common_dset]['dset_name'] = common_dset


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

        self.data = data


class Himawari(Satellite):

    def __init__(self, *args, **kwargs):
        Satellite.__init__(self, *args, **kwargs)

        geocat_reflectance_channels = [1, 2, 3, 4, 5, 6, 7]
        geocat_emissivity_channels = [7]
        geocat_btemperature_channels = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

        data = {}

        for chan in geocat_reflectance_channels:
            data['ahi_channel_{}_reflectance'.format(chan)] = copy(Satellite.refl_dict)
            data['ahi_channel_{}_reflectance'.format(chan)]['name'] = 'Channel {} reflectance'.format(chan)
            data['ahi_channel_{}_reflectance'.format(chan)]['dset_name'] = 'ahi_channel_{}_reflectance'.format(chan)
        for chan in geocat_emissivity_channels:
            data['ahi_channel_{}_emissivity'.format(chan)] = copy(Satellite.emis_dict)
            data['ahi_channel_{}_emissivity'.format(chan)]['name'] = 'Channel {} emissivity'.format(chan)
            data['ahi_channel_{}_emissivity'.format(chan)]['dset_name'] = 'ahi_channel_{}_emissivity'.format(chan)
        for chan in geocat_btemperature_channels:
            data['ahi_channel_{}_brightness_temperature'.format(chan)] = copy(Satellite.bt_dict)
            data['ahi_channel_{}_brightness_temperature'.format(chan)]['name'] = 'Channel {} brightness temperature'.format(chan)
            data['ahi_channel_{}_brightness_temperature'.format(chan)]['dset_name'] = 'ahi_channel_{}_brightness_temperature'.format(chan)

        data.update(Satellite.common_data)

        for common_dset in Satellite.common_data.keys():
            data[common_dset]['dset_name'] = common_dset

        self.data = data


class Navigation(Satellite):

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
