#!/usr/bin/env python
# encoding: utf-8
"""
geocat_l2_data.py

Purpose: Provide required data for GEOCAT Level-2 products.

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

from copy import copy
from matplotlib import cm as cm
import geocat_colormaps as g2_cmaps
import colormaps as new_cm


class Satellite(object):

    def factory(type):
        if type == "GOES_NOP": 
            return GOES_NOP()
        if type == "Himawari": 
            return Himawari()
        assert 0, "Bad satellite creation: " + type
 
    factory = staticmethod(factory)

    common_data = {}

    common_data['baseline_cmask_goes_nop_cloud_mask'] = {
            'name':'Cloud Mask',
            'quantity': 'cloud mask',
            'discrete':True,
            'values':(0,1,2,3),
            'mask_values':[],
            'units': None,
            'fill_boundaries':[-0.5,0.5,1.5,2.5,3.5],
            'fill_colours':['w','#ff0000','#00ffff','#00ff00'],
            'tick_names':['confident cloudy','probably cloudy','probably clear','confident clear'],
            'cmap':None,
            'meridian_color':'black',
            'coastline_color':'black',
            'country_color':'black'
            }
    common_data['baseline_cmask_seviri_cloud_mask'] = common_data['baseline_cmask_goes_nop_cloud_mask']


    common_data['eps_cmask_ahi_cloud_mask'] = {
            'name':'Cloud Mask',
            'quantity': 'cloud mask',
            'discrete':True,
            'values':(0,1,2,3),
            'mask_values':[],
            'units': None,
            'fill_boundaries':[-0.5,0.5,1.5,2.5,3.5],
            'fill_colours':['w','#ff0000','#00ffff','#00ff00'],
            'tick_names':['confident cloudy','probably cloudy','probably clear','confident clear'],
            'cmap':None,
            'meridian_color':'black',
            'coastline_color':'black',
            'country_color':'black'
            }

    common_data['eps_cmask_ahi_cld_probability'] = {
                'name':'Cloud Probability',
                'quantity':'cloud probability',
                'discrete':False,
                'values':[0.,1.],
                'mask_ranges':[],
                'logscale':False,
                'units': None,
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                #'cmap':cm.gray,
                'cmap':g2_cmaps.magma,
                #'cmap':g2_cmaps.inferno,
                #'cmap':g2_cmaps.plasma,
                #'cmap':g2_cmaps.viridis,
                'n_levels':256
                }

    common_data['goesnp_ctype_cloud_phase'] = {
            'name':'Cloud Phase',
            'quantity':'cloud phase',
            'discrete':True,
            'values':(0,1,2,3,4,5),
            'mask_values':[],
            'units': None,
            'fill_boundaries':[-0.5,  0.5,  1.5,  2.5,  3.5,  4.5,  5.5],
            'fill_colours':['#000000','#02BCFC','#01FB00','#006701','#FAF900','#F600FD'],
            'tick_names':['clear','water','supercooled','mixed','ice','unknown'],
            'cmap':None,
            'coastline_color':'cyan',
            'country_color':'magenta',
            'meridian_color':'yellow'
            }
    common_data['enterprise_cldphase_10_11_13_14_15_cloud_phase'] = common_data['goesnp_ctype_cloud_phase']

    common_data['goesnp_ctype_cloud_type'] = {
            'name':'Cloud Type',
            'quantity':'cloud type',
            'discrete':True,
            'values':(0,1,2,3,4,5,6,7,8,9),
            'mask_values':[],
            'units': None,
            'fill_boundaries':[-0.5,  0.5,  1.5,  2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5,  9.5],
            'fill_colours':['#000000','#0002F5','#02BCFC','#01FB00','#006701','#FAF900','#F90101','#FE8802','#8F9490','#F600FD'],
            'tick_names':['clear','spare','water','SC','Mixed','Thick Ice','Thin Ice','Multilay','Spare','Uncertain'],
            'cmap':None,
            'coastline_color':'cyan',
            'country_color':'magenta',
            'meridian_color':'yellow'
            }
    common_data['enterprise_cldphase_10_11_13_14_15_cloud_type'] = common_data['goesnp_ctype_cloud_type']

    common_data['ACHA_mode_7_goes_cloud_emissivity'] = {
                'name':'Cloud Emissivity',
                'quantity':'emissivity',
                'discrete':False,
                'values':[0.,1.],
                'mask_ranges':[],
                'logscale':False,
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }
    common_data['eps_cmask_ahi_emiss11_high'] = common_data['ACHA_mode_7_goes_cloud_emissivity']

    common_data['ACHA_mode_7_goes_cloud_optical_depth_vis'] = {
                'name':'Cloud Optical Depth (visible)',
                'quantity':'optical depth',
                'discrete':False,
                'values':[0.,12.],
                'mask_ranges':[],
                'logscale':False,
                'units': None,
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }

    common_data['ACHA_mode_7_goes_cloud_particle_effective_radius'] = {
                'name':'Cloud Particle Effective Radius',
                'quantity':'effective radius',
                'discrete':False,
                'values':[0.,40.],
                'mask_ranges':[],
                'logscale':False,
                'units': '$\mu\mathrm{m}$',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }

    common_data['ACHA_mode_7_goes_cloud_top_height'] = {
                'name':'Cloud Top Height',
                'quantity':'height',
                'discrete':False,
                'values':[0.,18000.],
                'mask_ranges':[],
                'logscale':False,
                'units': '$\mathrm{m}$',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }
    common_data['ACHA_mode_8_cloud_top_height'] = common_data['ACHA_mode_7_goes_cloud_top_height']

    common_data['ACHA_mode_7_goes_cloud_top_pressure'] = {
                'name':'Cloud Top Pressure',
                'quantity':'pressure',
                'discrete':False,
                'values':[100.,1100.],
                'mask_ranges':[],
                'logscale':True,
                'units': '$\mathrm{hPa}$',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis_r,
                'n_levels':256
                }
    common_data['ACHA_mode_8_cloud_top_pressure'] = common_data['ACHA_mode_7_goes_cloud_top_pressure']

    common_data['ACHA_mode_7_goes_cloud_top_temperature'] = {
                'name':'Cloud Top Temperature',
                'quantity':'temperature',
                'discrete':False,
                'values':[180.,295.],
                'mask_ranges':[],
                'logscale':False,
                'units': '$\mathrm{K}$',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis_r,
                'n_levels':256
                }
    common_data['ACHA_mode_8_cloud_top_temperature'] = common_data['ACHA_mode_7_goes_cloud_top_temperature']

    common_data['DCOMP_mode_3_cloud_albedo'] = {
                'name':'Cloud Albedo',
                'quantity':'albedo',
                'discrete':False,
                'values':[0.,1.],
                'mask_ranges':[],
                'logscale':False,
                'units': None,
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }

    common_data['DCOMP_mode_3_cloud_spherical_albedo'] = {
                'name':'Cloud Spherical Albedo',
                'quantity':'albedo',
                'discrete':False,
                'values':[0.,1.],
                'mask_ranges':[],
                'logscale':False,
                'units': None,
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }

    common_data['DCOMP_mode_3_cloud_ice_water_path'] = {
                'name':'Cloud Ice Water Path',
                'quantity':'ice water path',
                'discrete':False,
                'values':[0., 1000.],
                'mask_ranges':[],
                'logscale':False,
                'units': '$\mathrm{g m}^{-2}$',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }

    common_data['DCOMP_mode_3_cloud_liquid_water_path'] = {
                'name':'Cloud Liquid Water Path',
                'quantity':'liquid water path',
                'discrete':False,
                'values':[0., 1000.],
                'mask_ranges':[],
                'logscale':False,
                'units': '$\mathrm{g m}^{-2}$',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }

    common_data['DCOMP_mode_3_cloud_optical_depth_vis'] = {
                'name':'Cloud Optical Depth (visible)',
                'quantity':'optical depth',
                'discrete':False,
                'values':[0.,30.],
                'mask_ranges':[],
                'logscale':False,
                'units': None,
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }

    common_data['DCOMP_mode_3_cloud_particle_effective_radius'] = {
                'name':'Cloud Particle Effective Radius',
                'quantity':'effective radius',
                'discrete':False,
                'values':[0.,40.],
                'mask_ranges':[],
                'logscale':False,
                'units': '$\mu\mathrm{m}$',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }

    common_data['goesr_fog_IFR_fog_probability'] = {
                'name':'IFR Fog Probability',
                'quantity':'probability',
                'discrete':False,
                'values':[0.,100.],
                'mask_ranges':[],
                'logscale':False,
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.geocat_colormaps.get_cmap_fog(cmap_break=0.37),
                'n_levels':256
                }
    common_data['goesr_fog_bridge_IFR_fog_probability'] = common_data['goesr_fog_IFR_fog_probability']

    common_data['goesr_fog_LIFR_fog_probability'] = {
                'name':'LIFR Fog Probability',
                'quantity':'probability',
                'discrete':False,
                'values':[0.,100.],
                'mask_ranges':[],
                'logscale':False,
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.geocat_colormaps.get_cmap_fog(cmap_break=0.30),
                'n_levels':256
                }
    common_data['goesr_fog_bridge_LIFR_fog_probability'] = common_data['goesr_fog_LIFR_fog_probability']

    common_data['goesr_fog_MVFR_fog_probability'] = {
                'name':'MVFR Fog Probability',
                'quantity':'probability',
                'discrete':False,
                'values':[0.,100.],
                'mask_ranges':[],
                'logscale':False,
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.geocat_colormaps.get_cmap_fog(cmap_break=0.55),
                'n_levels':256
                }
    common_data['goesr_fog_bridge_MVFR_fog_probability'] = common_data['goesr_fog_MVFR_fog_probability']

    common_data['goesr_fog_bridge_fog_depth'] = {
                'name':'Fog Depth',
                'quantity':'depth',
                'discrete':False,
                'values':[],
                'mask_ranges':[],
                'logscale':False,
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':g2_cmaps.geocat_colormaps.get_cmap_fog(cmap_break=0.55),
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
                'cmap':g2_cmaps.viridis,
                'n_levels':256
                }


class GOES_NOP(Satellite):

    data = {}
    data.update(Satellite.common_data)

class Himawari(Satellite):

    data = {}
    data.update(Satellite.common_data)




#class Dataset_Options:
    #"""
    #This class contains static data for the interpretation of the various  
    #geocat discrete and continuous datasets.
    #"""

    #data = {}

    #data['baseline_cmask_goes_nop_cloud_mask'] = {
            #'name':'Cloud Mask',
            #'quantity': 'cloud mask',
            #'discrete':True,
            #'values':(0,1,2,3),
            #'mask_values':[],
            #'units': None,
            #'fill_boundaries':[-0.5,0.5,1.5,2.5,3.5],
            #'fill_colours':['w','#ff0000','#00ffff','#00ff00'],
            #'tick_names':['confident cloudy','probably cloudy','probably clear','confident clear'],
            #'cmap':None,
            #'meridian_color':'black',
            #'coastline_color':'black',
            #'country_color':'black'
            #}
    #data['baseline_cmask_seviri_cloud_mask'] = data['baseline_cmask_goes_nop_cloud_mask']


    #data['goesnp_ctype_cloud_phase'] = {
            #'name':'Cloud Phase',
            #'quantity':'cloud phase',
            #'discrete':True,
            #'values':(0,1,2,3,4,5),
            #'mask_values':[],
            #'units': None,
            #'fill_boundaries':[-0.5,  0.5,  1.5,  2.5,  3.5,  4.5,  5.5],
            #'fill_colours':['#000000','#02BCFC','#01FB00','#006701','#FAF900','#F600FD'],
            #'tick_names':['clear','water','supercooled','mixed','ice','unknown'],
            #'cmap':None,
            #'coastline_color':'cyan',
            #'country_color':'magenta',
            #'meridian_color':'yellow'
            #}

    #data['goesnp_ctype_cloud_type'] = {
            #'name':'Cloud Type',
            #'quantity':'cloud type',
            #'discrete':True,
            #'values':(0,1,2,3,4,5,6,7,8,9),
            #'mask_values':[],
            #'units': None,
            #'fill_boundaries':[-0.5,  0.5,  1.5,  2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5,  9.5],
            #'fill_colours':['#000000','#0002F5','#02BCFC','#01FB00','#006701','#FAF900','#F90101','#FE8802','#8F9490','#F600FD'],
            #'tick_names':['clear','spare','water','SC','Mixed','Thick Ice','Thin Ice','Multilay','Spare','Uncertain'],
            #'cmap':None,
            #'coastline_color':'cyan',
            #'country_color':'magenta',
            #'meridian_color':'yellow'
            #}

    #data['ACHA_mode_7_goes_cloud_emissivity'] = {
                #'name':'ACHA mode 7 GOES Cloud Emissivity',
                #'quantity':'emissivity',
                #'discrete':False,
                #'values':[0.,1.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': '%',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}

    #data['ACHA_mode_7_goes_cloud_optical_depth_vis'] = {
                #'name':'ACHA mode 7 GOES Cloud Optical Depth (visible)',
                #'quantity':'optical depth',
                #'discrete':False,
                #'values':[0.,12.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': None,
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}

    #data['ACHA_mode_7_goes_cloud_particle_effective_radius'] = {
                #'name':'ACHA mode 7 GOES Cloud Particle Effective Radius',
                #'quantity':'effective radius',
                #'discrete':False,
                #'values':[0.,40.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': '$\mu\mathrm{m}$',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}

    #data['ACHA_mode_7_goes_cloud_top_height'] = {
                #'name':'ACHA mode 7 GOES Cloud Top Height',
                #'quantity':'height',
                #'discrete':False,
                #'values':[100.,10000.],
                #'mask_ranges':[],
                #'logscale':True,
                #'units': '$\mathrm{m}$',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}

    #data['ACHA_mode_7_goes_cloud_top_pressure'] = {
                #'name':'ACHA mode 7 GOES Cloud Top Pressure',
                #'quantity':'pressure',
                #'discrete':False,
                #'values':[0.,1100.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': '$\mathrm{hPa}$',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix_r,
                #'n_levels':256
                #}

    #data['ACHA_mode_7_goes_cloud_top_temperature'] = {
                #'name':'ACHA mode 7 GOES Cloud Top Temperature',
                #'quantity':'temperature',
                #'discrete':False,
                #'values':[180.,295.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': '$\mathrm{K}$',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix_r,
                #'n_levels':256
                #}

    #data['DCOMP_mode_3_cloud_albedo'] = {
                #'name':'DCOMP mode 3 Cloud Albedo',
                #'quantity':'albedo',
                #'discrete':False,
                #'values':[0.,1.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': None,
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}

    #data['DCOMP_mode_3_cloud_spherical_albedo'] = {
                #'name':'DCOMP mode 3 Cloud Spherical Albedo',
                #'quantity':'albedo',
                #'discrete':False,
                #'values':[0.,1.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': None,
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}

    #data['DCOMP_mode_3_cloud_ice_water_path'] = {
                #'name':'DCOMP mode 3 Cloud Ice Water Path',
                #'quantity':'ice water path',
                #'discrete':False,
                #'values':[0., 1000.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': '$\mathrm{g m}^{-2}$',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}

    #data['DCOMP_mode_3_cloud_liquid_water_path'] = {
                #'name':'DCOMP mode 3 Cloud Liquid Water Path',
                #'quantity':'liquid water path',
                #'discrete':False,
                #'values':[0., 1000.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': '$\mathrm{g m}^{-2}$',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}

    #data['DCOMP_mode_3_cloud_optical_depth_vis'] = {
                #'name':'DCOMP mode 3 Cloud Optical Depth (visible)',
                #'quantity':'optical depth',
                #'discrete':False,
                #'values':[0.,30.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': None,
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}

    #data['DCOMP_mode_3_cloud_particle_effective_radius'] = {
                #'name':'DCOMP mode 3 Cloud Particle Effective Radius',
                #'quantity':'effective radius',
                #'discrete':False,
                #'values':[0.,40.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': '$\mu\mathrm{m}$',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}

    #data['goesr_fog_IFR_fog_probability'] = {
                #'name':'goesr fog IFR Fog Probability',
                #'quantity':'probability',
                #'discrete':False,
                #'values':[0.,100.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': '%',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':g2_cmaps.geocat_colormaps.get_cmap_fog(cmap_break=0.37),
                #'n_levels':256
                #}

    #data['goesr_fog_LIFR_fog_probability'] = {
                #'name':'goesr fog LIFR Fog Probability',
                #'quantity':'probability',
                #'discrete':False,
                #'values':[0.,100.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': '%',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':g2_cmaps.geocat_colormaps.get_cmap_fog(cmap_break=0.30),
                #'n_levels':256
                #}

    #data['goesr_fog_MVFR_fog_probability'] = {
                #'name':'goesr fog MVFR Fog Probability',
                #'quantity':'probability',
                #'discrete':False,
                #'values':[0.,100.],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': '%',
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':g2_cmaps.geocat_colormaps.get_cmap_fog(cmap_break=0.55),
                #'n_levels':256
                #}

    #data['unknown'] = {
                #'name':None,
                #'quantity':None,
                #'discrete':False,
                #'values':[None,None],
                #'mask_ranges':[],
                #'logscale':False,
                #'units': None,
                #'coastline_color':'cyan',
                #'country_color':'magenta',
                #'meridian_color':'yellow',
                #'cmap':cm.cubehelix,
                #'n_levels':256
                #}
