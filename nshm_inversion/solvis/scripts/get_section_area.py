#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 10:12:19 2024

WRITTEN TO RETURN THE AREA OF EACH NSHM SECTION WITHIN OTAGO POLYGON
RETURNS AREA BASED ON RUPTURE SURFACE GEOMETRY PROJECTED THROUGH DIP
ALSO RETURNS LIST OF RUPTURE IDS IN OTAGO FROM ALL INVERSION SOLUTIONS:
"filtered_rupture_ids.csv"

@author: wilja48p
"""

import nzshm_model as nm
import json
import pathlib

from shapely.geometry import shape
from shapely.geometry.polygon import Polygon 

from shapely.ops import transform #added by JW
import pyproj #added by JW
import numpy as np #added by JW

import csv
import os

from solvis import CompositeSolution, FaultSystemSolution
from solvis.filter import FilterRuptureIds

# load the search polygon to select ruptures
# load the search polygon to select ruptures

os.chdir("..") #Go up one level in directory
os.chdir("..") #Go up one level in directory
geojson = json.load(open('gis_files/orb_area_polygon.geojson'))
polygon: Polygon = shape(geojson['features'][0]['geometry'])

os.chdir("nshm_inversion/solvis") 

# load the model
# Can still use unscaled rate full version of the model, as not interested in rates
slt = nm.get_model_version("NSHM_v1.0.4").source_logic_tree
comp = CompositeSolution.from_archive("WORK/NSHM_v1.0.4_CompositeSolution.zip", slt)

# get the crustal fault solutions
fss: FaultSystemSolution = comp._solutions['CRU']  # NB this API call will change in the future

# get the rupture_ids that intersect and within the selected polygon
# This includes all logic tree branches
rupture_ids = FilterRuptureIds(fss).for_polygon(polygon)
rupture_ids_list=rupture_ids.tolist()

#Create csv file that lists all rupture ids from complete inversion solution, needed to run combine_branches.py
with open('OtagoFaults/filtered_rupture_ids.csv', 'w', newline='') as file:

    writer = csv.writer(file)
    for ii in rupture_ids_list:
        writer.writerow([ii])
        
fsr = fss.model.fault_sections_with_rupture_rates
#get fault sections of intersecting ruptures
rup_fs = fsr[fsr['Rupture Index'].isin(rupture_ids)]

sec_area_np=np.zeros((len(rup_fs),2))

count=0

#loop through each Otago rupture
for rr in range(len(rupture_ids_list)):
#for rr in range(10):       
    rupture_surface_gdf=fss.rupture_surface(rupture_ids_list[rr])
    rupture_sections_info=rupture_surface_gdf[['key_0','section','geometry','DipDeg']]
    rupture_polygon=rupture_sections_info['geometry'].tolist()
    section_dip=rupture_sections_info['DipDeg'].tolist()

    
    #for each section within rupture
    for ss in range(len(rupture_polygon)):
        
        section_polygon=rupture_polygon[ss] #define the geometry of each rupture section
        
        # define a projection from lat/lon to NZTM
        projected_crs = pyproj.CRS.from_epsg(2193)
        projector = pyproj.Transformer.from_crs('epsg:4326', projected_crs, always_xy=True).transform
        
        #project lat/lon coordinates of polygon to NZTM
        projected_polygon = transform(projector, section_polygon)
        
        #derive area of polgon in m^2
        tmp_area=projected_polygon.area #area of section polygon with no info about dip
        section_area=tmp_area/np.cos(np.radians(section_dip[ss])) #section polygon area once projected through fault dip
        #sec_area_np[count,:]=[rupture_sections_info['section'][ss],section_area]
        sec_area_np[count,:]=[rupture_sections_info['section'][ss],section_area]
        count=count+1

#remove duplicate sections
sec_area_np_unique=np.unique(sec_area_np,axis=0)   

np.savetxt("OtagoFaults/otago_section_area.csv", sec_area_np_unique,  
              delimiter = ",")

