#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 14:34:54 2024

@author: wilja48p
"""

##### SELECT 2022 NSHM RUPTURE SETS THAT INTERSECT AND WITHIN SEARCH POLYGON ######
### THEN RETURN A GEOJSON FILE FOR RUPTURES THAT INTERSECT THE POLYGON  #######

#A paticipation rate is also calculated for each section. For rupture sections within the 
#polygon, this should be equal to the full partipation rate of that section. For rupture sections
#outside the polygon, it is the rate for ruptures that intersect the polygon only

import nzshm_model as nm
import json
import pathlib
import geopandas 

from shapely.geometry import shape
from shapely.geometry.polygon import Polygon 

from shapely.ops import transform #added by JW
import pyproj #added by JW
import numpy as np #added by JW
import os #added by JW


from solvis import CompositeSolution, FaultSystemSolution, export_geojson

#Set logic tree exploration option
# 1 = weighted mean rates from all branches
# 2 = select certain branch

logic_tree_opt=1

#select branch. 
if logic_tree_opt==2:
    #uncomment to print branch options
    #ref_table=open("nshm_logictreebranch_lookuptable.txt",'r')
    #branch_ref_table=ref_table.read()
    #print(branch_ref_table) 
    branch="U2NhbGVkSW52ZXJzaW9uU29sdXRpb246MTIwNzI4"


# load the search polygon to select ruptures

os.chdir("..") #Go up one level in directory
os.chdir("..") #Go up one level in directory
geojson = json.load(open('gis_files/orb_area_polygon.geojson'))
polygon: Polygon = shape(geojson['features'][0]['geometry'])

# load the model
os.chdir('nshm_inversion')
slt = nm.get_model_version("NSHM_v1.0.4").source_logic_tree()
comp = CompositeSolution.from_archive("NSHM_v1.0.4_CompositeSolution.zip", slt)

#Go back to solvis directory
os.chdir('solvis')

# get the crustal fault solutions
fss: FaultSystemSolution = comp._solutions['CRU']  # NB this API call will change in the future

# get the rupture_ids that intersect or are within the selected polygon
rupture_ids: list = fss.get_ruptures_intersecting(polygon)
rupture_ids_list=rupture_ids.tolist()


#get rupture rates that are within search polygon based on logic tree opt

if logic_tree_opt==1: #Use weighted rupture rates from all logic tree branches
    rr = fss.ruptures_with_rupture_rates
    rup_info: list = rr[rr['Rupture Index'].isin(rupture_ids)]
    rup_info_list=rup_info['Rupture Index'].tolist()

elif logic_tree_opt==2: #Use rupture rates from specific logic tree branch
    cr = fss.composite_rates 
    #select rupture rates by logic tree branch
    rup_info=cr[(cr["Rupture Index"].isin(rupture_ids)) & (cr.solution_id==branch)]
    rup_info_list=rup_info['Rupture Index'].tolist()  
    #select only rupture ids that are in logic tree branch solution
    rupture_ids=rupture_ids[rupture_ids.isin(rup_info_list)]
    rupture_ids_list=rupture_ids.tolist()


#get fault sections of ruptures
#from line 285 at: https://github.com/GNS-Science/solvis-graphql-api/blob/181a42ee119f664a7fba16ed55c067073512935a/solvis_graphql_api/composite_solution/cached.py#L246
fsr = fss.fault_sections_with_rupture_rates

#get fault sections of intersecting ruptures
rup_fs = fsr[fsr['Rupture Index'].isin(rupture_ids)]
rup_fs_list=rup_fs['Rupture Index'].tolist()
rup_sec=rup_fs['section'].tolist() #all unique fault sections in rup_fs
unique_rup_sec=np.unique(rup_fs['section']) #all unique fault sections in rup_fs

tmp_indx1=[0]*len(unique_rup_sec)
sec_rate=[0]*len(unique_rup_sec) 
    
#loop through each unique section
for rid in range(unique_rup_sec.size):
    
    tmp_indx1[rid]=rup_sec.index(unique_rup_sec[rid])
    
    #for each section, find all the ruptures its associated with in rup_fs
    tmp_indx2 = np.where(np.isin(rup_fs['section'], unique_rup_sec[rid]))[0].tolist()
    rate=[0] * len(tmp_indx2) 
    
    for jj in range(len(tmp_indx2)):
        
        #for each section's rupture, index the weighted mean rupture rate
        if logic_tree_opt==1:
            rate[jj]=rup_fs['rate_weighted_mean'][tmp_indx2[jj]]
        
        #for each section's rupture, index the rupture in rup_info and then extract rate
        elif logic_tree_opt==2:        
            tmp_indx3=rup_info_list.index(rup_fs_list[tmp_indx1[jj]])
            rate[jj]=rup_info['Annual Rate'][tmp_indx3]
     
        #combine the section's rate for all ruptures it particpated in    
        sec_rate[rid]=sum(rate) 
 

#remove rows in rup_fs with duplicate sections
rup_fs_unique=rup_fs.iloc[np.array(tmp_indx1)]

#add section prticipation rates to dataframe
rup_fs_unique.insert(len(rup_fs_unique.columns),"Participation_Rate",sec_rate)

#convert section particpation rates to geojson file
gdf = geopandas.GeoDataFrame(rup_fs_unique, crs="EPSG:4326")

if logic_tree_opt ==1:
    export_geojson(gdf,"WORK/otago_sections_weighted.geojson") 
elif logic_tree_opt ==2:
    filename_list=["WORK/otago_sections_" + branch + ".geojson"]
    tmp_str=""
    filename=tmp_str.join(filename_list)
    export_geojson(gdf,filename)


