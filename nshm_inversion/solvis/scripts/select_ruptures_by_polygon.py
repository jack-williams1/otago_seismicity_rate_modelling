##### SELECT 2022 NSHM RUPTURE SETS THAT INTERSECT AND WITHIN SEARCH POLYGON ######
### THEN RETURN THE RUPTURES'S MOMENT MAGNITUDE THAT IS RELEASED WITHIN THE POLOGON ###

import json
import pathlib

from shapely.geometry import shape
from shapely.geometry.polygon import Polygon 

from shapely.ops import transform #added by JW
import pyproj #added by JW
import numpy as np #added by JW
import os

from solvis import InversionSolution
from solvis.filter import FilterRuptureIds

def section_intersecting(soln: InversionSolution, rupture_id: int, rupture_ids_list, rup_info_list):

    
    #get index of intersecting rupture
    #Rupture Index in rupture_ids and rup_info are NOT in the same order
    #so need seperate indicies
    tmp_idx1=rupture_ids_list.index(rupture_id)
    tmp_idx2=rup_info_list.index(rupture_id)
    
    #get surfaces of each rupture section
    #rupture_surface_gdf=fss.rupture_surface(rupture_ids[tmp_idx1])
    rupture_surface_gdf=soln.rupture_surface(rupture_ids_list[tmp_idx1])
    rupture_sections_info=rupture_surface_gdf[['key_0','section','geometry','DipDeg','Magnitude','Area (m^2)']]
    rupture_polygon=rupture_sections_info['geometry'].tolist()
    section_dip=rupture_sections_info['DipDeg'].tolist()
    
    contained_indx=[0] * len(rupture_polygon) 
    weighted_section_area=[0] * len(rupture_polygon)
    
    #index the ruture sections individual geometry and area
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
        
        #return value of 1 if rupture section entirely within search polygon
        if polygon.contains(section_polygon):
            contained_indx[ss]=1
        else:
            #return value of 0.5 if rupture section partly within polygon
            if section_polygon.intersects(polygon): 
                contained_indx[ss]=0.5
            
            #return value of 0 if rupture section not within polygon
            else:
                contained_indx[ss]=0
        
        #return section area weighted by whether it crosses polygon
        weighted_section_area[ss]=contained_indx[ss]*section_area
    
    #if rupture only partly within polygon, derive new area and magnitude
    if  min(contained_indx)<1:   
        weighted_rupture_area[tmp_idx2]=sum(weighted_section_area) #Total area of rupture within polygon
        rupture_area=(rupture_sections_info['Area (m^2)'].tolist())[0] #Total area of rupture
        
        rupture_mw[tmp_idx2]=(rupture_sections_info['Magnitude'].tolist())[0]#get total magnitude of rupture
        rupture_mo=10**(rupture_mw[tmp_idx2]*1.5+9.05) #Total seismic moment of rupture
        weighted_rupture_mo=rupture_mo*weighted_rupture_area[tmp_idx2]/rupture_area #seismic moment released within Otago polygon
        weighted_rupture_mw[tmp_idx2]=(np.log10(weighted_rupture_mo)-9.05)/1.5 #moment magnitude of polygon released within Otago
    
    #if rupture entirely within polygon no need to change area and magnitude
    else: 
        rupture_mw[tmp_idx2]=(rupture_sections_info['Magnitude'].tolist())[0]
        weighted_rupture_area[tmp_idx2]=rupture_sections_info['Area (m^2)'].tolist()[0]    
        weighted_rupture_mw[tmp_idx2]=rupture_mw[tmp_idx2]
    
    return weighted_rupture_area, rupture_mw, weighted_rupture_mw 

logic_tree_branches = {
    "geol_b0.823": "SW52ZXJzaW9uU29sdXRpb246MTEzMDUz",
    "geol_b0.959": "SW52ZXJzaW9uU29sdXRpb246MTEzMDMy",
    "geol_b1.089": "SW52ZXJzaW9uU29sdXRpb246MTEzMDM5",
    "geod_b0.823": "SW52ZXJzaW9uU29sdXRpb246MTEzMDY3",
    "geod_b0.959": "SW52ZXJzaW9uU29sdXRpb246MTEzMDgw",
    "geod_b1.089": "SW52ZXJzaW9uU29sdXRpb246MTEzMDYz",
}



# ================ MODIFY THESE LINES ================ 
# set your path
work_dir = pathlib.Path("WORK/")
output_dir= pathlib.Path("OtagoFaults/")
# select branch
branch = "geod_b0.823"
# ================ MODIFY THESE LINES ================ 

soln_id = logic_tree_branches[branch]

# load the search polygon to select ruptures

os.chdir("..") #Go up one level in directory
os.chdir("..") #Go up one level in directory
geojson = json.load(open('gis_files/orb_area_polygon.geojson'))
polygon: Polygon = shape(geojson['features'][0]['geometry'])

os.chdir("nshm_inversion/solvis") 

# load the inversion solution
soln_filepath = work_dir / f"{soln_id}.zip"
soln = InversionSolution.from_archive(soln_filepath)


# get the rupture_ids that intersect and within the selected polygon
rupture_ids = FilterRuptureIds(soln).for_polygon(polygon)
rupture_ids_list=rupture_ids.tolist()

# ruptures = soln.model.ruptures_with_rupture_rates
rr = soln.model.ruptures_with_rupture_rates
rup_info: list = rr[rr['Rupture Index'].isin(rupture_ids)]
rup_info_list=rup_info['Rupture Index'].tolist()

#set dummy variables 
rupture_mw=[0] * len(rupture_ids_list) 
weighted_rupture_area=[0] * len(rupture_ids_list) 
weighted_rupture_mw=[0] * len(rupture_ids_list) 

#for rid in rupture_ids[:len(rupture_ids)]: 
for rid in rupture_ids_list[:len(rupture_ids_list)]:    
    contain=section_intersecting(soln, rid, rupture_ids_list, rup_info_list)

#Reshape rows in contain to columns and add to new rup_info_w dataframe
rup_info_w=rup_info.assign(weighted_area=np.reshape(contain[0:1],(len(rupture_ids),1)))
rup_info_w=rup_info_w.assign(magnitude=np.reshape(contain[1:2],(len(rupture_ids),1)))
rup_info_w=rup_info_w.assign(weighted_magnitude=np.reshape(contain[2:],(len(rupture_ids),1)))

ruptures_csv_filename = f"{soln_id}_Mmod_Rtrim.csv"
ruptures_csv_filepath = output_dir / ruptures_csv_filename
rup_info_w.to_csv(ruptures_csv_filepath)
    
