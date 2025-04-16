# NSHM INVERSION FAULT MODEL ANALYSIS WITH SOLVIS [RUNS IN PYTHON]

Scripts for analysing seismicity forecasts for Otago from the NZ NSHM 2022 Inversion Fault Model(Gerstenberger et al 2024)

Last updated by JW 16/04/25

## DEPENDENCIES 

### From the GNS Science NSHM project

-solvis library for geospatial analysis of the OpenSHA inversion solutions from the NZ NSHM 2022: https://github.com/GNS-Science/solvis/

-nzshm-model containing meta data about ther NSHM model logic tree branches: https://github.com/GNS-Science/nzshm-model

All project dependencies are listed in the file pyproject.toml

### Included in this folder

scripts/...: Folder containing python scripts for analysing the IFM results

OtagoFaults/... : Folder containing outputs to script runs

### Source data

-NSHM_v1.0.4_CompositeSolution.zip: inversion results from v1.0.4 of the NZ NSHM 2022. WARNING. Do not use the rupture rates in this file, as they been scaled by the non-stationary moment rate weighting, and for ruptures <Mw 8, by the 0.8 scaling for combination with the NZ NSHM 2022 DSM in fault polygons (see Gerstenberger et al 2024).
-"branch_id".zip where 'branch_id' corresponds to the following logic tree branches:
    "geol_b0.823": "SW52ZXJzaW9uU29sdXRpb246MTEzMDUz",
    "geol_b0.959": "SW52ZXJzaW9uU29sdXRpb246MTEzMDMy",
    "geol_b1.089": "SW52ZXJzaW9uU29sdXRpb246MTEzMDM5",
    "geod_b0.823": "SW52ZXJzaW9uU29sdXRpb246MTEzMDY3",
    "geod_b0.959": "SW52ZXJzaW9uU29sdXRpb246MTEzMDgw",
    "geod_b1.089": "SW52ZXJzaW9uU29sdXRpb246MTEzMDYz",


- All source files should be located in the folder WORK/.... To acquire a copy of these files, please send an email to: nshm@gnz.cri.nz

## INSTILLATION

This project uses poetry to install and manage dependencies. To install poetry, see https://python-poetry.org/docs/#installing-with-the-official-installer

The project requires python 3.10 or above. Use your preferred python environment manager to provide a clean python3.x , then run poetry install which will install the required libraries.


## SCRIPTS

Scripts select_ruptures_by_polygon.py, get_rupture_traces.py, and get_section_area.py need to be run in solvis directory, and through command line using poetry: poetry run python scripts/xxxxxx.py 

Files need to be run in this order:

1.) nshm_inversion/solvis/scripts/select_ruptures_by_polygon.py: Select ruptures from a given logic tree branch's inversion solution that are within or intersect a spatial polygon (here set to 'orb_area_polygon'). Returns .csv file that collates the rupture info. Calculates a down-weighted rupture magnitude based on how many of its sections are within the spatial polygon. See "MODFIY THESE LINES" section for where to select logic tree branch

2.) nshm_inversion/solvis/scripts/get_section_area.py: Returns the area of each fault section within a spatial polygon: 'otago_section_area.csv'. This is necessary for finding a fault's partial moment within a given rupture (as used in 'nshm_inversion_results_by_fault.m').
Also returns a list of Otago ruptures ids from ALL logic tree branch solutions (i.e. from NSHM_v1.0.4_CompositeSolution). This is necessary to run combine_branches.py

3.) nshm_inversion/solvis/scripts/combine_branches.py: Collates the rupture rates from the 3 geologic and 3 geodetic logic tree branch solutions, and returns the weighted average rate for each rupture in files 'geo_rates_w.csv' and 'ged_rates_w.csv.' This is necessary for running get_rupture_traces.py and the scripts in nshm_inversion/mfd_analysis. This can be run through a Python IDE (e.g., Spyder). 

4.) nshm_inversion/solvis/scripts/get_rupture_traces.py: Select all fault sections of ruptures that are within or intersect a spatial polygon (here set to orb_area_polygon) for a given, logic tree branch, or from the 3 geologic and geodetic branches analysed in this study. Then for each fault section, sum the rates for all ruptures that are associated with the section (i.e. 'participation rate'), and return this with the section info as a .geojson file. These options are selected in the coded bracketed by "MODIFY THESE LINES." These geojson files can be used to plot maps of fault sections within traces coloured by rate (equivalent to Figure S3 in manuscript). Also option to create geojson files for each rupture in the inversion solution from the 3 geologic logic tree branches considered in this study (only if opt=2).

5.) nshm_inversion/solvis/scripts/get_fault_rupture_patches.py: For faults listed in 'fault_select', will find all ruptures that they participate in and their rates for the 3 geologic deformation logic tree branches explored here. Will then plot the fault's rupture patches in figure , coloured semi-quantitatively by the rupture's rate. Plot is equivalent to Figure S10 in the manuscript. Requires that get_rupture_traces.py with opt=2 is run first so that rupture surfaces have been created in folder 'OtagoFaults/rupture_surfaces' and 'fault_sections_geo_w.csv' lookup table is created. This script can be run through a Python IDE (e.g., Spyder). 

## REFERENCES

Gerstenberger, M. C., Van Dissen, R., Rollins, C., DiCaprio, C., Thingbaijim, K. K., Bora, S., ... & Williams, C. (2024). The seismicity rate model for the 2022 Aotearoa New Zealand national seismic hazard model. Bulletin of the Seismological Society of America, 114(1), 182-216.

## TROUBLE SHOOTING

-If solvis cannot import modules, in command line, and in the correct environment and working directory, run: poetry install [https://github.com/python-poetry/poetry/issues/1868] 

-If return errors such as "EOFError: marshal data too short" ensure ALL __pycache__ folders are empty and/or check directory to solvis that is included in them 
