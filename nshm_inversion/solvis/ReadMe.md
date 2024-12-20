# NSHM INVERSION FAULT MODEL ANALYSIS WITH SOLVIS [RUNS IN PYTHON]

Scripts for analysing seismicity forecasts for Otago from the NZ NSHM 2022 Inversion Fault Model(Gerstenberger et al 2024)

## DEPENDENCIES ##

From the GNS Science NSHM project:

 - [solvis](https://github.com/GNS-Science/solvis/) library for geospatial analysis of the OpenSHA inversion solutions from the NZ NSHM 2022.
 - [nzshm-model](https://github.com/GNS-Science/nzshm-model) containing meta data about ther NSHM model logic tree branches.

All project dependencies are listed in the file [pyproject.toml](pyproject.toml)

### Included in this folder:

 - [nshm_logictreebranch_lookuptable.txt](nshm_logictreebranch_lookuptable.txt): Lookup table for parameters that each inversion fault model logic tree branch explores. NB this information is also available in the nzshm-model package.

 - scripts/...: Folder containing python scripts for analysing the IFM results

 - WORK/... : Folder containing outputs to script runs

- included in parent directory (nshm_inversion):


## Source data

The file **NSHM_v1.0.4_CompositeSolution.zip** should be located in the parent folder `nshm_inversion`. This file contains the  inversion results from v1.0.4 of the NZ NSHM 2022. To acquire a copy, please send an email to: nshm@gnz.cri.nz

## INSTALLATION

This project uses **poetry** to install and manage dependencies. To install poetry, see https://python-poetry.org/docs/#installing-with-the-official-installer

The project requires python 3.10 or above. Use your preferred python environment manager to provide a clean python3.x , then run `poetry install` which will install the required libraries.

## SCRIPTS

The scripts `select_ruptures_by_polygon.py`, `get_rupture_traces.py`, and `get_section_area.py` need to be run in solvis directory, and through command line using poetry e.g.
`poetry run python scripts/xxxxxx.py`

- `nshm_inversion/solvis/scripts/select_ruptures_by_polygon.py`: Select ruptures from the inversion solution that are within or intersect a spatial polygon (here set to 'orb_area_polygon'). Returns .csv file that collates the rupture info. Calculates a down-weighted rupture magnitude based on how many of its sections are within the spatial polygon. Analysis can be done for weighted average of all inversion logic tree branches or individual branches. See 'nshm_logictreebranch_lookuptable.txt' for info on what each branch means. Can also return 'otago_fault_sections.csv' which collates information about the individual ruptures that each section participates in (this is currently commented out). This is necessary for running scripts in nshm_inversion/mfd_analysis

- `nshm_inversion/solvis/scripts/get_section_area.py`: Returns the area of each fault section within a spatial polygon. This is necessary for finding a fault's partial moment within a given rupture (as used in 'nshm_inversion_results_by_fault.m'). Also option to return .geojson file for each rupture that an Otago fault participates in. This is necessary to run get_fault_rupture_patches.py

- `nshm_inversion/solvis/scripts/get_rupture_traces.py`: Select all fault sections of ruptures that are within or intersect a spatial polygon (here set to orb_area_polygon) for all, or a given, logic tree branch. Then for each fault section, sum the rates for all ruptures that are associated with the section (i.e. 'participation rate'), and return this with the section info as a .geojson file. Can be used to plot maps of fault sections within traces coloured by rate (equivalent to Figure S3 in manuscript).

- `nshm_inversion/solvis/scripts/get_fault_rupture_patches.py`: For some given faults, will find all ruptures that they participate in and their rates for the 3 geologic deformation logic tree branches explored here. Will then plot the fault's rupture patches in map, coloured semi-quantitatively by the rupture's rate. Plot is equivalent to Figure S10 in the manuscript. Requires that rupture surfaces have been created in get_section_area.py, and rupture rates are collated from nshm_otago_inversion_results_by_fault.m Needs to be run from the working directory '...nshm_inversion/solvis. This can be run through a Python IDE (e.g., Spyder).

## REFERENCES

Gerstenberger, M. C., Van Dissen, R., Rollins, C., DiCaprio, C., Thingbaijim, K. K., Bora, S., ... & Williams, C. (2024). The seismicity rate model for the 2022 Aotearoa New Zealand national seismic hazard model. Bulletin of the Seismological Society of America, 114(1), 182-216.

## TROUBLE SHOOTING

- If solvis cannot import modules, in command line, and in the correct environment and working directory, run: poetry install [https://github.com/python-poetry/poetry/issues/1868]
- If return errors such as "EOFError: marshal data too short" ensure ALL __pycache__ folders are empty and/or check directory to solvis that is included in them
