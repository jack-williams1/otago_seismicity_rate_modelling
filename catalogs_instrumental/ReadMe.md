## ANALYSIS OF INSTRUMENTAL CATALOGS [RUNS IN MATLAB]

Analysis of Otago events in the Augmented New Zealand earthquake catalog (Rollins et al 2022).

## DEPENDENCIES ##

Download the Augmented New Zealand earthquake catalog from: https://geodata.nz/geonetwork/srv/api/records/00876699-4791-4528-a826-1d7ac3c77516 AND place in this folder under the name 'NZNSHM2022_augmentedEQcatalogue_annotated_2Aug2022.txt'

## SCRIPTS ##

catalogs_instrumental/nz_augmented_catalog.m: select events from the New Zealand Augmented Earthquake Catalog that are located within the Otago search polygon ('gis_files/orb_area_polygon.shp'). Selects: (1) all events, (2) events 1951-2021, and (3) events 1971-2021. Options to plot MFD, models catalog G-R relationship using maximum likelihood function: GRrelation_MLEWeichert_EQMATca.m (Tinti and Mulargia 1987), and return sampled catalog as .mat variable and csv file. These output files are needed for subsequent analysis in all_catalog_comparison.m

## REFERENCES ##

-Rollins C, Thingbaijam KK, Hutchinson, J, Gerstenberger M, Christophersen A, Eberhart-Phillips D, Rastin SJ, Van Dissen R. (2022). An augmented New Zealand earthquake catalogue, event classifications, and models of the depth distribution of shallow earthquakes in the greater New Zealand region Lower Hutt (NZ): GNS Science. 83 p. (GNS Science report; 2021/58).

-Tinti, S., & Mulargia, F. (1987). Confidence intervals of b values for grouped magnitudes. Bulletin of the Seismological Society of America, 77(6), 2125-2134.
