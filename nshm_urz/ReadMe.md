## EXTRACT GEODETIC SLIP RATES OF OTAGO FAULTS AND COMPARES TO GEOLOGIC [RUNS IN MATLAB]

## PREAMBLE

Analyse earthquake forecasts developed from the NZ NSHM 2022 Distributed Seismicity Model (DSM). The DSM has two components: Multiplicative hybrid forecast (Rastin et al 2024), and Uniform Rate Zones (Iturrieta et al 2022, 2024). These have then been combined using a floor ensemble method (Iturrieta et al 2024). These codes provide the ability to extract and analyse forecast rates for spatial cells contained within the Otago study area

## INPUTS

The following DSM forecasts are included here, which were downloaded from: https://github.com/pabloitu/nz_nshm2022_nonpoisson/tree/main:

-nshm_urz/inputs/fe_m.csv: Poisson floor ensemble forecast

-nshm_urz/inputs/m.csv: Hybrid model forecast rate

-nshm_urz/inputs/npfe_ao.csv: Negative binomial forecast- additive optimized

-nshm_urz/inputs/npfe_m.csv: Negative binomial forecast- multiplicative


Shape files for the different URZ spatial bins available at: https://github.com/pabloitu/nonpoisson_nz_src/tree/main/results/spatial:

nshm_urz/inputs/hw_final_j2_3... : Discretization of J2 strain rate tensor in NZ into 3 spatial bins

nshm_urz/inputshw_final_j2_2..: Discretization of J2 strain rate tensor in NZ into 2 spatial bins

nshm_urz/inputs/deg2utm.m: function to convert geographic coordinates (Palacios 2023)

## SCRIPTS

nshm_urz/find_centre_point.m: code to extract the centre point of each forecasts' spatial cell (where each spatial cell is a 0.1 x 0.1 long-lat box). The centre point and its forecast rate are then written into a .txt file, which is useful for plotting the rates as a raster file in (for example) QGIS. See for example, Figure 9a-c in the manuscript. [note this file has been run already]

nshm_urz/nshm_urz_analysis.m: First selects all spatial cells centred within Otago study area. The for each forecast, extracts rate, with option to plot rates, and save values in a .txt file. By summing the rates from all Otago cells, it is then possible to get the total regional rate (once renormalised for absolute rates of M>4.95 earthquakes in NZ hybrid volume. Using the crustal b-value then possible to construct MFD's for forecast, compare to last 50 years of Otago seismicity (using output of catalogs_instrumental/nz_augmented_catalog.m), and randomly simulate 1000x alternative forecasts for each recurrence mode;, as shown in Figure 9 in the main manuscript


## REFERENCES

-Iturrieta, P., Gerstenberger, M. C., Rollins, C., Van Dissen, R., Wang, T., & Schorlemmer, D. (2024). Accounting for the variability of earthquake rates within low‐seismicity regions: Application to the 2022 Aotearoa New Zealand National Seismic Hazard Model. Bulletin of the Seismological Society of America, 114(1), 217-243.

-Iturrieta, P., M. Gerstenberger, C. Rollins, R. J. Van Dissen, T. Wang, and D. Schorlemmer. "Accounting for earthquake rates’ temporal and spatial variability through least-information uniform rate zone forecasts." GNS Science Rept 14 (2022).

-Palacios, R. (2023). deg2utm (https://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm), MATLAB Central File Exchange. 

-Rastin, S. J., Rhoades, D. A., Rollins, C., Gerstenberger, M. C., Christophersen, A., & Thingbaijam, K. K. (2024). Spatial Distribution of Earthquake Occurrence for the New Zealand National Seismic Hazard Model 2022. Bulletin of the Seismological Society of America, 114(5), 2767-2788.
