## COMPILE EVENT TIMINGS FROM OTAGO'S PALEOSEISMIC RECORD AND USE TO EXTRACT EVENT RATE ## [RUNS IN MATLAB]


## INPUTS ##

paleoseismic_constraints/OtagoEarthquakeTimings.xlsx: Spreadsheet with event timings for all known surface ruptures within the Otago study area. See 'Refs' column and Table 3 in manuscript for references for these estimates.
paleoseismic_constraints/assymmetric_dist/: Folder containing text files that describe the event timings where the event timing's probability distribution is assymmetric (hence ages cannot be treated as a normal distribution

## SCRIPTS ##

paleoseismic_constraint/otago_eq_timings_simulations.m: For each rupture timing's probability distribution, randomly sample 10,000 earthquake ages. Then separately, randomly generate 10,000x 10,000-year random time windows between 0-20 ka. For each time window, derive how many earthquakes occurred within it from one of the 10,000 randomly sampled earthquake timings. Corrections made for the equivocal history of the Dunstan Fault's paleoseismic record. From the total number of events in a 10,000-year time record, an event rate can be generated and compared to the SRMs. Stored as 'otago_paleoseismic_constraint.'
paleoseismic_constraint/slip_rate_timeframe.m: For each Otago fault in the study area, plot a histogram to show the distribution of time intervals represented by their slip rate estimate in the NZ CFM (i.e., the SRT_pref or SRT_gen parameters in the NZ CFM; Seebeck et al 2024). Used to plot Figure S1 in the manuscript

## REFERENCES ## 

-Seebeck, H., Van Dissen, R., Litchfield, N., Barnes, P. M., Nicol, A., Langridge, R., ... & Lee, J. (2024). The New Zealand Community Fault Model–version 1.0: An improved geological foundation for seismic hazard modelling. New Zealand Journal of Geology and Geophysics, 67(2), 209-229.