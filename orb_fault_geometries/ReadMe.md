## EXTRACT OTAGO FAULTS FROM NZ CFM GEOMETRIC MODEL ## [RUNS IN MATLAB]

!First run function 'ReadAndSaveCfm.m' in MATLAB command window (Meade 2004). This will convert all individual .ts files  (Tsurf) for Otago faults in this folder from the NZ CFM (https://www.gns.cri.nz/data-and-resources/new-zealand-community-fault-model/; Seebeck et al 2022,2024) and convert them to a .mat variable 'Cfm1.2Pre_V7.mat'(note this has been done already)

orb_fault_geometries/faultgeometries.m: Loads 'Cfm1.2Pre_V7.mat' and for each fault, extracts its area as given in the NZ CFM. This is then written into 'OtagoRangeBasinFaults.xlsx' spreadsheet, so that the areas can be used for calculating fault's moment rate in the stochastic event catalogs (this has already been done). Note, the area of each fault is multiplied by 0.8. This is for consistency with how the bottom fault depth (Dfc) was handled in the NSHM Geologic Deformation model (van Dissen et al 2024).

## REFERENCES ##

-Meade (2004): script writer for the SCEC CFM to convert.ts files to MATLAB variable. Available at: https://github.com/kaeufl/ReadAndSaveCfm/blob/master/ReadAndSaveCfm/ReadAndSaveCfm.m
-Seebeck, H., R. J. Van Dissen, N. J. Litchfield, P. M. Barnes, A. Nicol, R. M. Langridge, D. J. A. Barrell, P. Villamor, S. M. Ellis, M. S. Rattenbury, et al. (2022). New Zealand Community Fault Model–version 1.0, GNS Science report 2021/57, GNS Science, Lower Hutt, New Zealand, 97 pp., doi: 10.21420/GA7S-BS61.
-Seebeck, H., Van Dissen, R., Litchfield, N., Barnes, P. M., Nicol, A., Langridge, R., ... & Lee, J. (2024). The New Zealand Community Fault Model–version 1.0: An improved geological foundation for seismic hazard modelling. New Zealand Journal of Geology and Geophysics, 67(2), 209-229.
-Van Dissen, R. J., Johnson, K. M., Seebeck, H., Wallace, L. M., Rollins, C., Maurer, J., ... & DiCaprio, C. J. (2024). Upper Plate and Subduction Interface Deformation Models in the 2022 Revision of the Aotearoa New Zealand National Seismic Hazard Model. Bulletin of the Seismological Society of America, 114(1), 37-56.