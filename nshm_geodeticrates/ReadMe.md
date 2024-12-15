## EXTRACT GEODETIC SLIP RATES OF OTAGO FAULTS AND COMPARES TO GEOLOGIC [RUNS IN PYTHON]

## DEPENDENCIES

GeoPandas: https://geopandas.org/en/stable/

Pandas: https://pandas.pydata.org

Shapely: https://pypi.org/project/shapely/

Matplotlib: https://matplotlib.org

Os: https://github.com/python/cpython/blob/3.13/Lib/os.py

This requires that "nshm_inversion/solvis/scripts/get_section_are.py" is run first so that can read 'otago_section_area.csv.'

## SCRIPTS

nshm_geodetic_rates/extract_otago_sections.py: extract slip rate estimates from the NZ NSHM 2022 geodetic model (Johnson et al 2022, 2024) for faults within the Otago study area polygon (orb_area_polygon.geojson). Data extracted from file 'fault_sections.geojson' and used to create file 'otago_sections_geodetic_rates.geojson' which is necessary for running the geodetic based stochastic catalogs. Code also used to compare geodetic and geologic (i.e., NZ CFM, Seebeck et al 2024, van Dissen et al 2024) slip rate estimates for Otago faults (Figure 12 in manuscript)

## REFERENCES

-Johnson, K., L. Wallace, J. Maurer, I. Hamling, C. Williams, C. Rollins, M. Gerstenberger, and R. Dissen (2022). Geodetic Deformation Model for the 2022 Update of the New Zealand National Seismic Hazard Model, GNS Science Report, GNS Science, Lower Hutt, New Zealand

-Johnson, K. M., Wallace, L. M., Maurer, J., Hamling, I., Williams, C., Rollins, C., ... & Van Dissen, R. (2024). Inverting geodetic strain rates for slip deficit rate in complex deforming zones: An application to the New Zealand plate boundary. Journal of Geophysical Research: Solid Earth, 129(3), e2023JB027565.
-Seebeck, H., Van Dissen, R., Litchfield, N., Barnes, P. M., Nicol, A., Langridge, R., ... & Lee, J. (2024). The New Zealand Community Fault Model–version 1.0: An improved geological foundation for seismic hazard modelling. New Zealand Journal of Geology and Geophysics, 67(2), 209-229.

-Van Dissen, R. J., Johnson, K. M., Seebeck, H., Wallace, L. M., Rollins, C., Maurer, J., ... & DiCaprio, C. J. (2024). Upper Plate and Subduction Interface Deformation Models in the 2022 Revision of the Aotearoa New Zealand National Seismic Hazard Model. Bulletin of the Seismological Society of America, 114(1), 37-56.
