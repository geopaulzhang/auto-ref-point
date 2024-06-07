[![Language](https://img.shields.io/badge/python-3.8%2B-blue?style=flat-square)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue?style=flat-square)](https://github.com/geopaulzhang/auto-ref-point/blob/main/LICENSE)
[![Citation](https://img.shields.io/badge/DOI-10.1109%2FLGRS.2024.3390568-blue?style=flat-square)](https://doi.org/10.1109/LGRS.2024.3390568)

## Automated InSAR Reference Point Selection

This repository contains a python script that automatically selects the optimal reference point for SBAS InSAR analysis. When provided multiple Areas of Interest (AOIs, e.g., polygons for a shapefile), the script outputs optimal reference points for each polygon. It requires users to install MintPy, gdal, and other basic Python packages as prerequisites. Details of the procedure of reference point selection consisting of five steps are described in the following paper:

+ Zhang, B., Hestir, E., Yunjun, Z., Reiter, M.E., Viers, J.H., Schaffer-Smith, D., Sesser, K. and Oliver-Cabrera, T., 2024. Automated Reference Points Selection for InSAR Time Series Analysis on Segmented Wetlands. _IEEE Geoscience and Remote Sensing Letters, 21,_ 4008705. https://ieeexplore.ieee.org/document/10504695.

#### Input data for the script include:

1. a shapefile for the Area of Interest (AOI) (may include one or multiple polygons) with a WGS84 geographic coordinate system
2. the products from TOPS stack processing from ISCE, i.e., coregistered spatial coherence and connected component

#### Final output products include:

1. A list of optimal reference points for each shapefile. By the default, the first reference point is selected for the rest of the SBAS procedure (Fig. 3f).

#### Intermediate output products for each of the five steps:

1. [Step 2] A map of percentage of InSAR pairs satisfying the criteria with a coherence greater than 0.9 (Fig. 3a)
2. [Step 3] For each AOI pixel, a map showing the percentage of a non-AOI pixel connected to the AOI pixel (Fig. 3b)
3. [Step 4a] A map of reference candidates (Fig. 3c)
4. [Step 4b] A connected component map for the reference candidates representing their spatial distribution (Fig. 3d)
5. [Step 5] A connected component map representing the coherence path to the AOI (Fig. 3e)

The script will be incorporated into MintPy in the near future. 

If any questions, please contact the developer Boya ("Paul") Zhang (bzhang64@ucmerced.edu). 




