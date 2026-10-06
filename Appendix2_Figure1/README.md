# Appendix 2 – Figure 1

This folder contains the MATLAB source code and associated source data used for DNA-PAINT nanocluster analysis shown in Appendix 2 – Figure 1.

## Files

- `Appendix2_Figure1_SourceCode1_DBSCAN.m`  
  MATLAB script used for DBSCAN-based identification and quantification of DNA-PAINT localization clusters.

- `Appendix2_Figure1_SourceData1.xlsx`  
  Source data containing cluster- and ROI-level measurements used in the analysis.

## Analysis

DNA-PAINT localization coordinates are extracted from localization files and DBSCAN clustering is applied to the XY coordinates within focal-adhesion (FA) and non-FA regions.

The analysis calculates measurements including:

- number of localizations per cluster;
- equivalent cluster diameter;
- cluster area;
- localization density;
- ROI area;
- number of clusters;
- cluster density per unit area.

DBSCAN parameters used in the analysis:

- Pixel size: 160 nm
- Epsilon: 0.10 pixels
- MinPts: 10

## Software

MATLAB (MathWorks).
