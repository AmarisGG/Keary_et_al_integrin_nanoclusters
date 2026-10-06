# Keary et al. – Differential spatial regulation and activation of integrin nanoclusters inside focal adhesions

This repository contains source code and associated source data used for quantitative image analysis and figure generation in https://doi.org/10.7554/eLife.105270.2

The repository includes analyses associated with:

- Appendix 1 – Figure 1
- Appendix 1 – Figure 2
- Appendix 2 – Figure 1

The source data associated with each analysis are provided together with the corresponding code.

## Repository contents

```text
Keary_et_al_integrin_nanoclusters/
│
├── README.md
│
├── Appendix1_Figure1/
│   ├── Appendix1_Figure1_SourceCode1.m
│   └── Appendix1_Figure1_SourceData1.xlsx
│
├── Appendix1_Figure2/
│   ├── Appendix1_Figure2_SourceCode1_Fiji.zip
│   ├── Appendix1_Figure2_SourceCode2_MATLAB.zip
│   └── Appendix1_Figure2_SourceData1.xlsx
│
└── Appendix2_Figure1/
    ├── Appendix2_Figure1_SourceCode1_DBSCAN.m
    └── Appendix2_Figure1_SourceData1.xlsx
```

---

## Appendix 1 – Figure 1

### STED fluorescence-intensity profiles

`Appendix1_Figure1_SourceCode1.m` generates the STED fluorescence-intensity line profiles associated with Appendix 1 – Figure 1.

The analysis includes profiles for:

- α5 integrin
- active β1 integrin (9EG7)
- αVβ3 integrin
- paxillin

Profiles corresponding to Glass, noFA, and FA regions are read directly from the associated Excel source-data file and plotted without smoothing or normalization.

### Files

- `Appendix1_Figure1_SourceCode1.m` – MATLAB source code
- `Appendix1_Figure1_SourceData1.xlsx` – source data used for figure generation

---

## Appendix 1 – Figure 2

### STED radial-profile and focal-adhesion edge-enrichment analysis

The analysis for Appendix 1 – Figure 2 follows a Fiji/ImageJ preprocessing workflow, followed by MATLAB analysis and figure generation.

### Fiji/ImageJ analysis

The Fiji macros generate consecutive 60-nm radial regions extending inward from a manually defined focal-adhesion boundary and measure fluorescence intensity within each radial region.

Radial profiles are generated up to 1000 nm from the FA edge.

The Fiji workflow contains:

- `01_generate_radial_bins_60nm_to1000nm.ijm`
- `02_measure_radial_mean_intensity_C2_C4.ijm`

The resulting radial-intensity measurements are exported for subsequent MATLAB analysis.

### MATLAB analysis

The MATLAB workflow:

- processes individual FA radial-intensity profiles;
- performs minimum subtraction for each individual profile;
- normalizes each profile to its maximum intensity;
- interpolates profiles onto a common radial-distance grid;
- calculates mean normalized intensity and SEM across FAs;
- reports the final radial profiles up to 600 nm from the FA edge;
- calculates an FA-level edge-enrichment metric.

Edge enrichment is defined as:

```text
mean normalized intensity at 60–180 nm
---------------------------------------
mean normalized intensity at 180–600 nm
```

The final figure-generation script produces:

- normalized radial fluorescence-intensity profiles shown as mean ± SEM;
- edge-enrichment plots displaying individual FA measurements;
- statistical comparisons using two-sided Welch's t-tests.

The source-data workbook contains analyses for paxillin, α5 integrin, active β1 integrin, and Tensin-3.

### Files

- `Appendix1_Figure2_SourceCode1_Fiji.zip` – Fiji/ImageJ macros
- `Appendix1_Figure2_SourceCode2_MATLAB.zip` – MATLAB processing and figure-generation scripts
- `Appendix1_Figure2_SourceData1.xlsx` – processed source data used for figure generation

---

## Appendix 2 – Figure 1

### DNA-PAINT nanocluster analysis

`Appendix2_Figure1_SourceCode1_DBSCAN.m` performs DBSCAN-based cluster identification and quantification of DNA-PAINT localization data within focal-adhesion (FA) and non-FA regions.

Localization coordinates are extracted from DNA-PAINT localization files and DBSCAN clustering is applied to the XY coordinates.

The analysis calculates cluster- and ROI-level measurements including:

- number of localizations per cluster;
- cluster equivalent diameter;
- cluster area;
- localization density;
- ROI area;
- number of clusters;
- cluster density per unit area.

The DBSCAN parameters used for this analysis were:

```text
Pixel size = 160 nm
Epsilon    = 0.10 pixels
MinPts     = 10
```

The associated source-data workbook contains the resulting cluster metrics used for the analysis presented in Appendix 2 – Figure 1.

### Files

- `Appendix2_Figure1_SourceCode1_DBSCAN.m` – MATLAB DBSCAN analysis
- `Appendix2_Figure1_SourceData1.xlsx` – source data containing cluster and ROI-level measurements

---

## Software

The analyses in this repository were performed using:

- MATLAB (MathWorks)
- Fiji/ImageJ

Individual scripts contain additional information regarding required input formats, analysis parameters, and generated outputs.

## Data availability

The processed source data required to reproduce the analyses and plots described above are provided within this repository.

Additional experimental and raw data associated with the study are available through the data repository associated with the publication.

## Citation

If you use the code or data contained in this repository, please cite the corresponding publication:
https://doi.org/10.7554/eLife.105270.2

## Contact

For questions regarding the analyses or source code contained in this repository, please contact the corresponding author. 
