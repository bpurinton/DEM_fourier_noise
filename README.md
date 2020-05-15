---
title: "Testing for High Frequency Noise in DEMs Using FFTs"
date: "January 2017"
author: "Ben Purinton ([purinton@uni-potsdam.de](purinton@uni-potsdam.de))"
---

# CODE MIGRATION TO PYTHON :)

**Note: This analysis has been migrated to Python as of May 2020: [https://github.com/bpurinton/DEM-FFT](https://github.com/bpurinton/DEM-FFT)**  

**Proceed below if you want to work with MATLAB<sup>TM</sup>.**

___



## Matlab code for DEM noise analysis using 2D DFT

This MATLAB<sup>TM</sup> function is intended for the spectral analysis of gridded topographic data (DEMs) for analysis of high-frequency, low-wavelength noise as presented in:

Purinton, B., and Bookhagen, B.: Validation of digital elevation models (DEMs) and geomorphic metrics on the southern Central Andean Plateau, Earth Surface Dynamics, 2017.
(https://doi.org/10.5194/esurf-5-211-2017)


## Running it

Prior to analysis the user must download and set path to a few required functions:

1. The 2DSpecTools package from T. Perron available [here](http://web.mit.edu/perron/www/downloads.html)
    * For background on this spectral analysis procedure and the paper that spurred this analysis refer to: Perron, J. T., Kirchner, J. W., and Dietrich, W. E.: Spectral signatures of characteristic spatial scales and nonfractal structure in landscapes, Journal of Geophysical Research, 113, 2008. [https://doi.org/10.1029/2007JF000866](https://doi.org/10.1029/2007JF000866)

2. The TopoToolbox package from W. Schwanghart available [here](https://github.com/csdms-contrib/topotoolbox)
    * Schwanghart, W., and Scherler, D.: Short Communication: TopoToolbox 2 – MATLAB-based software for topographic analysis and modeling in Earth surface sciences, Earth Surface Dynamics, 2, 1-7, 2014. [https://doi.org/10.5194/esurf-2-1-2014](https://doi.org/10.5194/esurf-2-1-2014)


