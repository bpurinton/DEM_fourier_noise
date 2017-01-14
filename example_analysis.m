%% Spectral Analysis of Topographic Noise
%
% Ben Purinton (purinton[at]uni-potsdam.de), January 2017
% Universitaet Potsdam, Institute for Earth and Environmental Science,
% Germany
%
%
% This example script and the associated functions demonstrate the 2D
% discrete Fourier transform (2D DFT) and associated 1D analysis of DEM
% noise in:
% Purinton, B., and Bookhagen, B.: Validation of digital elevation models
% (DEMs) and geomorphic metrics on the southern Central Andean
% Plateau, Earth Surface Dynamics, 2017.
%
% The simple function developed for this analysis is:
% _FFTnoise.m_
%
%
% Much of the processing steps and the functions for detrending of the
% elevation matrix and windowing are taken from work by:
% Perron, J. T., Kirchner, J. W., and Dietrich, W. E.: Spectral signatures
% of characteristic spatial scales and nonfractal structure in landscapes,
% Journal of Geophysical Research, 113, 2008.
%
% These functions and example scripts are available online at
% http://web.mit.edu/perron/www/downloads.html, and include:
% _Detrend.m_, _lsplane.m_ (Obtained from http://www.eurometros.org/),
% _fft2D.m_, and _Hann2D.m_
%
%
% Functions for loading and displaying the DEMs are from the TopoToolbox
% suite in:
% Schwanghart, W., and Scherler, D.: Short Communication: TopoToolbox 2 –
% MATLAB-based software for topographic analysis and modeling in Earth
% surface sciences, Earth Surface Dynamics, 2, 1-7, 2014.
%
% These functions are available online at
% https://github.com/csdms-contrib/topotoolbox, and include: 
% _GRIDobj.m_ and _imageschs.m_
%
%
% Regressions and statistical tests rely on functions from Matlab (TM)
% Statistics and Machine Learning Toolbox.
%
%
% This script is intended to take any rectangular DEM clip and output a
% figure with subplots showing 1) hillshade, 2) 1D plot of spectral power,
% 3) 2D plot of normalized spectral power, and 4) normalized 1D plot of
% spectral power

%% Setup

% Be sure to download and set the path to all the dependencies (Perron et
% al., 2008; Schwanghart and Scherler, 2014) listed above from:

% http://web.mit.edu/perron/www/downloads.html ('2DSpecTools')

% and

% https://github.com/csdms-contrib/topotoolbox ('TopoToolbox')

%% Set Parameters

clear
clc

% For this example script we include a 30 m SRTM-C clip from Purinton and
% Bookhagen (2017). 
%
% Must be rectangular GeoTiff with no NaNs and with units in UTM meters

data = 'srtmc_30m_fft_clip.tif'; % data = 'C:\path\to\dem\dem.tif';

% Choose the necessary function parameters
rotate = 1; % 0 (default) to analyze the matrix as is, 1 to rotate the matrix by 90 degrees
nbins = 20; % number of bins for regression and normalization (default 20)
envelope_percentile = 99.9; % percentile envelope for 1D normalized power spectrum (default 99.9 highlights the outliers)

%% Call Function and Generate Plot
FFTnoise(data, rotate, nbins, envelope_percentile);
