function h = FFTnoise(data, rotate, nbins, envelope_percentile)

% Generate plot with Fourier spectral analysis of DEM noise
%
% Syntax
%
%     FFTnoise(dem, rotate, nbins, envelope_percentile)
%
% Description
%
%     FFTnoise takes a rectangular, evenly space grid of elevation values
%     and returns a figure with a 2D and 1D Fourier spectral analysis of
%     the matrix in order to assess noise that may be present at high
%     frequencies and low wavelengths.
%
% Input:
%    data    - full path to rectangular elevation matrix in GeoTIFF format
%    rotate  - set to 1 to rotate the matrix 90 degress prior to analysis, default is 0 (don't rotate)
%    nbins   - number of bins to use for normalization regression and envelope (default is 20)
%    envelope_percentile - choose a percentile for generating envelope in 1D normalized plot (default is 99.9th percentile)
%
% Output:
%    h    - 2 by 2 figure with hillshade and spectral analysis results 
%
% This function requires dependencies from:
% http://web.mit.edu/perron/www/downloads.html
% and
% https://github.com/csdms-contrib/topotoolbox
%
% Author: Ben Purinton (purinton[at]uni-potsdam.de)
% Date: January 2017
%
% References:
% Purinton, B., and Bookhagen, B.: Validation of digital elevation models
% (DEMs) and geomorphic metrics on the southern Central Andean
% Plateau, Earth Surface Dynamics, 2017.
%
% Perron, J. T., Kirchner, J. W., and Dietrich, W. E.: Spectral signatures
% of characteristic spatial scales and nonfractal structure in landscapes,
% Journal of Geophysical Research, 113, 2008.
%
% Schwanghart, W., and Scherler, D.: Short Communication: TopoToolbox 2 –
% MATLAB-based software for topographic analysis and modeling in Earth
% surface sciences, Earth Surface Dynamics, 2, 1-7, 2014.

%% Setup

% check the inputs and set rotation to default if not selected
if rotate > 1, error('rotate must be 1 or 0'), end
if nbins > 50, error('select less bins for regression (try 20)'), end
if (nargin < 2), rotate = 0; end
if (nargin < 3), nbins = 20; end
if (nargin < 4), envelope_percentile = 99.9; end

% pull out matrix of elevation values
dem = GRIDobj(data);
Z = dem.Z;
% rotate the matrix
if rotate == 1, Z = rot90(Z); end
% get dx and dy
dx = dem.cellsize;
dy = dem.cellsize;
[Ny Nx] = size(Z);


% generate a full screen figure for subplots
h = figure;
set(h, 'units', 'normalized', 'position', [0 0 1 1], 'color', 'w')

% subplot for hillshade
subplot(2,2,1)
imageschs(dem, dem, 'colormap', 'gray')
c = colorbar;
ylabel(c, 'elevation [m]')
xlabel('Easting [m]')
ylabel('Northing [m]')
title('DEM Clip')


%% Detrend the matrix

Zo = Z; % Save the original elevations
Z = Detrend(Z);

%% Take 2D FFT (function from http://web.mit.edu/perron/www/downloads.html)

% we always using padding a windowing, see Perron et al. (2008) for details
pad = 1; window = 1;
[Pm fm Pv fv] = fft2D(Z, dx, dy, pad, window);
% outputs periodogram (P) and frequencies (f) as 2D matrix (m) and 1D vector (v)

wv = 1./fv; % convert frequencies to wavelengths

%% Background Spectrum

% bin the 1D data
x = log10(wv(:));
y = log10(Pv(:));
% sort x and y by x
sorted = sortrows([x y],1);
x = sorted(:,1); y = sorted(:,2);
% find the extrema of x
xmin = x(1); xmax = x(end);
xrange = xmax - xmin;
% determine the window width
w = xrange/nbins; 
% Allocate memory for the binned data
B = zeros(nbins,6);
% loop through the bins
for i=1:nbins
    % determine min and max x values of current window position
    xlo = xmin+(i-1)*w; % for windows with no overlap
    xhi = xlo+w;
    % find min and max indices of x vector corresponding to this range
	window = find((x >= xlo) & (x <= xhi));
    mini = min(window); maxi = max(window);
    % calculate [center of bin, mean, std, n, max, min] points that fall
    % within this window, but watch out for windows with only one point:
    if isempty(window)
        B(i,:) = [mean([xlo xhi]) 0 0 0 0 0];
    else
        B(i,:) = [mean([xlo xhi]) mean(y(mini:maxi)) ...
                 std(y(mini:maxi)) maxi-mini+1 ...
                 max(y(mini:maxi)) min(y(mini:maxi))];
    end
end

% generate a regression through the background spectrum bins
fit = robustfit(B(:,1),B(:,2), 'bisquare');

% subplot with non-normalized 1D power spectrum and background spectrum
subplot(2,2,2)
loglog(wv,Pv,'or','markersize',1) % plot the DFT elements
hold on
loglog(10.^B(:,1),10.^B(:,2), 'ok','markerfacecolor','w') % plot binned values
loglog(10.^B(:,1),10^fit(1)*(10.^B(:,1)).^fit(2),'k') % plot regression
text(10^B(1,1) + 1000, 10^B(1,2), sprintf('regression slope = %0.2f', fit(2)), 'FontSize', 10, 'BackgroundColor', 'w') % add regression slope
legend('DFT Element', 'Log-Bins for Regression', 'Regression for Normalization', 'Location', 'SouthWest')
set(gca, 'xdir', 'reverse')
grid on
xlabel('Wavelength (m)')
ylabel('DFT mean squared amplitude (m^2)')
title('Non-normalized 1D Power Spectrum')

%% Normalized 1D and 2D Spectra

% Use the 1D fit to normalize the 1D and 2D spectrum for plotting
Pmn = Pm./(10^fit(1)*fm.^fit(2));
Pvn = Pv./(10^fit(1)*wv.^fit(2));

%% Plot normalized 2D power spectrum

subplot(2,2,3)
imagesc(log10(Pmn))
axis image
c = colorbar;
colormap(c, 'gray')
ylabel(c, 'Normalized Spectral Power')
xlabel('x frequency (m^{-1})')
ylabel('y frequency (m^{-1})')
[nfy nfx] = size(fm);
nyq = fm(nfy/2+1,1); % the Nyquist frequency in the x direction
numticks=8;
set(gca, 'XTick', linspace(1,nfx,numticks), 'XTickLabel',...
    sprintf('%0.3f|', linspace(-nyq,nyq,numticks)));
set(gca, 'YTick', linspace(1,nfy,numticks), 'YTickLabel',...
    sprintf('%0.3f|', linspace(nyq,-nyq,numticks)));
set(gca,'TickDir','out')
set(gca, 'Fontsize', 8)
grid on
title('Normalized 2D Spectral Power')

%% Plot normalized 1D power spectrum

% first generate bins for 99.9th percentile envelope
x = log10(wv(:));
y = Pvn(:);
nbin = 50; % using 50 log bins to generate the envelope
% sort x and y by x
sorted = sortrows([x y],1);
x = sorted(:,1); y = sorted(:,2);
% find the extrema of x
xmin = x(1); xmax = x(end);
xrange = xmax - xmin;
% determine the window width
w = xrange/nbin; 
% Allocate memory for the binned data
B = zeros(nbin,2);
% loop through the bins
for i=1:nbin
    % determine min and max x values of current window position
    xlo = xmin+(i-1)*w; % for windows with no overlap
    xhi = xlo+w;
    % find min and max indices of x vector corresponding to this range
	window = find((x >= xlo) & (x <= xhi));
    mini = min(window); maxi = max(window);
    % calculate [center of bin, and choosen envelope percentile] points
    % that fall within this window, but watch out for windows with only one
    % point:
    if isempty(window)
        B(i,:) = [mean([xlo xhi]) 0];
    else
        B(i,:) = [mean([xlo xhi]) prctile(y(mini:maxi), envelope_percentile)];
    end
end

% plot
subplot(2,2,4)
semilogx(wv,Pvn,'or','markersize',2);
hold on
semilogx(10.^B(:,1),B(:,2),'k-', 'LineWidth', 2)
set(gca, 'xdir', 'reverse')
grid on
legend('DFT Element', sprintf('%0.2f Percentile Envelope', envelope_percentile), 'Location', 'NorthWest')
ylabel('Normalized Spectral Power')
xlabel('Wavelength (m)')
title('Normalized 1D Spectral Power')

