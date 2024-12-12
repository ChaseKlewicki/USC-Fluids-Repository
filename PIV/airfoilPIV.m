% This script shows examples of how to use PIVlab functions for flow around
% a fixed wing airfoil. Requires PIVlab app from Matlab Exchange!

% Authors: Chase Klewicki (main)

clear
close all
clc

%% PIV analysis with PIVlab functions

% folder containing images
folder = ["Example Data\"];

% Camera frames per second
fps = 50;

% Chord length (m)
c = 0.2;

% Angle of Attack (rads)
alpha = 4 * pi / 180;

%% Calibration Preprocessing

% Load calibration image
calibration = imread(folder+"cal.tif");

% Change from 12 bit to 16 bit
calibration = bitshift(calibration, 4);

% Rotation (deg)
rot = -2;

% Rotate Image
calibration = imrotate(calibration, rot);

figure(1)
imshow(calibration)
h = gca;
h.Visible = 'On';
xlabel('Pixels')
ylabel('Pixels')

%% Calibration ROI

% Define ROI
% bw = roipoly(calibration); %% Manually with mouse

xMin = 250;
xMax = 800;
yMin = 50;
yMax = 550;

xPix = [xMin, xMax, xMax, xMin, xMin];
yPix = [yMin, yMin, yMax, yMax, yMin];
bw = poly2mask(xPix, yPix, size(calibration, 1), size(calibration, 2));

% Change values outside mask to max value
calibration(~bw) = 2^16;


% Plot 12 bit image
figure(1)
imshow(calibration)
h = gca;
h.Visible = 'On';
xlabel('Pixels')
ylabel('Pixels')

%% Locate calibration circles

% Calibration plate circle radii range
radius = [5, 10];

% Find white circles on calibration plate
[centers, radii, metric] = imfindcircles(calibration, radius, ...
    "Sensitivity", .925, 'ObjectPolarity', 'bright', ...
    'EdgeThreshold', 0.05);

% Filter circles based on radii
% centers = centers(radii - mean(radii)<std(radii)*1.5,:);

% number of rows of circles on calibration plate
bins = 21;

% Determine size of each bin
binSize = (max(centers(:, 2)) - min(centers(:, 2))) / (bins - 1);

% Edges of bins
edges = min(centers(:, 2)) - binSize / 2:binSize:max(centers(:, 2)) + ...
    binSize / 2;

% bin data based on its y location
bin = discretize(centers(:, 2), edges);

% Plot identified circles based on their bin
colors = lines;
figure(1)
for i = 1:bins
    viscircles(centers(bin == i, :), radii(bin == i), "color", ...
        colors(i, :));
end

% Sort based on y bin and arrange in ascending x order
sorted = sortrows([centers, bin, radii], [3, 1]);
centers = sorted(:, 1:2);
radii = sorted(:, 4);
bin = sorted(:, 3);

% Function for computing median euclidean distance
euc = @(x) median(sqrt(diff(x(:, 1)).^2+diff(x(:, 2)).^2), 'omitnan');


% Compute mean of median euclidian distance between points in each bin
% to determine pixels between dots in same columns
pixel = mean(splitapply(euc, centers, bin));

%% Default calibration for regualr 2-d PIV
% Physical distance between dots in the same column (m)
dist = .01;

% PIV scaling pix / m
scaling = pixel / dist;

%% Mask

% Load NACA 65 (1) 412 airfoil
airfoil = load(['C:\Users\charl\Google Drive\School', ...
    '\Graduate\PhD\', ...
    'Experimentation\Data\NACA 65 (1) 412', ...
    '\NACA 65(1)412 contour.txt']);

% Split into x and y coords
xAirfoil = airfoil(:, 1);
yAirfoil = airfoil(:, 2);

% Add translation from origin (Manually determined from processing image
% in PIVLab and plotting on image)
xOffset = 0.023;
yOffset = 0.028;

% Rotate coordinates based on AoA
xAirfoilRot = xAirfoil * cos(alpha-0.5*pi/180) + ...
    yAirfoil * sin(alpha-0.5*pi/180);
yAirfoilRot = -xAirfoil * sin(alpha-0.5*pi/180) + ...
    yAirfoil * cos(alpha-0.5*pi/180);

% Translate coordinates and scale for 10 cm chord
xAirfoilTranslated = xAirfoilRot * c + xOffset;
yAirfoilTranslated = yAirfoilRot * c + yOffset;

% Extract coordinates from top of airfoil
xAirfoilTop = flip(xAirfoilTranslated(1:26, :));
yAirfoilTop = flip(yAirfoilTranslated(1:26, :));

% Create mask (m)
%     xAirfoilTop(1) * scaling
%     xAirfoilTop(end) * scaling
%     xAirfoilTop(1) * scaling
xMask = [; ...
    0; ... % In front of airoil
    xAirfoilTop * scaling; ... % Airfoil profile
    size(calibration, 2); ...
    0; ...
    0; ... % Bottom left corner of photo
    ];

%     0
%     yAirfoilTop(end)/3 * scaling
yMask = [; ...
    yAirfoilTop(1) / 1.5 * scaling; ...
    yAirfoilTop * scaling; ...
    0; ...
    0; ...
    yAirfoilTop(1) / 1.5 * scaling; ...
    ];


image = imread(folder+"Img000000.tif");
% image = imrotate(image, rot);

% Select ROI
xROI = 300;
yROI = 25;
xWidth = 406;
yHeight = 1235;

% Plot first image
figure(2)
imshow(image);
hold on
plot(xMask, yMask, 'r')
plot(xAirfoilTop*scaling, yAirfoilTop*scaling, 'm')

%% Compute Mean Intensity

% Create list of image names by listing all files in folder starting with
% img and ending with .tif
imageFiles = dir(folder+'img*.tif');

% Determine number of images
nFiles = length(imageFiles);

% Compute mean intensity

% Initalize matrix of images
pics = zeros([size(image), nFiles]);

% Loop through file names
for ii = 1:nFiles
    % Get current file name
    currentFilename = imageFiles(ii).name;
    % Read in image
    pics(:, :, ii) = imread(folder+currentFilename);
end

% Compute mean image intensity
meanIntensity = mean(pics, 3);

% clear pics variable for memory
clear pics

%% Test Settings

tic

% image preprocessing settings

% Region of interest: [x,y,width,height] in pixels, may be left empty
% ROI=[xROI-16, yROI-16 xWidth+32, yHeight+32];
ROI = []
% 1 = enable CLAHE (contrast enhancement), 0 = disable
CLAHE = 1;
% CLAHE window size
CLAHE_size = 64;
% 1 = enable highpass, 0 = disable
Highpass = 0;
% highpass size
Highpass_size = 15;
% 1 = enable clipping, 0 = disable
Clipping = 0;
% 1 = enable Wiener2 adaptive denoise filter, 0 = disable
Wiener = 1;
% Wiener2 window size
Wiener_size = 2;
% Minimum intensity of input image (0 = no change)
Minimum_intensity = 0;
% Maximum intensity on input image (1 = no change)
Maximum_intensity = 0.5;

% Test settings
index = 3;
% load image A
imageA = imread(folder+imageFiles(index).name);
% imageA = imrotate(imread(folder + imageFiles(index).name), rot);
% load image B
imageB = imread(folder+imageFiles(index+1).name);
% imageB = imrotate(imread(folder + imageFiles(index+1).name), rot);
% Preprocess images
image1 = PIVlab_preproc(imageA-uint16(meanIntensity), ROI, CLAHE, ...
    CLAHE_size, Highpass, Highpass_size, Clipping, Wiener, Wiener_size, ...
    Minimum_intensity, Maximum_intensity);
image2 = PIVlab_preproc(imageB-uint16(meanIntensity), ROI, CLAHE, ...
    CLAHE_size, Highpass, Highpass_size, Clipping, Wiener, Wiener_size, ...
    Minimum_intensity, Maximum_intensity);


% Show the preprocessed image with mean intensity subtracted
figure(3)
imshow(image1)

% PIV Settings

% window size of first pass
Int_area_1 = 64;
% step of first pass
Step_size_1 = 32;
% Subpixel interpolation: 1 = 3point Gauss, 2 = 2D Gauss
Subpix_finder = 1;
% Mask
Mask = {xMask, yMask};
%
% 1-4 nr. of passes
Num_of_passes = 4;
% second pass window size
Int_area_2 = 32;
% third pass window size
Int_area_3 = 32;
% fourth pass window size
Int_area_4 = 32;
% '*spline' is more accurate, but slower
Window_deformation = '*spline';
% Repeat the correlation four times and multiply the correlation matrices
% 0 or 1 : Slows down processing by factor of 4-5
Repeated_Correlation = 0;
% 0 or 1 : Disable Autocorrelation in the first pass.
Disable_Autocorrelation = 0;
% 0 or 1 : Use circular correlation (0) or linear correlation (1).
Correlation_style = 0;
% Return Correlation Matrices
do_correlation_matrices = 1;
% 0 or 1 : Repeat the last pass of a multipass analyis
Repeat_last_pass = 1;
% Repetitions of last pass will stop when the average difference to
% the previous pass is less than this number.
Last_pass_quality_slope = 0.025;

% PIV
[x, y, u, v, typevector, correlation_map, ...
    correlation_matrices] = ...
    piv_FFTmulti(image1, image2, Int_area_1, Step_size_1, ...
    Subpix_finder, Mask, ROI, Num_of_passes, Int_area_2, Int_area_3, ...
    Int_area_4, Window_deformation, Repeated_Correlation, ...
    Disable_Autocorrelation, Correlation_style, ...
    do_correlation_matrices, Repeat_last_pass, Last_pass_quality_slope);

% Validation Settings
% calibrated velocities
calu = u / scaling * fps;
calv = u / scaling * fps;
% Valid velocity limits [min max]
valid_vel = [];
% Validate based on std (0 off, 1 on)
do_stdev_check = 1;
% Number of std valid data resides within
stdthresh = 8;
% Validate based on local median (0 off, 1 on)
do_local_median = 0;
% Number of neighbors to include in local median
neigh_thresh = 4;
% Validate based on correlation coeff (0 off, 1 on)
do_corr_check = 0;
% Correlation Threshold
corr_thresh = 0;

% Velocity validation
[u_valid, v_valid] = PIVlab_postproc(u, v, calu, calv, valid_vel, ...
    do_stdev_check, stdthresh, do_local_median, neigh_thresh);

% Correlation validation
if do_corr_check
    u_valid(correlation_map < corr_thresh) = nan;
    v_valid(correlation_map < corr_thresh) = nan;
end

% Interpolate values not in mask
uInterp = inpaint_nans(u_valid, 4);
vInterp = inpaint_nans(v_valid, 4);

% Check if a mask exists
if ~isempty(Mask)
    % Remove interpolated values in mask
    [in, on] = inpolygon(x, y, xMask, yMask);
    uInterp(in) = nan;
    vInterp(in) = nan;
end


% Plot vector fields on image
figure(4)
imshow(image1)
hold on
% quiver(x, y, uInterp, vInterp,'r')
quiver(x, y, u_valid, v_valid, 'g')

plot(xAirfoilTop*scaling, yAirfoilTop*scaling, 'm')
% h = images.roi.Polygon(gca,'Position',[xMask,yMask]);
toc

figure(5)
contourf(x, y, correlation_map, 'linestyle', 'none');
colorbar
% set(gca, 'clim', [0, 1])

%% Analyze all Images

% Numbe of images to increment by
step = 1;

% Initialize outputs
uOrig = zeros([size(u), nFiles - step]);
vOrig = zeros([size(v), nFiles - step]);
correlation = zeros([size(u), nFiles - step]);
uFilt = zeros([size(u), nFiles - step]);
vFilt = zeros([size(v), nFiles - step]);

% Loop through images
for ii = 1:1:nFiles - step
    tic
    % load image A
    imageA = imread(folder+imageFiles(ii).name);
    % imageA = imrotate(imread(folder + imageFiles(index).name), rot);
    % load image B
    imageB = imread(folder+imageFiles(ii+1).name);
    % imageB = imrotate(imread(folder + imageFiles(index+1).name), rot);
    % Preprocess images
    image1 = PIVlab_preproc(imageA-uint16(meanIntensity), ROI, CLAHE, ...
        CLAHE_size, Highpass, Highpass_size, Clipping, Wiener, ...
        Wiener_size, ...
        Minimum_intensity, Maximum_intensity);
    image2 = PIVlab_preproc(imageB-uint16(meanIntensity), ROI, CLAHE, ...
        CLAHE_size, Highpass, Highpass_size, Clipping, Wiener, ...
        Wiener_size, ...
        Minimum_intensity, Maximum_intensity);

    % PIV%
    [x, y, u, v, typevector, correlation_map, ...
        correlation_matrices] = ...
        piv_FFTmulti(image1, image2, Int_area_1, Step_size_1, ...
        Subpix_finder, Mask, ROI, Num_of_passes, ...
        Int_area_2, Int_area_3, ...
        Int_area_4, Window_deformation, Repeated_Correlation, ...
        Disable_Autocorrelation, Correlation_style, ...
        do_correlation_matrices, Repeat_last_pass, ...
        Last_pass_quality_slope);

    % Velocity validation
    [u_valid, v_valid] = PIVlab_postproc(u, v, calu, calv, valid_vel, ...
        do_stdev_check, stdthresh, do_local_median, neigh_thresh);

    % Correlation validation
    if do_corr_check
        u_valid(correlation_map < corr_thresh) = nan;
        v_valid(correlation_map < corr_thresh) = nan;
    end

    % Interpolate values not in mask
    uInterp = inpaint_nans(u_valid, 4);
    vInterp = inpaint_nans(v_valid, 4);

    % Check if a mask exists
    if ~isempty(Mask)
        % Remove interpolated values in mask
        [in, on] = inpolygon(x, y, xMask, yMask);
        uInterp(in) = nan;
        vInterp(in) = nan;
    end

    % Save data
    uOrig(:, :, ii) = u;
    vOrig(:, :, ii) = v;
    uFilt(:, :, ii) = uInterp;
    vFilt(:, :, ii) = vInterp;
    correlation(:, :, ii) = correlation_map;

    % Display progress
    disp(num2str(ii/nFiles))
    toc

    % Update figure as processing
    % (Slows down processing mucho! Comment out when not using)
    %     figure(40)
    %     clf
    %     imshow(image1)
    %     hold on
    %     quiver(x, y, uInterp, vInterp,'r')
    %     quiver(x, y, u_valid, v_valid,'g')
    %     plot(xAirfoilTop * scaling, yAirfoilTop * scaling,'m')
    %     drawnow
end
disp(num2str((ii + 1)/nFiles))

% Determine x and y values within mask
[in, ~] = inpolygon(x, y, xMask, yMask);

%% Animate correlation values

figure(3)
for i = 1:size(uOrig, 3)
    clf
    contourf(x, y, correlation(:, :, i), linspace(-1, 1, 10), ...
        'LineStyle', 'None')
    clim([0, 1])
    hold on
    fill(xAirfoilTranslated*scaling, yAirfoilTranslated*scaling, 'k')
    drawnow
end

%% Save data

data.x = x;
data.y = y;
data.uOrig = uOrig;
data.vOrig = vOrig;
data.correlation = correlation;
data.vFilt = vFilt;
data.uFilt = uFilt;
data.xMask = xMask;
data.yMask = yMask;
data.scaling = scaling;
data.fps = fps;
data.chord = c;
data.AoA = alpha;
data.xAirfoil = xAirfoilTranslated;
data.yAirfoil = yAirfoilTranslated;

save('data.mat', 'data')