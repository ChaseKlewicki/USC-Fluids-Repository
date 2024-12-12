% This script shows examples of how to use the functions in the
% Vector Field Analysis directory.

% Authors: Chase Klewicki (main)
clear
close all
clc

%% Load Data
load('Example Data/data.mat');
load('Example Data/airfoil.mat');

%% coordinateTransform

% Wall normal height for curvilinear transform (in normalized units)
height = 0.2;

% Number of points in wall normal direction (attempt to make same
% resolution as PIV grid)
Nh = 20;

% Number of points in wall parallel direction (attempt to make same
% resolution as PIV grid)
Nw = 40;

% Time average vector fields
uBar = mean(uMat, 3);
vBar = mean(vMat, 3);

% Transform coordinate system to be normal to airfoil profile
[xw, yw, uw, vw, xgr, ygr, ugr, vgr] = coordinateTransform(xAirfoil, ...
    yAirfoil, xMat, yMat, uBar, vBar, height, Nh);

% plot normalized time averaged velocties with curvilinear mesh

figure(6)
set(6, 'position', [1, 1, 800, 500])
subplot(1, 2, 1)
contourf(xMat, yMat, uBar, 'linestyle', 'none')
hold on
fill(xMask, yMask, 'k')
plot(xAirfoil, yAirfoil, 'r', 'LineWidth', 3)
mesh(xgr, ygr, xgr*0, 'EdgeColor', 'k', 'FaceAlpha', '0');
axis equal
ylim([min(yMat, [], 'all'), max(yMat, [], 'all')])
xlim([min(xMat, [], 'all'), max(xMat, [], 'all')])
xlabel('x/c')
ylabel('y/c')
a = colorbar;
ylabel(a, 'u/u_\infty')
set(gca, 'clim', [min(vBar, [], 'all'), max(uBar, [], 'all')])

subplot(1, 2, 2)
contourf(xMat, yMat, vBar, 'linestyle', 'none')
hold on
fill(xMask, yMask, 'k')
plot(xAirfoil, yAirfoil, 'r', 'LineWidth', 3)
mesh(xgr, ygr, xgr*0, 'EdgeColor', 'k', 'FaceAlpha', '0');
axis equal
ylim([min(yMat, [], 'all'), max(yMat, [], 'all')])
xlim([min(xMat, [], 'all'), max(xMat, [], 'all')])
xlabel('x/c')
ylabel('y/c')
a = colorbar;
ylabel(a, 'v/u_\infty')
set(gca, 'clim', [min(vBar, [], 'all'), max(uBar, [], 'all')])

% Plot of contours of time averaged velocities in curvilinear coordinates

figure(7)
set(7, 'position', [1, 1, 800, 500])
subplot(2, 1, 1)
contourf(xw, yw, uw, 'linestyle', 'none')
hold on
mesh(xw, yw, xw*0, 'EdgeColor', 'k', 'FaceAlpha', '0.25', 'edgeAlpha', '0.25');
axis equal
xlabel('x/c')
ylabel('y/c')
a = colorbar;
ylabel(a, 'u/u_\infty')
set(gca, 'clim', [min(vBar, [], 'all'), max(uBar, [], 'all')])


subplot(2, 1, 2)
contourf(xw, yw, vw, 'linestyle', 'none')
hold on
mesh(xw, yw, xw*0, 'EdgeColor', 'k', 'FaceAlpha', '0.25', 'edgeAlpha', '0.25');
axis equal
xlabel('x/c')
ylabel('y/c')
a = colorbar;
ylabel(a, 'v/u_\infty')
set(gca, 'clim', [min(vBar, [], 'all'), max(uBar, [], 'all')])

%% POD

dx = abs(mean(diff(xMat, [], 2), 'all'));
dy = abs(mean(diff(yMat, [], 1), 'all'));
dt = abs(mean(diff(tVec, [], 1), 'all'));

% Compute POD
[Zeta, Sigma, Xi] = POD(uMat, vMat, dx, dy, dt);

% Compute approximations at various ranks
U01 = Zeta(:, 1) * Sigma(1, 1) * Xi(:, 1)';
U03 = Zeta(:, 1:3) * Sigma(1:3, 1:3) * Xi(:, 1:3)';
U10 = Zeta(:, 1:10) * Sigma(1:10, 1:10) * Xi(:, 1:10)';
U500 = Zeta(:, 1:500) * Sigma(1:500, 1:500) * Xi(:, 1:500)';

%% DMD

% Number of DMD mode to compute 
r=25;

[eigenvalueMat, modeMat, uModeMat, vModeMat] = DMD(uMat, vMat, ...
    dx, dy, dt, r);

% Plot eigenvalues
figure
plot(diag(eigenvalueMat), '.')
hold on
% Plot unit circle
th = 0:pi/50:2*pi;
plot(cos(th),sin(th), '--k')

%% FTLE

% Sort in ascending order
[yMat_sorted, index] = sort(yMat);

% Create x and y vectors
xVec = xMat(1, :);
yVec = yMat_sorted(:, 1);

% Rearrange velocities for yMat sort
uMesh = flipud(uMat);
vMesh = flipud(vMat);

% FTLE parameters (See function for description)
tLength = -85;
tStep = -2;
xMinROI = min(xVec);
xMaxROI = max(xVec);
yMinROI = min(yVec);
yMaxROI = max(yVec);
ROIx = 200;
ROIy = 100;
method = 'Euler';

% First Frame to compute FTLE Field
frameStart = length(tVec);
% Last Frame to compute FTLE Field
frameEnd = length(tVec) - 85;
% Frame Increment
frameInc = -1;
% Frame loop vector
fLoop = frameStart:frameInc:frameEnd;

% Preallocate space for sigma
sigma = zeros([ROIx, ROIy, length(fLoop)]);

for tStart = fLoop
    [sigma(:, :, tStart), xPos, yPos] = FTLE(uMesh, vMesh, xVec, yVec, ...
        tStart, tLength, tStep, dt, ...
        xMinROI, xMaxROI, yMinROI, yMaxROI, ...
        ROIx, ROIy, method, 'xMask', xMask, 'yMask', yMask-dy);
end

% Remove empty fields in sigma
sigma = sigma(:, :, any(sigma, [1, 2]));

%% Animate FTLE field

figure(9)
% Color map
colormap(hot)
% Color limits
cmin = 5;
cmax = 10;

for i = 1:size(sigma, 3)
    clf
    % Plot FTLE field
    contourf(xPos(:, :, 1), yPos(:, :, 1), sigma(:, :, i), ...
        linspace(cmin, cmax, 10), 'linestyle', 'none');
    colorbar
%     clim([cmin cmax])
    hold on
    fill(xMask, yMask, 'k')
    plot(xAirfoil, yAirfoil, 'w', 'LineWidth', 2)
    axis equal
    axis([xMinROI, xMaxROI, yMinROI, yMaxROI]);
    xlabel('x/c')
    ylabel('y/c')
    drawnow
end

% plot trajectory of last frame
figure(10)
for j = 1:size(xPos, 3)
    plot(xPos(:, :, 1), yPos(:, :, 1), 'b.')
    hold on
    plot(xPos(:, :, j), yPos(:, :, j), 'r.')
    hold off
    axis equal
    axis([xMinROI, xMaxROI, yMinROI, yMaxROI]);
    drawnow;
end
