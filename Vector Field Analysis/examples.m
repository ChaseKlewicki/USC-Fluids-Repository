% This script shows examples of how to use the functions in the
% Vector Field Analysis directory.
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
set(6, 'position', [1, 1, 1600, 1000])
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
set(7, 'position', [1, 1, 1600, 1000])
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

%% FTLE

xVec = xMat(1, :);
yVec = yMat(:, 1);

tStart = length(tVec);
tLength = -85;
tStep = -2;
xMinROI = min(xVec);
xMaxROI = max(xVec);
yMinROI = min(yVec);
yMaxROI = max(yVec);
ROIx = 100;
ROIy = 50;
method = 'Euler';

mymapred = [(0:0.1:1)', zeros(11, 2)];

for i = tStart:round(tLength/20):tStart + tLength

    % Compute FTLE
    [sigma, xPos, yPos] = FTLE(uMat, vMat, xVec, yVec, ...
        i, tLength, tStep, dt, ...
        xMinROI, xMaxROI, yMinROI, yMaxROI, ...
        ROIx, ROIy, method);

    % Plot FTLE field
    figure(4)
    contourf(xPos(:, :, 1), yPos(:, :, 1), sigma, 10);
    title('FTLE plot');
    colormap(mymapred)
    hold on
    fill(xMask, yMask, 'k')
    plot(xAirfoil, yAirfoil, 'w', 'LineWidth', 3)
    axis equal
    clim([4, 6])
    axis([xMinROI, xMaxROI, yMinROI, yMaxROI]);
    xlabel('x/c')
    ylabel('y/c')
    drawnow

end
