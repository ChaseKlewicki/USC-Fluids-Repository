clc; clear; close all; 

% Written by Oliver Geoffrey Khan on 4 December 2023
%
% This script will provide a brief example of how to use the function
% computeSVPM.m to find the inviscid velocity and pressure distribution
% around an airfoil. 


% The airfoil coordinates will be contained in a .txt file that should be
% in the same folder as this script and the function computeSVPM.m. Once 
% this is the case, we can call computeSVPM.m directly on the airfoil file name.

airfoilFile = 'airfoils/NACA65_1412.txt';

% Also define the other necessary parameters for computation. 

alpha = 0; % degrees 

% Define uniform grid away from the airfoil at which to find the velocity
% and pressure.
numXGrid = 140; % Number of grid points in x direction.
numYGrid = 140; % Number of grid points in y direction. 

numPanels = 200; 

plotting = 1; % Make computeSVPM.m produce plots. 

%% Call the primary function. 
[uMat, vMat, CpMat, xGrid, yGrid, CpAirfoil, XC, YC, Cl, Cd, Cm, SVec, midIndS] =...
    computeSVPM(airfoilFile, alpha, numXGrid, numYGrid, numPanels, plotting);
%%
dCpds = gradient(CpAirfoil(midIndS+1 : end), cumsum(SVec(midIndS+1 : end)));
SVecTop = SVec(midIndS+1 : end);
figure 
for i = 2:length(SVecTop)
    plot(sum(SVecTop(1:i)), trapz(dCpds(1:i),cumsum(SVecTop(1:i))),'.')
    hold on
end
%% Print the lift coefficient. 
fprintf('At alpha = %.2f, Cl = %.3f\n', alpha, Cl)

%% Plot the velocity vector field around the airfoil. 

figure(8)

quiver(xGrid, yGrid, uMat, vMat)
hold on 

fill(XC,YC,'k') % Plot airfoil  

% Make the plot pretty. 
axis equal
ylim([-0.4, 0.4])
xlim([-0.1, 1.1])
set(gca, 'FontSize', 15)
xlabel('x/c', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('y/c', 'Interpreter', 'latex', 'FontSize', 20)
title('Velocity Field', 'Interpreter', 'latex', 'FontSize', 20)

%% Plot Cl vs. alpha. 

alphaVec = -2:2:10; % Define angle of attack array.
ClVec = zeros(1, length(alphaVec)); % Initialize Cl array.

% For each angle of attack, find the lift coefficient. (This might take
% some time since the function is not perfectly optimized to do this.)
for j = 1:length(alphaVec)
    
    [~, ~, ~, ~, ~, ~, ~, ~, ClVec(j), ~, ~, SVec] =...
    computeSVPM(airfoilFile, alphaVec(j), numXGrid, numYGrid, numPanels, 0);

    % Make sure to turn plotting off when calling the function in this loop.
end 

figure(9)

% Generate Cl vs. alpha graph 
plot(alphaVec, ClVec, 'k-', 'LineWidth', 1.5)


% Make the plot pretty. 
set(gca, 'FontSize', 15)

xlabel('$\alpha$', 'Interpreter','latex', 'FontSize', 20)

ylabel('$C_l$', 'Interpreter', 'latex', 'FontSize', 20)







