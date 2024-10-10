function [uMat, vMat, CpMat, xGrid, yGrid, CpAirfoil, XC, YC, Cl, Cd, Cm, SVec, midIndS] =...
    computeSVPM(airfoilFile, alpha, numXGrid, numYGrid, numPanels, plotting)

% Code organized by: Oliver Geoffrey Khan on 16 Nov 2023
%
% COMPUTESVPM is a function that will use the Source-Vortex Panel Method to 
% determine the velocity and pressure around an airfoil at a given 
% angle of attack. The flow is assumed to be inviscid. A given airfoil is
% divided up into panels, and along each panel, there is a source sheet of 
% constant strength and a vortex sheet of constant strength. While the
% source sheet strengths can vary from panel to panel, the vortex sheet 
% strength is the same for all panels. 
% 
% Credit for all of the algorithms goes to jte0419 on GitHub, whose
% original code can be found here: 
%   
% https://github.com/jte0419/Panel_Methods
%
% INPUTS:
%   - airfoilFile: Name of the .txt file of the airfoil coordinates 
%                  (eg. 'NACA65.txt'). This file should have no header and
%                  contain numbers only. 
%
%   - alpha: airfoil angle of attack with the freestream velocity [degrees].
%
%   - numXGrid: number of grid points in the x-direction at which to 
%               compute u, v, and Cp.
%               A good number is around 150.
%
%   - numYGrid: number of grid points in the y-direction at which to 
%               compute u, v, and Cp.
%               A good number is around 150.
%
%   - numPanels: The number of panels to model around the airfoil
%                surface. A good number is around 200.
%
%   - plotting: If plotting == 1, this function will plot diagrams of the 
%               panel normal vectors, streamlines and a pressure coefficient 
%               contour. 
%
%               If plotting == 0, this function will not make any plots. 
%
% OUTPUTS: 
%   - uMat: [numXGrid x numYGrid] matrix containing the normalized 
%           x-component of velocity at all grid points in xGrid, yGrid. 
%
%   - vMat: [numXGrid x numYGrid] matrix containing the normalized 
%           y-component of velocity at all grid points in xGrid, yGrid. 
%
%   - CpMat: [numXGrid x numYGrid] matrix containing the pressure coefficient 
%           at all grid points in xGrid, yGrid.
%
%   - xGrid: [numXGrid x numYGrid] matrix in meshgrid style containing 
%            normalized x-coordinates of all grid points in the flow. 
%
%   - yGrid: [numXGrid x numYGrid] matrix in meshgrid style containing 
%            normalized y-coordinates of all grid points in the flow.
%
%   - CpAirfoil: [numPanels x 1] vector containing the values of Cp around
%                the airfoil surface. The j-th element of CpAirfoil is the
%                value of Cp at the control point of the j-th panel. 
%
%   - XC: [numPanels x 1] vector containing the x-coordinates of panel
%         control points. This output can come in handy for producing
%         additonal plots of the pressure coefficient over the airfoil surface. 
%
%   - YC: [numPanels x 1] vector containing the y-coordinates of panel
%         control points. 
%
%   - Cl: 2D lift coefficient.
%
%   - Cd: 2D pressure drag coefficient. 
%
%   - Cm: 2D moment coefficient.
%
% ADDITIONAL INFORMATION:
%   - The velocities here are assumed to be normalized by a freestream
%     velocity, V_∞.
%   
%   - The x and y coordinates are assumed to be normalized by the airfoil 
%     chord length, c.
% 
%   - The pressure coefficient, Cp, solves Cp = (p - p_∞) / (0.5 ρ V_∞)^2.


%% Load the airfoil file

%  This program assumes that the airfoil file contains just the raw numbers
%  of the airfoil coordinates, and nothing else. 

airfoil = importdata(airfoilFile);

%% Increase airfoil resolution 

% Compute comulative distance along airfoil contour.
dist = vecnorm(diff(airfoil, 1, 1), 2, 2);

% Parametrize x and y coordinates. 
abscissa = [0; cumsum(dist)];

% Compute spline coefficients. 
pp = spline(abscissa, flip(airfoil', 2));

% New abscissa.
abscissaNew = linspace(abscissa(1), abscissa(end), numPanels + 1);

% Resample airfoil contour with custom distribution. 
airfoilRefined = ppval(pp, abscissaNew);

%Take the transpose to be consistent with original airfoil coordinate file. 
airfoilRefined = airfoilRefined'; 


%% Define the panel boundary points (XB, YB) and control points (XC, YC). 

% To begin, let's find the x-coordinates of the panel boundary points. 

% We want the panels to progress from the trailing edge, along the bottom 
% of the airfoil, to the leading edge, along the top of the airfoil, back 
% to the trailing edge. This orientation allows us to define normal vectors 
% properly. 

XB = airfoilRefined(:,1);

% Now, let's define the y-coordinates of panel boundary points.
YB = airfoilRefined(:,2);


%% Verfity panel orientation and reorient panels if they are not in the 
%% correct order.

% Check for direction of points
edge = zeros(numPanels,1);                                                     

for index = 1:numPanels                                                          
   
    edge(index) = ( XB(index+1)-XB(index) ) * ( YB(index+1)+YB(index) ); 
    
    % If the panels go from TE, along pressure surface, to LE, along suction
    % surface, back to TE, then this quantity above should be positive.

end

sumEdge = sum(edge);                                                        

% If panels are oriented counterclockwise, then flip them. 
% If they are oriented clockwise, then do not flip them. 

if sumEdge < 0                                                          
    XB = flipud(XB);                                                        
    YB = flipud(YB);                                                       
end

%% Define panel geometry parameters.

% Initialize all relevant arrays.
XC   = zeros(numPanels,1); % x-coordinate of control points                                                   
YC   = zeros(numPanels,1); % y-coordinate of control points                                                 
SVec    = zeros(numPanels,1); % Panel lengths                                               
phiDeg = zeros(numPanels,1);  % Angle between panels and the +x-axis [deg.]                                                

% Take the control points to be the midpoints of the panels, take SVec(j) 
% to be the length of panel j, and phiDeg(j) to be value of phi for panel j.

for j = 1:numPanels                                                          
    XC(j)   = XB(j) + 0.5 * (XB(j+1) - XB(j));                                         
    YC(j)   = YB(j) + 0.5 * (YB(j+1) - YB(j)); 
    
    %Find x and y displacements along panel.
    deltaX      = XB(j+1)-XB(j);                                                
    deltaY      = YB(j+1)-YB(j); 
    
    %Define S and phi values. 
    SVec(j)    = (deltaX^2 + deltaY^2)^0.5;                                          
    phiDeg(j) = atan2d(deltaY, deltaX); 

    %If phiDeg(j) is negative, then make it positive. 
    if (phiDeg(j) < 0)                                                    
        
        phiDeg(j) = phiDeg(j) + 360; 
        
        %Recall that all angles here are in degrees.
    end
end

% Compute angles of panel normal vectors with respect to the positive
% x-axis, and then compute the angles of panel normal vectors with the
% oncoming freestream flow. 

% Angle from positive x-axis to outward normal vector [degrees]
deltaDeg = phiDeg + 90; 

% Angle between freestream vector and outward normal vector [degrees]
betaDeg = deltaDeg - alpha;                                        

% Make all panel angles between 0 and 360 [degrees]
betaDeg(betaDeg > 360) = betaDeg(betaDeg > 360) - 360;                              

% Convert angles from [degrees] to [radians]
phi  = phiDeg.*(pi/180); 
delta = deltaDeg * (pi/180);
beta = betaDeg.*(pi/180);                                                 


%% Determine source panel strengths and vortex panel strength.

% Set the normalized freestream velocity equal to 1. 
Vinf = 1;

% Find the normal (I) and tangential (J) geometric integrals for the source
% panel method. 
[I,J] = computeIJ_SPM(XC,YC,XB,YB,phi,SVec);  

% Find the normal (K) and tangential (L) geometric integrals for the vortex
% panel method. 
[K,L] = computeKL_VPM(XC,YC,XB,YB,phi,SVec);              

% Construct the primary A matrix. 
A = I + pi * eye(numPanels,numPanels);

% Populate the right column of the A matrix. 
for j = 1:numPanels                                                  
    A(j,numPanels+1) = -sum(K(j,:));      
end

% Enforce the Kutta condition in the bottom row of the A matrix.
for j = 1:numPanels 
    
    % Account for the contribution of the sources to the Kutta condition.
   
    A(numPanels+1,j) = J(1,j) + J(numPanels,j) ;                                

end

%Account for the contribution of the vortex sheets to the Kutta condition.

A(numPanels+1,numPanels+1) = -sum(L(1,:) + L(numPanels,:)) + 2*pi;                  

% Define the b vector 
b = -Vinf * 2 * pi * cos(beta);

% Modify the last element of the b array to satisfy the Kutta condition.
b(numPanels+1) = -2 * pi * Vinf *( sin(beta(1)) + sin(beta(numPanels)) );          

% Compute the result array that contains the vortex strength and
% the source strengths for each of the panels. 
resultArray = A\b;                                                               

% Extract source strengths and vortex strength from the result array.
lambda = resultArray(1:end-1);                                          
gamma  = resultArray(end);                         

%% Determine the velocity and pressure coefficient along each panel. 

%Initialize tangential velocity array and pressure coefficient array.
Vtan = zeros(numPanels,1);                                               
CpAirfoil = zeros(numPanels,1);                                                  


for j = 1:numPanels
    % Determine tangential velocity and pressure coefficient for panel j.

    % Contribution from uniform flow. 
    term1 = Vinf*sin(beta(j));

    % Contribution from the source sheets on the other panels.
    term2 = (1/(2*pi))*sum(lambda .* J(j,:)');                              
    
    % Contribution from vortex sheet on panel j.
    term3 = gamma/2;   

    % Contribution from vortex sheets on the other panels. 
    term4 = -(gamma/(2*pi)) * sum(L(j,:));                                   
    
    % Define tangential velocity on panel j.
    Vtan(j) = term1 + term2 + term3 + term4;                                 
    
    %Define pressure coefficient on panel j. 
    CpAirfoil(j) = 1 - ( Vtan(j)/Vinf )^2;                                              
end


%% Compute the lift and moment coefficients.

% Determine the normal force coefficients on all the panels. Recall that
% the chord line is parallel to the x-axis and that delta is the angle 
% between a panel normal vector and the +x-axis.
CN = -CpAirfoil .* SVec .* sin(delta);                                                    

% Determine the axial force coefficients on all the panels. 
CA = -CpAirfoil .* SVec .* cos(delta);                                                     

% Compute lift, drag, and moment coefficients. 
Cl = sum(CN .* cosd(alpha)) - sum(CA .* sind(alpha));  

Cd = sum(-CN .* sind(alpha)) + sum(CA .* cosd(alpha));  

Cm = sum(CpAirfoil .* (XC - 0.25) .*SVec .* cos(phi));    


%% Calculate the streamlines away from the airfoil. 

% Define grid parameters.                                                       
xVals  = [min(XB)-0.5, max(XB)+0.5]; % x range                                
yVals  = [min(YB)-0.4, max(YB)+0.4]; % y range                          
    

% Define streamline parameters:

% Step size for propagation of streamlines. 
stepsize = 0.01;                                                       

% Maximum number of vertices. 
maxVert  = numXGrid*numYGrid*100; 

% Percentage of streamlines to be plotted on the grid.  
streamlinePercentage = 25;                                                       

% Create an array of streamline starting points parallel to the y-axis. 
yStreamLineStarts = linspace(yVals(1), yVals(2), ...
    floor( (streamlinePercentage/100)*numYGrid) )';     
    
% Generate the evenly-spaced grid points.
Xgrid   = linspace(xVals(1), xVals(2), numXGrid)'; 

Ygrid   = linspace(yVals(1), yVals(2), numYGrid)';                         

[xGrid, yGrid] = meshgrid(Xgrid, Ygrid);                                       
    
% Initialize matrices containing the normalized velocities. 
uMat = zeros(numXGrid, numYGrid);                                              
vMat = zeros(numXGrid, numYGrid);                                          
    
% Solve for grid point x and y velocities.
for m = 1:1:numXGrid
    for n = 1:1:numYGrid

        % Extract current grid point location. 
        XP      = xGrid(m,n);                                             
        YP      = yGrid(m,n);  

        % Compute Source Panel Method streamline geometric integrals. 
        [Mx,My] = streamline_SPM(XP,YP,XB,YB,phi,SVec); 

        % Compute Vortex Panel Method streamline geometric integrals. 
        [Nx,Ny] = streamline_VPM(XP,YP,XB,YB,phi,SVec);    


        % If the current grid point is on the airfoil, then set the velocity
        % equal to zero here. 
        [in,on] = inpolygon(XP,YP,XB,YB);

        if (in == 1 || on == 1)                                       
            uMat(m,n) = 0;                                               
            vMat(m,n) = 0; 
        
        % If the grid point is off the airfoil, then compute the true x-
        % and y-components of velocity.
       
        else                                                           
            uMat(m,n) = Vinf * cosd(alpha) + sum(lambda.*Mx./ (2*pi)) + ...    
                            sum(-gamma .*Nx ./ (2*pi));
                
            vMat(m,n) = Vinf * sind(alpha) + sum(lambda.*My./ (2*pi)) + ...    
                            sum(-gamma .* Ny ./ (2*pi));
        end
    end
end
    
% Determine the magnitude of velocity at the grid point and use that
% magnitude to compute the pressure coefficient. 
VMag  = sqrt(uMat.^2 + vMat.^2);                                      

CpMat = 1-(VMag ./ Vinf).^2;                                       


%% If plotting == 1, then plot streamlines and Cp contours.
%% Also plot a diagram of the panel normal vector orientations. 

if plotting  == 1

    % Plot the airfoil with representations of the normal vectors on each
    % panel

    figure();                                                              
    
    cla; hold on; grid off;                                                 
   
    set(gcf,'Color','White');                                        
    set(gca,'FontSize',15);                           
    
    fill(XB,YB,'k');                                                  
    
    for j = 1:numPanels 
        % Determine beginning and starting points of panel normal vectors
        % and then plot. 
        Xnormal(1) = XC(j);                                                      
        Xnormal(2) = XC(j) + SVec(j)*cosd(betaDeg(j)+alpha);                     
        
        Ynormal(1) = YC(j);                                                    
        Ynormal(2) = YC(j) + SVec(j)*sind(betaDeg(j)+alpha); 

        plot(Xnormal,Ynormal,'r','LineWidth',2);                                       
    end
   
    title('Panel Normal Vector Diagram', 'Interpreter','latex', 'FontSize', 20)
    
    xlabel('$x/c$', 'Interpreter','latex', 'FontSize', 20)                                                      
    ylabel('$y/c$', 'Interpreter','latex', 'FontSize', 20)                                                    
	
    xlim('auto')                                                          
    ylim('auto')                                                        
    axis equal                                                      
    zoom reset                                                           

    % Plot streamlines around the airfoil. 

    figure()

    cla; hold on; grid on;                                                 
    
    set(gcf,'Color','White');                                               
    
    set(gca,'FontSize',15);                                                
   
    % For each streamline starting point, plot streamline trajectory. 
    for index = 1:length(yStreamLineStarts)                                                 
       
        % Start each streamline from the left side of the grid and advance
        % rightward. 
        slHandle = streamline(xGrid, yGrid, uMat, vMat, xVals(1), ...
            yStreamLineStarts(index), [stepsize, maxVert] ); 
        
        % Customize streamline line width. 
        set(slHandle,'LineWidth',2);                                              
   
    end
    
    % Plot airfoil.
    fill(XB,YB,'k');                                                       
    
    title('Streamlines', 'Interpreter','latex', 'FontSize', 20)
    
    xlabel('$x/c$', 'Interpreter','latex', 'FontSize', 20);                                                      
    ylabel('$y/c$', 'Interpreter','latex', 'FontSize', 20);                                                     
    
    xlim([xVals(1), xVals(end)]);  

    axis equal;                                                            
    ylim([yVals(1), yVals(end)]);                                                          
    zoom reset;                                                            

   % Now plot the pressure coefficient contour. 

    figure()
   
    cla; hold on; grid on;                                                  
    
    set(gcf,'Color','White');                                               
    
    set(gca,'FontSize',15);                                                 
    
    % Plot Cp contour. 
    contourf(xGrid, yGrid, CpMat, 100, 'EdgeColor', 'none');                          
    
    barHandle = colorbar;

    ylabel(barHandle, '$c_{p}$', 'FontSize', 20, 'interpreter', 'latex')

    % Plot the airfoil.
    
    fill(XB,YB,'k');  
    
    xlabel('$x/c$', 'interpreter','latex', 'FontSize', 20);                                                      
    
    ylabel('$y/c$', 'interpreter', 'latex', 'FontSize', 20); 

    xlim([xVals(1), xVals(end)]);                                                          
	
    axis equal;                                                            
    
    ylim([yVals(1), yVals(end)]);                                                           
    
    zoom reset;    

    % Also plot pressure coefficient vectors along the airfoil surface 
    
    figure();

    cla
    hold on
    grid on                                              
    
    set(gcf,'Color','White');                                    
    set(gca,'FontSize',15);                                      
    
    % Scale pressure coefficients to make the plot prettier.
    CpScaled = CpAirfoil * 0.08;

    % Take the absolute value of all pressure coefficients. 
    CpAbs = abs(CpScaled); 

    for j = 1:length(CpAbs)
        
        % Extract control point x-coordinate.
        X(1) = XC(j);  
        
        % Extract ending x value based on the magnitude of Cp.
        X(2) = XC(j) + CpAbs(j) * cosd(betaDeg(j) + alpha);                           
        
        % Extract control point y-coordinate.
        Y(1) = YC(j); 

        % Define an ending y-value based on the magnitude of Cp.
        Y(2) = YC(j) + CpAbs(j) * sind(betaDeg(j) + alpha);                         
        
        if (CpAirfoil(j) < 0) 
            
            % Plot the Cp vector as a blue line if Cp is negative. 
            p{1} = plot(X, Y, 'b-', 'LineWidth', 2);                               
        
        elseif (CpAirfoil(j) >= 0)  
            
            % Plot the Cp vector as a red line if Cp is positive.
            p{2} = plot(X, Y, 'r-', 'LineWidth', 2);                                  
        
        end
    end

    % Plot airfoil. 
    fill(XB,YB,'k');                                                      
    
    % Define legend.
    legend([p{1},p{2}],{'C_{p} < 0','C_{p} > 0'});                      
    
    title('$C_{p}$ Vector Diagram', 'Interpreter', 'latex', 'FontSize', 20)
    xlabel('$x/c$', 'Interpreter','latex', 'FontSize', 20);                                                   
    ylabel('$y/c$', 'Interpreter','latex', 'FontSize', 20);                                                     
    
    xlim('auto')                                                 
    ylim('auto')                                                          
    axis equal                                                        
    zoom reset  


    % Now, plot the Cp distribution around the airfoil. 

    figure()
    cla
    hold on
    grid on

    set(gcf,'Color','White')                                               
    set(gca,'FontSize',15)   

    
    % Determine index of airfoil LE, and call it "midIndS," because it will
    % most likely be the index of the middle value in an array that
    % contains the cumulative arclength around the airfoil surface. 
    minXControlPoint = min(XC);

    midIndS = find(XC == minXControlPoint);

    % Plot Cp vs. x for the suction surface. 
    CpUpper = plot(XC(midIndS+1:end), CpAirfoil(midIndS+1 : end),'b-',...
        'LineWidth', 1.5);
    
    % Plot Cp vs. x for the pressure surface. 
    CpLower = plot(XC(1:midIndS),CpAirfoil(1:midIndS),'r-','LineWidth', 1.5);


    % Create legend and format the rest of the plot.
    legend([CpUpper,CpLower], {'Suction Surface', 'Pressure Surface'}, ...
        'FontSize', 12);
    
    title('$C_{p}$ vs. $x/c$', 'Interpreter','latex', 'FontSize', 20)
    
    xlabel('$x/c$', 'Interpreter','latex', 'FontSize', 20);                                              
    ylabel('$C_{p}$', 'Interpreter','latex', 'FontSize', 20);                                                          
    xlim([0 1]);                                                           
    ylim('auto');                                                          
    set(gca,'Ydir','reverse')  

    zoom reset; 

    % Finally, plot the Cp gradient around the airfoil. 
    figure()
    
    % We compute dCp/ds, the derivative of Cp with respect to arc length. 
    % The built-in gradient function is sufficient to do this, as we will
    % differentiate Cp with a uniform grid spacing of:
    % abscissaNew(2) - abscissaNew(1),
    % which represents the approximate arclength distance between control points. 
    
    % First separate pressure coefficient values along the suction surface
    % from pressure coefficient values along the pressure surface. 
    CpSuctionSurface = CpAirfoil(midIndS+1 : end);

    CpPressureSurface = flip(CpAirfoil(1:midIndS)); %Orient from LE to TE

    % We can now evaluate dCp/ds as we walk along the suction surface from
    % the leading edge to the trailing edge.
    dCpdsSuctionSurface = gradient(CpSuctionSurface,...
        abscissaNew(2)-abscissaNew(1));

   % Recall that abscissaNew(2)-abscissaNew(1) is the approximate arc
   % length spacing between panel control points. 


    % We can also evaluate dCp/ds as we walk along the pressure surface from
    % the leading edge to the trailing edge.

     dCpdsPressureSurface = gradient(CpPressureSurface,...
         abscissaNew(2)-abscissaNew(1));

    % Plot dCp/ds vs. x/c for the suction surface. 
    xSuctionSurface = XC(midIndS+1:end); 

    plot(xSuctionSurface,dCpdsSuctionSurface,'b-',...
        'LineWidth', 1.5);
    hold on 

    % Plot dCp/ds vs. x/c for the pressure surface. 
    xPressureSurface = flip(XC(1:midIndS));

    plot(xPressureSurface,dCpdsPressureSurface,'r-','LineWidth', 1.5);

    % Create legend and format the rest of the plot.
    legend('Suction Surface', 'Pressure Surface','FontSize', 12);

    set(gca, 'FontSize', 12)
    
    title('$\frac{dC_{p}}{ds}$ vs. $x/c$', 'Interpreter', 'latex', 'FontSize', 20)
    
    xlabel('$x/c$', 'Interpreter', 'latex', 'FontSize', 20);                                              
    ylabel('$\frac{dC_{p}}{ds}$', 'Interpreter', 'latex', 'FontSize', 25);                                                          
    
    xlim([0 1])
    ylim([-50,30])                                                        
    zoom reset


    % Now compute and plot dCpdx. This can be done through a quick
    % modification of our work above. 

    dCpdxSuctionSurface = gradient(CpSuctionSurface, xSuctionSurface);
    dCpdxPressureSurface = gradient(CpPressureSurface, xPressureSurface);
    
    % Now plot dCp/d(x/c) vs. x/c
    figure()

    plot(xSuctionSurface,dCpdxSuctionSurface,'b-',...
        'LineWidth', 1.5);
    hold on 

    plot(xPressureSurface,dCpdxPressureSurface,'r-','LineWidth', 1.5);

    % Create legend and format the rest of the plot.
    legend('Suction Surface', 'Pressure Surface','FontSize', 12);

    set(gca, 'FontSize', 12)
    
    title('$\frac{dC_{p}}{d(x/c)}$ vs. $x/c$', 'Interpreter', 'latex', 'FontSize', 20)
   
    xlabel('$x/c$', 'Interpreter', 'latex', 'FontSize', 20);                                              
    ylabel('$\frac{dC_{p}}{d(x/c)}$', 'Interpreter', 'latex', 'FontSize', 25);                                                          
    
    xlim([0 1])
    ylim([-50, 30])                                                         
    zoom reset
end 


end 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define auxillary geometric integral functions. 

function [Mx,My] = streamline_SPM(XP,YP,XB,YB,phi,SVec)

% STREAMLINE_SPM is a function that will compute the geometric integrals at
% a point (XP, YP) away from an airfoil due to source panels on the
% airfoil. 
%
% Credit for this algorithm goes to jte0419 on GitHub, whose original 
% code can be found here: 
%   
% https://github.com/jte0419/Panel_Methods
%
% A YouTube video explaining the derivations behind the following formulae
% can be found at: 
%
% https://www.youtube.com/watch?v=BnPZjGCatcg
%
% INPUTS:
% - XP: x-coordinate of computation point.
%
% - YP: y-coordinate of computation point.
%
% - XB: x-coordinate of boundary points.
%
% - YB: y-coordinate of boundary points.
%
% - phi: [numPanels x 1] column vector containing the angles between the
%        positive x-axis and interiors of the panels [radians].
%
% - SVec: [numPanels x 1] column vector containing the lengths of all the
%         panels. 
% 
% OUTPUTS:
% - Mx: [numPanels x 1] vector of x-direction source panel geometric integrals. 
%
% - My: [numPanels x 1] vector of y-direction source panel geometric integrals. 

% Extract the number of panels.
numPanels = length(XB)-1;                                                      

% Initialize Mx and My geometric integral arrays.
Mx = zeros(numPanels,1);                                                      
My = zeros(numPanels,1);                                                       

% Compute Mx and My
for j = 1:numPanels                                                      
    % Compute the appropriate terms that appear in the derivations of Mx
    % and My. 
    A  = -(XP-XB(j)) * cos(phi(j)) - (YP-YB(j))*sin(phi(j));                    
    B  = (XP-XB(j))^2 + (YP-YB(j))^2;                                        
    
    % Define C and D terms for x-direction geometric integral.
    Cx = -cos(phi(j));                                                     
    Dx = XP - XB(j);

    % Define C and D terms for y-direction geometric integral. 
    Cy = -sin(phi(j));                                                      
    Dy = YP - YB(j);                                                   
    
    E  = sqrt(B-A^2); 

    % Set E to zero if E is not real. 
    if (~isreal(E))
        E = 0;
    end
    
    % Compute Mx
    term1 = 0.5 * Cx * log((SVec(j)^2 + 2*A*SVec(j) + B)/B); 

    term2 = ((Dx-A*Cx)/E)*(atan2((SVec(j)+A),E) - atan2(A,E));                 
    
    % Compose geometric integral for x-direction. 
    Mx(j) = term1 + term2;                                                  
    
    % Compute My
    term1 = 0.5*Cy*log((SVec(j)^2 + 2*A*SVec(j) + B) / B); 

    term2 = ((Dy-A*Cy)/E) * ( atan2( (SVec(j)+A),E ) - atan2(A,E) );                 
    
    % Compose geometric integral for y-direction.
    My(j) = term1 + term2;                                                  
    
    % Set NANs, INFs, or imaginary numbers to zero in Mx and My arrays 
    if (isnan(Mx(j)) || isinf(Mx(j)) || ~isreal(Mx(j)))
        Mx(j) = 0;
    end
    
    if (isnan(My(j)) || isinf(My(j)) || ~isreal(My(j)))
        My(j) = 0;
    end
end

end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Nx,Ny] = streamline_VPM(XP,YP,XB,YB,phi,SVec)

% STREAMLINE_VPM is a function that will compute the geometric integrals at
% a point (XP, YP) away from an airfoil due to vortex panels on the
% airfoil. 
%
% Credit for this algorithm goes to jte0419 on GitHub, whose original 
% code can be found here: 
%   
% https://github.com/jte0419/Panel_Methods
%
% A YouTube video explaining the derivations behind the following formulae
% can be found at: 
%
% https://www.youtube.com/watch?v=TBwBnW87hso
%
% INPUTS
% - XP: x-coordinate of computation point.
%
% - YP: y-coordinate of computation point.
%
% - XB: x-coordinate of boundary points.
%
% - YB: y-coordinate of boundary points.
%
% - phi: A [numPanels x 1] column vector containing the angles between the 
%        positive x-axis and interior of panel [radians].
%
% - SVec: A [numPanels x 1] column vector containing the lengths of all the
%         panels. 
% 
% OUTPUTS
% - Nx: [numPanels x 1] array of x-direction vortex panel geometric integrals. 
%
% - Ny: [numPanels x 1] array of y-direction vortex panel geometric integrals. 


% Extract the number of panels.
numPanels = length(XB)-1;                                             

% Initialize Nx and Ny geometric integral arrays. 
Nx = zeros(numPanels,1);                                                   
Ny = zeros(numPanels,1);                                                  

% Calculate Nx and Ny using the derived formulae. 
for j = 1:1:numPanels                                                         
    % Compute intermediate values.
    A  = -(XP-XB(j)) * cos(phi(j)) - (YP-YB(j))*sin(phi(j));                   
    B  = (XP-XB(j))^2 + (YP-YB(j))^2;                                        
    
    % Define C and D terms for x-direction integrals.
    Cx = sin(phi(j));                                                       
    Dx = -(YP-YB(j));                                                      

    % Define C and D terms for y-direction integrals. 
    Cy = -cos(phi(j));                                                     
   
    Dy = XP-XB(j);                                                          
   
    E  = sqrt(B-A^2);                                                      
    
    % Make E = 0 if E is not real. 
    if (~isreal(E))
        E = 0;
    end
    
    % Compute Nx
    term1 = 0.5*Cx*log( (SVec(j)^2 + 2*A*SVec(j)+B)/B );                              
    term2 = ((Dx-A*Cx)/E) * ( atan2((SVec(j)+A),E) - atan2(A,E) );                 
   
    % Formulate Nx geometric integral. 
    Nx(j) = term1 + term2;                                                  
    
    % Compute Ny
    term1 = 0.5*Cy*log((SVec(j)^2+2*A*SVec(j)+B)/B);                             
    term2 = ( (Dy-A*Cy)/E )*( atan2((SVec(j)+A),E) - atan2(A,E) );                 
    
    % Formulate Ny geometric integral. 
    Ny(j) = term1 + term2;                                                  
    
	% Turn NANs, INFs, or imaginary numbers to zeros. 
    if (isnan(Nx(j)) || isinf(Nx(j)) || ~isreal(Nx(j)))
        Nx(j) = 0;
    end
    if (isnan(Ny(j)) || isinf(Ny(j)) || ~isreal(Ny(j)))
        Ny(j) = 0;
    end
end

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I,J] = computeIJ_SPM(XC,YC,XB,YB,phi,SVec)

% COMPUTEIJ_SPM is a function that will compute the normal and tangential
% geometric integral matrices for the Source Panel Method. These quantities 
% tell us how the velocity at one source panel along the airfoil is 
% influenced by the other source panels. 
%
% Credit for this algorithm goes to jte0419 on GitHub, whose original 
% code can be found here: 
%   
% https://github.com/jte0419/Panel_Methods
%
% A YouTube video explaining the derivations behind the following formulae
% can be found at: 
%
% https://www.youtube.com/watch?v=JRHnOsueic8
%
% INPUTS
% - XC: x-coordinate of control points.
%
% - YC: y-coordinate of control points.
%
% - XB: x-coordinate of boundary points.
%
% - YB: y-coordinate of boundary points.
%
% - phi: [numPanels x 1] array containing the angles between the 
%        positive x-axis and interior of panel [radians]
%
% - SVec: [numPanels x 1] array containing the lengths of the source
%         panels. 
% 
% OUTPUTS
% - I: [numPanels x numPanels] matrix of panel-normal Source Panel Method 
%      geometric integrals. 
%
% - J: [numPanels x numPanels] matrix of panel-tangential Source Panel Method 
%      geometric integrals. 

% Extract the number of panels.
numPanels = length(XC);                                                        

% Initialize I and J matrices. 
I = zeros(numPanels,numPanels);                                                   
J = zeros(numPanels,numPanels);                                                  

% Compute I and J matrices of integrals
for i = 1:1:numPanels                                                        
    for j = 1:1:numPanels                                                
        if (j ~= i)                                            
            % Compute intermediate values
            A  = -(XC(i)-XB(j)) * cos(phi(j))-(YC(i)-YB(j)) * sin(phi(j));   
            B  = (XC(i)-XB(j))^2+(YC(i)-YB(j))^2;

            % Define C and D coefficients for the panel-normal direction. 
            Cn = sin(phi(i)-phi(j));                                     
            Dn = -(XC(i)-XB(j))*sin(phi(i))+(YC(i)-YB(j)) * cos(phi(i)); 

            % Define C and D coefficients for the panel-tangential
            % direction.
            Ct = -cos(phi(i)-phi(j));                                    
            Dt = (XC(i)-XB(j))*cos(phi(i))+(YC(i)-YB(j)) * sin(phi(i));     
            
            E  = sqrt(B-A^2);                                               
            
            % Set E = 0 if E is not real. 
            if (~isreal(E))
                E = 0;
            end
            
            % Determine I (needed for normal velocity).
            term1  = 0.5*Cn*log((SVec(j)^2+2*A*SVec(j)+B)/B);                     
            term2  = ( (Dn-A*Cn)/E )*( atan2((SVec(j)+A),E) - atan2(A,E) );        
            
            I(i,j) = term1 + term2;                                        
            
            % Determine J (needed for tangential velocity).
            term1  = 0.5*Ct*log((SVec(j)^2+2*A*SVec(j)+B)/B);                     
            term2  = ( (Dt-A*Ct)/E )*( atan2((SVec(j)+A),E) - atan2(A,E) );        
            
            J(i,j) = term1 + term2;                                       
        end
        
        % Change NANs, INFs, or imaginary numbers to zeros.
        if (isnan(I(i,j)) || isinf(I(i,j)) || ~isreal(I(i,j)))
            I(i,j) = 0;
        end

        if (isnan(J(i,j)) || isinf(J(i,j)) || ~isreal(J(i,j)))
            J(i,j) = 0;
        end
    end
end

end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K,L] = computeKL_VPM(XC,YC,XB,YB,phi,SVec)

% COMPUTEKL_VPM is a function that will compute the normal and tangential
% geometric integral matrices for the Vortex Panel Method. These quantities 
% tell us how the velocity at one vortex panel along the airfoil is 
% influenced by the other vortex panels. 
%
% Credit for this algorithm goes to jte0419 on GitHub, whose original 
% code can be found here: 
%   
% https://github.com/jte0419/Panel_Methods
%
% A YouTube video explaining the derivations behind the following formulae
% can be found at: 
%
% https://www.youtube.com/watch?v=IxWJzwIG_gY
%
% INPUTS
% - XC: x-coordinate of control points.
%
% - YC: y-coordinate of control points.
%
% - XB: x-coordinate of boundary points.
%
% - YB: y-coordinate of boundary points.
%
% - phi: [numPanels x 1] array containing the angles between the
%        positive x-axis and interior of panel [radians].
%
% - SVec: [numPanels x 1] array containing the lengths of the source
%         panels. 
% 
% OUTPUTS
% - K: [numPanels x numPanels] matrix of panel-normal Vortex Panel Method 
%      geometric integrals. 
%
% - L: [numPanels x numPanels] matrix of panel-tangential Vortex Panel 
%       Method geometric integrals. 


% Extract the number of panels
numPanels = length(XC);                                                      

% Initialize arrays
K = zeros(numPanels, numPanels);                                                   
L = zeros(numPanels, numPanels);                                                  

% Compute integral
for i = 1:1:numPanels                                                          
    for j = 1:1:numPanels                                                     
        if (j ~= i)                                                         
            A  = -(XC(i)-XB(j))*cos(phi(j))-(YC(i)-YB(j))*sin(phi(j));      
            B  = (XC(i)-XB(j))^2+(YC(i)-YB(j))^2;                           
           
            % Determine normal C and D terms. 
            Cn = -cos(phi(i)-phi(j));                                    
            Dn = (XC(i)-XB(j))*cos(phi(i))+(YC(i)-YB(j))*sin(phi(i));      
            
            % Determine tangential C and D terms. 
            Ct = sin(phi(j)-phi(i));                                       
            Dt = (XC(i)-XB(j))*sin(phi(i))-(YC(i)-YB(j))*cos(phi(i));     
            
            E  = sqrt(B-A^2);                                             
            
            % If E is not real, then set E = 0.
            if (~isreal(E))
                E = 0;
            end
            
            % Compute K matrix element
            term1  = 0.5*Cn * log((SVec(j)^2+2*A*SVec(j)+B)/B);                     
            term2  = ((Dn-A*Cn)/E) * (atan2( (SVec(j)+A),E) - atan2(A,E) );         
            K(i,j) = term1 + term2;                                        
            
            % Compute L matrix element
            term1  = 0.5*Ct*log((SVec(j)^2+2*A*SVec(j)+B)/B);                     
            term2  = ((Dt-A*Ct)/E)*(atan2((SVec(j)+A),E)-atan2(A,E));         
            L(i,j) = term1 + term2;                                   
        end
        
        % Set NANs, INFs, or imaginary numbers to zero.
        if (isnan(K(i,j)) || isinf(K(i,j)) || ~isreal(K(i,j)))
            K(i,j) = 0;
        end

        if (isnan(L(i,j)) || isinf(L(i,j)) || ~isreal(L(i,j)))
            L(i,j) = 0;
        end
    end
end
end 

