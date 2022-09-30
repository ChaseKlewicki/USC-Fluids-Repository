function [sigma, xPos, yPos] = FTLE(uMesh, vMesh, xVec, yVec, ...
    tStart, tLength, tStep, ...
    dt, xMinROI, xMaxROI, ...
    yMinROI,yMaxROI, ...
    nx ,ny, method, xMask, yMask)
% Function for computing 2-D finite time Lypanov exponents from a
% time series of vector fields. At the start time the flow map is
% intialized and its deformation is computed by integrating the vector
% field time series. The FTLEs are determined based on the
% stretching of the of the flow map.
%(The flow map can be thought of as seeding the flow with particles at
% the start time and observing their trajectory as time passes)

% Inputs
% xVec: Vector of x grid values for vector field [n x 1]
% yVec: Vector of y grid values for vector field [m x 1]
% uMesh: Matrix of x-component of velocity [m x n x p] Meshgrid style
% vMesh: Matrix of y-component of velocity [m x n x p] Meshgrid style
% t_start: starting time step [scalar] (index)
% t_length: length of integration in time steps [scalar] (index length)
% NOTE: Must be negative for backward integration and less than or equal to
% p dimension
% t_step: increment of time steps (index steps) NOTE: Must be
% negative for backward integration
% dt: time between time steps (Seconds)
% xMinROI: Minimum x value for bounding box of region of interest
% xMaxROI: Maximum x value for bounding box of region of interest
% yMinROI: Minimum y value for bounding box of region of interest
% yMaxROI: Maximum y value for bounding box of region of interest
% nx: Number of grid points in ROI x-direction
% ny: Number of grid points in ROI y-direction
% Method: Integration method for determining particle trajectories
% (see trajectory function)
% xMask: vector of x locations for closed polygon mask
% yMask: vector of y locations for closed polygon mask

% Outputs
% sigma: Field of FTLE values. [nx x ny x t_length]
% xPos: X position of particles as they are advected.
% To plot FTLE values use first time index [nx x ny x t_length]
% yPos: Y position of particles as they are advected.
% To plot FTLE values use first time index [nx x ny x t_length]

% Author: Chase Klewicki (Derived from Dabiri Labs LCS Matlab Kit)


% Find bounds of vector field data
xMin=min(xVec);
xMax=max(xVec);
yMin=min(yVec);
yMax=max(yVec);

% define the initial region of the fluid to track
sxcoor=linspace(xMinROI,xMaxROI,nx);
sycoor=linspace(yMinROI,yMaxROI,ny);

% initial position
[xPosTemp, yPosTemp]=ndgrid(sxcoor,sycoor);

% Check for mask
if exist('xMask','var') && exist('yMask','var')
    % find positions inside mask and set to nan
    [in, ~] = inpolygon(xPosTemp, yPosTemp, xMask, yMask);
    xPosTemp(in) = nan;
    yPosTemp(in) = nan;
else
    % If no mask exists set all points to be outside mask
    in = zeros(size(xPosTemp));
end

% initialize FTLE value
sigma=zeros(nx,ny);

% integration time length
tSpan=abs(tLength)*dt;

% Initialize flags
calcFTLE=ones(nx,ny,'logical');
outDomain=zeros(nx,ny,'logical');

% Create time vector for loop
tLoop=tStart:tStep:tStart + tLength;

% Initalize position matrix
xPos = repmat(zeros(nx,ny),1,1,length(tLoop));
yPos = repmat(zeros(nx,ny),1,1,length(tLoop));

% Define time vector
tVec = (0:sign(tStep)*1:tStep)*dt;

for t = tLoop

    tIndex = find(tLoop==t);

    xPos(:,:,tIndex) = xPosTemp;
    yPos(:,:,tIndex) = yPosTemp;

    %%%%%% if the point is inside, calculate it's trajectory
    [xPosTemp(~outDomain), yPosTemp(~outDomain)] = ...
        trajectory(xVec, yVec, tVec, ...
        uMesh(:,:,t:sign(tStep)*1:t + tStep), ...
        vMesh(:,:,t:sign(tStep)*1:t + tStep), ...
        xPosTemp(~outDomain), yPosTemp(~outDomain), method);


    % If the mask exists
    if exist('xMask','var') && exist('yMask','var')
        % Check for positions inside mask
        [in, ~] = inpolygon(xPosTemp, yPosTemp, xMask, yMask);
    end

    % Determine particles that have trajectories outside the domain
    index = (xPosTemp-xMin).*(xPosTemp-xMax)>=0 | ...
        (yPosTemp-yMin).*(yPosTemp-yMax)>=0 | in;

    % Change flag to indicate they have departed
    outDomain(index) = true;

    % Find adjacent points
    M = zeros(size(xPosTemp));
    M(index) = 1;
    index = conv2(M,[1,1,1;1,1,1;1,1,1],'same')>0;

    % Find row and column of points
    [ix,iy] = find(index);

    % Loop through points that have exited domain and their neighbors
    for i = 1:numel(ix)
        % Compute FTLE
        sigma(ix(i),iy(i)) = calculateFTLE(ix(i),iy(i), ...
            xPosTemp,yPosTemp, ...
            tSpan, sxcoor, sycoor);
    end

    % Change flag to indicate these points FTLEs have been computed
    calcFTLE(index) = false;

    % Print progress
    c_proc = strcat( 'process accomplished :  ', ...
        num2str( 100 * tIndex/abs(tLength/tStep), ...
        ' %03.0f' ),'/100');
    disp( c_proc );
end

% calculate FTLE for all the points inside the domain at the last
% time step.

% Find row and column of points
[ix,iy] = find(calcFTLE);

% Loop through points that have exited domain and their neighbors
for i = 1:numel(ix)
    % Compute FTLE
    sigma(ix(i),iy(i)) = calculateFTLE(ix(i),iy(i), ...
        xPosTemp,yPosTemp, ...
        tSpan, sxcoor, sycoor);
end

end


function [X,Y] = trajectory(xVec, yVec, tVec, uMesh, vMesh, x0, y0, method)
% Function which computes trajectory of particles in a grid by
% integrating vector field of velocity values.

% Inputs
% xVec: Vector of x grid values for vector field [n x 1]
% yVec: Vector of y grid values for vector field [m x 1]
% tVec: Vector of time values for vector field   [p x 1]
% NOTE: Time vector must start at zero and end at desired end time.
% (End time can be negative or positive.)
% uMesh: Matrix of x-component of velocity [m x n x p] Meshgrid style
% vMesh: Matrix of y-component of velocity [m x n x p] Meshgrid style
% x0: Matrix of intial x-position [m x n]
% y0: Matrix of intial y-position [m x n]
% method:
% 'Euler': A forward Euler integration scheme using cubic
% interpolation of 2-D velocity field and tVec as integration interval
% (Much faster and less accurate than 'RK45')
% 'RK45': Uses ode45 Matlab function and 3-D cubic interpolation to
% compute 4th order Runge Kutta integration.

% Outputs
% X: x-component of trajectory [m x n]
% Y: y-component of trajectory [m x n]

% Initialize length
X = zeros(size(x0));
Y = zeros(size(y0));


% NaNs propogate NanNs due to cubic interpolation! Highly suggested
% that all NaNs are removed from data prior to running this function.
% Replace NaNs with zeros.
uMesh(isnan(uMesh))=0;
vMesh(isnan(vMesh))=0;

% Create grid for interpolation
[xLoc,yLoc,tLoc] = meshgrid(xVec,yVec,tVec);


% Compute trajectory using Euler method with integration step
% equal to time vector and 2-D cubic interpolation
if strcmp(method,'Euler')
    for i = 1:length(tVec)
        x0 = x0 + tVec(i) * interp2(xLoc(:,:,i),yLoc(:,:,i), ...
            uMesh(:,:,i), x0, y0, 'cubic', 0);
        y0 = y0 + tVec(i) * interp2(xLoc(:,:,i),yLoc(:,:,i), ...
            vMesh(:,:,1), x0, y0, 'cubic', 0);
    end
    X = x0;
    Y = y0;
end

% Compute trajectory using ODE45 and 3-D cubic interpolation
if strcmp(method,'RK45')
    tic

    for ind = 1:numel(x0)
        % Create trajectory function with third order interpolation
        xTrajectory =@(t,x) interp3(xLoc,yLoc,tLoc,uMesh,...
            x0(ind),y0(ind),t,'cubic',0);
        yTrajectory =@(t,y) interp3(xLoc,yLoc,tLoc,...
            vMesh,x0(ind),y0(ind),t,'cubic',0);

        % Determine x and y trajectory Runge-Kutta 4th order
        % integration
        [~,x] = ode45(xTrajectory, tVec, x0(ind));
        [~,y] = ode45(yTrajectory, tVec, y0(ind));

        % assign trajectory at end time
        X(ind) = x(end);
        Y(ind) = y(end);
    end
    toc
end
end


%Function calculateFTLE

function out=calculateFTLE(ix, iy,flowmap_x,flowmap_y,tSpan, sxcoor,sycoor)
% Function to compute FTLE of specfic flow map location.

% Inputs
% ix: Row index of flow map
% iy: Column index of flow map
% flowmap_x: Matrix of x position of flow map [n x m]
% flowmap_y: Matrix of y position of flow map [n x m]
% tSpan: Length of integration (scalar time units (e.g. seconds))
% sxcoor: Vector of initial x positions [n x 1]
% sycoor: Vector of initial y positions [m x 1]

% Output
% out: Finite time Lyapunov exponent

% Determine size of flow map
[nx, ny]=size(flowmap_x);

% If the location is not on the edge of the flow map
if (ix-1)*(ix-nx)<0 && (iy-1)*(iy-ny)<0
    % Compute and assemble jacobian
    dPhixdX=(flowmap_x(ix+1,iy)-flowmap_x(ix-1,iy))...
        /(sxcoor(ix+1)-sxcoor(ix-1));
    dPhixdY=(flowmap_x(ix,iy+1)-flowmap_x(ix,iy-1))...
        /(sycoor(iy+1)-sycoor(iy-1));
    dPhiydX=(flowmap_y(ix+1,iy)-flowmap_y(ix-1,iy))...
        /(sxcoor(ix+1)-sxcoor(ix-1));
    dPhiydY=(flowmap_y(ix,iy+1)-flowmap_y(ix,iy-1))...
        /(sycoor(iy+1)-sycoor(iy-1));
    A=[dPhixdX dPhixdY;dPhiydX dPhiydY];
    % Determine max eigenvalue of streching and compute FTLE
    delta=A'*A;
    out=log(max(eig(delta)))/abs(tSpan);
    % If the location is on the edge of the flow map
else
    out=0;
end
end