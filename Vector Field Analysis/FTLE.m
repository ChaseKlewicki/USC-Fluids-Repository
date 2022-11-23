function [sigma, xPos, yPos] = FTLE(uMesh, vMesh, xVec, yVec, ...
    tStart, tLength, tStep, dt, xMinROI, xMaxROI, ...
    yMinROI, yMaxROI, nx, ny, method, options)
% Function for computing 2-D finite time Lypanov exponents from a
% time series of vector fields. At the start time the flow map is
% intialized and its deformation is computed by integrating the vector
% field time series. The FTLEs are determined based on the
% stretching of the of the flow map.
%(The flow map can be thought of as seeding the flow with particles at
% the start time and observing their trajectory as time passes)

% Required Inputs
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
% xMinROI: Minimum x value for flow map
% xMaxROI: Maximum x value for flow map
% yMinROI: Minimum y value for flow map
% yMaxROI: Maximum y value for flow map
% nx: Number of grid points in ROI x-direction
% ny: Number of grid points in ROI y-direction
% Method: Integration method for determining particle trajectories
% (see trajectory function)

% Optional inputs
% xMask: Vector of x locations for closed polygon mask
% yMask: Vector of y locations for closed polygon mask
% extrap: Logical indicating whether values outside vector domain
%   should be extrapolated
% uExtrap: Function handle for extrapolation function for u velocities.
%    Function should be in the form u = u(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest u velocity value
% vExtrap: Function handle for extrapolation function for v velocities.
%    Function should be in the form v = v(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest v velocity value

% Outputs
% sigma: Field of FTLE values. [nx x ny x t_length]
% xPos: X position of particles as they are advected.
% To plot FTLE values use first time index [nx x ny x t_length]
% yPos: Y position of particles as they are advected.
% To plot FTLE values use first time index [nx x ny x t_length]

% Authors: Chase Klewicki (main), Oliver Kahn (RKF45)


arguments
    uMesh double
    vMesh double
    xVec double{mustBeVector}
    yVec double{mustBeVector}
    tStart(1, 1) double
    tLength(1, 1) double
    tStep(1, 1) double
    dt(1, 1) double
    xMinROI(1, 1) double = min(xVec)
    xMaxROI(1, 1) double = max(xVec)
    yMinROI(1, 1) double = max(yVec)
    yMaxROI(1, 1) double = max(yVec)
    nx(1, 1) double{mustBeInteger} = 100
    ny(1, 1) double{mustBeInteger} = 100
    method(1, :) char ...
        {mustBeMember(method, {'Euler', 'RKF45', 'ODE45'})} = 'Euler'
    options.xMask double{mustBeVector}
    options.yMask double{mustBeVector}
    options.extrap(1, 1) {mustBeNumericOrLogical} = false
    options.uExtrap function_handle
    options.vExtrap function_handle
end

% If extrapolation is true
if options.extrap

    % Define the initial region of the fluid to track
    sxcoor = linspace(xMinROI, xMaxROI, nx);
    sycoor = linspace(yMinROI, yMaxROI, ny);

    % Initial position
    [xPos, yPos] = ndgrid(sxcoor, sycoor);

    % Check for mask
    if isfield(options, 'xMask') && isfield(options, 'yMask')
        % find positions inside mask and set to nan
        [in, ~] = inpolygon(xPos, yPos, options.xMask, options.yMask);
        xPos(in) = nan;
        yPos(in) = nan;
    end

    % Integration time length
    tSpan = abs(tLength) * dt;

    % Create time vector for loop
    tLoop = tStart:tStep:tStart + tLength;

    % Initalize position matrix
    xPos = repmat(xPos, 1, 1, length(tLoop));
    yPos = repmat(yPos, 1, 1, length(tLoop));

    % Define time vector
    tVec = (0:sign(tStep):tStep) * dt;

    for t = tLoop

        % Determine the index
        tIndex = find(tLoop == t);

        % Compute trajectories
        [xPos(:, :, tIndex+1), yPos(:, :, tIndex+1)] = ...
            trajectory(xVec, yVec, tVec, ...
            uMesh(:, :, t:sign(tStep)*1:t+tStep), ...
            vMesh(:, :, t:sign(tStep)*1:t+tStep), ...
            xPos(:, :, tIndex), yPos(:, :, tIndex), method, ...
            options.extrap, @options.uExtrap, @options.vExtrap);
    end

    % Determine components of Jacobian
    [dPhiXdY, dPhiXdX] = gradient(xPos(:, :, end), sycoor, sxcoor);

    [dPhiYdY, dPhiYdX] = gradient(yPos(:, :, end), sycoor, sxcoor);

    % Construct Jacobian for each element
    A1 = cat(3, dPhiXdX, dPhiXdY);

    A2 = cat(3, dPhiYdX, dPhiYdY);

    % Create 4-D matrix of jacobians
    B = cat(4, A1, A2);

    % Change order of dimensions for page transpose
    B = permute(B, [3, 4, 1, 2]);

    % Compute stretching
    delta = pagemtimes(B, pagetranspose(B));

    % Initialize FTLE value
    sigma = zeros(nx, ny);

    % Compute FTLE for each point
    for j = 1:numel(sigma)
        sigma(j) = log(max(eig(delta(:, :, j)))) / abs(tSpan);
    end


    % If extrapolation is false stop integrating trajectories when
    % particle exits vector field domain
else

    % Find bounds of vector field data
    xMin = min(xVec);
    xMax = max(xVec);
    yMin = min(yVec);
    yMax = max(yVec);

    % Define the initial region of the fluid to track
    sxcoor = linspace(xMinROI, xMaxROI, nx);
    sycoor = linspace(yMinROI, yMaxROI, ny);

    % Initial position
    [xPosTemp, yPosTemp] = ndgrid(sxcoor, sycoor);

    % Check for mask
    if isfield(options, 'xMask') && isfield(options, 'yMask')
        % find positions inside mask and set to nan
        [in, ~] = inpolygon(xPosTemp, yPosTemp, ...
            options.xMask, options.yMask);
        xPosTemp(in) = nan;
        yPosTemp(in) = nan;
    else
        % If no mask exists set all points to be outside mask
        in = zeros(size(xPosTemp));
    end

    % Initialize FTLE value
    sigma = zeros(nx, ny);

    % Integration time length
    tSpan = abs(tLength) * dt;

    % Initialize flags
    calcFTLE = ones(nx, ny, 'logical');
    outDomain = zeros(nx, ny, 'logical');

    % Create time vector for loop
    tLoop = tStart:tStep:tStart + tLength;

    % Initalize position matrix
    xPos = repmat(zeros(nx, ny), 1, 1, length(tLoop));
    yPos = repmat(zeros(nx, ny), 1, 1, length(tLoop));

    % Define time vector
    tVec = (0:sign(tStep) * 1:tStep) * dt;

    for t = tLoop

        tIndex = find(tLoop == t);

        xPos(:, :, tIndex) = xPosTemp;
        yPos(:, :, tIndex) = yPosTemp;

        % If the point is inside, calculate it's trajectory
        [xPosTemp(~outDomain), yPosTemp(~outDomain)] = ...
            trajectory(xVec, yVec, tVec, ...
            uMesh(:, :, t:sign(tStep)*1:t+tStep), ...
            vMesh(:, :, t:sign(tStep)*1:t+tStep), ...
            xPosTemp(~outDomain), yPosTemp(~outDomain), method);


        % If the mask exists
        if isfield(options, 'xMask') && isfield(options, 'yMask')
            % Check for positions inside mask
            [in, ~] = inpolygon(xPosTemp, yPosTemp, ...
                options.xMask, options.yMask);
        end

        % Determine particles that have trajectories outside the domain
        index = (xPosTemp - xMin) .* (xPosTemp - xMax) >= 0 | ...
            (yPosTemp - yMin) .* (yPosTemp - yMax) >= 0 | in;

        % Change flag to indicate they have departed
        outDomain(index) = true;

        % Find adjacent points
        M = zeros(size(xPosTemp));
        M(index) = 1;
        index = conv2(M, [1, 1, 1; 1, 1, 1; 1, 1, 1], 'same') > 0;

        % Find row and column of points
        [ix, iy] = find(index);

        % Loop through points that have exited domain and their neighbors
        for i = 1:numel(ix)
            % Compute FTLE
            sigma(ix(i), iy(i)) = calculateFTLE(ix(i), iy(i), ...
                xPosTemp, yPosTemp, ...
                tSpan, sxcoor, sycoor);
        end

        % Change flag to indicate these points FTLEs have been computed
        calcFTLE(index) = false;
    end

    % calculate FTLE for all the points inside the domain at the last
    % time step.

    % Find row and column of points
    [ix, iy] = find(calcFTLE);

    % Loop through points that have exited domain and their neighbors
    for i = 1:numel(ix)
        % Compute FTLE
        sigma(ix(i), iy(i)) = calculateFTLE(ix(i), iy(i), ...
            xPosTemp, yPosTemp, ...
            tSpan, sxcoor, sycoor);
    end

end
end

function [X, Y] = trajectory(xVec, yVec, tVec, uMesh, vMesh, x0, y0, ...
    method, extrapolate, uExtrap, vExtrap)
% Function which computes trajectory of particles in a grid by
% integrating vector field of velocity values.

% Inputs
% xVec: Vector of x grid values for vector field [n x 1] (Must be in
%   ascending order)
% yVec: Vector of y grid values for vector field [m x 1] (Must be in
%   ascending order)
% tVec: Vector of time values for vector field  [p x 1]
%   (End time can be negative or positive.)
% uMesh: Matrix of x-component of velocity [m x n x p] Meshgrid style
% vMesh: Matrix of y-component of velocity [m x n x p] Meshgrid style
% x0: Matrix of intial x-position [m x n]
% y0: Matrix of intial y-position [m x n]
% method:
%   'Euler': A forward Euler integration scheme using cubic
%   interpolation of 2-D velocity field and tVec as integration interval
%   (~5x Faster than 'RKF45' with error O(h^2))
%   'RKF45': Runge–Kutta–Fehlberg method, error O(h^5)
%   'ODE45': Uses ode45 Matlab function and 3-D cubic interpolation to
%   compute 4th order Runge Kutta integration.
%   (extremely slow, error O(h^5))
% extrapolate: Logical indicating whether values outside vector domain
%   should be extrapolated
% uExtrap: Function handle for extrapolation function for u velocities.
%    Function should be in the form u = u(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest u velocity value
% vExtrap: Function handle for extrapolation function for v velocities.
%    Function should be in the form v = v(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest v velocity value

%
% Outputs
% X: x-component of trajectory [m x n]
% Y: y-component of trajectory [m x n]

arguments
    xVec double{mustBeVector}
    yVec double{mustBeVector}
    tVec double{mustBeVector}
    uMesh double
    vMesh double
    x0 double
    y0 double
    method(1, :) char ...
        {mustBeMember(method, {'Euler', 'RKF45', 'ODE45'})} = 'Euler'
    extrapolate(1, 1) {mustBeNumericOrLogical} = false
    uExtrap function_handle = @NOP
    vExtrap function_handle = @NOP
end

% Initialize length
X = zeros(size(x0));
Y = zeros(size(y0));


% NaNs propogate NanNs due to cubic interpolation! Highly suggested
% that all NaNs are removed from data prior to running this function.
% Replace NaNs with mean velocity.
uMesh(isnan(uMesh)) = 0;
vMesh(isnan(vMesh)) = 0;

[tVecSorted, indices] = sort(tVec);

%Define interpolation functions
uInterp = griddedInterpolant({xVec, yVec, tVecSorted}, ...
    permute(uMesh(:, :, indices), [2, 1, 3]), ...
    'cubic', 'nearest');
vInterp = griddedInterpolant({xVec, yVec, tVecSorted}, ...
    permute(vMesh(:, :, indices), [2, 1, 3]), ...
    'cubic', 'nearest');


    function u = uFunction(t, x, y)

        % Interpolate to velocity find values
        % (nearest value if outside domain)
        u = uInterp(x, y, t);

        % Create polygon of domain
        domain = [xVec(1), yVec(1); ...
            xVec(1), yVec(end); ...
            xVec(end), yVec(end); ...
            xVec(end), yVec(1); ...
            xVec(1), yVec(1)];

        % Determine values inside domain
        [in, on] = inpolygon(x, y, domain(:, 1), domain(:, 2));

        % Conditions to be outside domain
        conditions = ~in & ~on & ~isnan(x);

        % If extrapolation is true
        if extrapolate
            u(conditions) = uExtrap(x(conditions), y(conditions), ...
                t(conditions), u(conditions));
            % Else do not advect particles outside vector field domain
        else
            u(conditions) = 0;
        end

    end

    function v = vFunction(t, x, y)
        % Interpolate to velocity find values
        % (nearest value if outside domain)
        v = vInterp(x, y, t);

        % Create polygon of domain
        domain = [xVec(1), yVec(1); ...
            xVec(1), yVec(end); ...
            xVec(end), yVec(end); ...
            xVec(end), yVec(1); ...
            xVec(1), yVec(1)];

        % Determine values inside domain
        [in, on] = inpolygon(x, y, domain(:, 1), domain(:, 2));

        % Conditions to be in domain
        conditions = ~in & ~on & ~isnan(x);

        % If extrapolation is true
        if extrapolate
            v(conditions) = vExtrap(x(conditions), y(conditions), ...
                t(conditions), v(conditions));
            % Else do not advect particle
        else
            v(conditions) = 0;
        end

    end


% Compute trajectory using Euler method with integration step
% equal to time vector and 2-D cubic interpolation
if strcmp(method, 'Euler')
    for i = 1:length(tVec)
        x0 = x0 + ...
            tVec(i) * uFunction(tVec(i)*ones(size(x0)), x0, y0);
        y0 = y0 + ...
            tVec(i) * vFunction(tVec(i)*ones(size(x0)), x0, y0);
    end
    X = x0;
    Y = y0;
end

% Compute trajectory using ODE45 and 3-D cubic interpolation
if strcmp(method, 'ODE45')
    tic

    for ind = 1:numel(x0)
        % Create trajectory function with third order interpolation
        xTrajectory = @(t, x) uFunction(t, x, y0(ind));
        yTrajectory = @(t, y) vFunction(t, x0(ind), y);

        % Determine x and y trajectory Runge-Kutta 4th order
        % integration
        [~, x] = ode45(xTrajectory, tVec, x0(ind));
        [~, y] = ode45(yTrajectory, tVec, y0(ind));

        % assign trajectory at end time
        X(ind) = x(end);
        Y(ind) = y(end);
    end
    toc
end

% Compute trajectory using ODE45 and 3-D cubic interpolation
if strcmp(method, 'RKF45')

    [m, n] = size(x0);

    % Flatten data
    x0 = x0(:);
    y0 = y0(:);

    % RKF 45 butcher tableau

    A = [0, 2 / 9, 1 / 3, 3 / 4, 1, 5 / 6];

    B = [0, 0, 0, 0, 0, 0; ...
        2 / 9, 0, 0, 0, 0, 0; ...
        1 / 12, 1 / 4, 0, 0, 0, 0; ...
        69 / 128, -243 / 128, 135 / 64, 0, 0, 0; ...
        -17 / 12, 27 / 4, -27 / 5, 16 / 15, 0, 0; ...
        65 / 432, -5 / 16, 13 / 16, 4 / 27, 5 / 144, 0];

    CH = [47 / 450, 0, 12 / 25, 32 / 225, 1 / 30, 6 / 25];

    %     CT = [-1/150, 0, 3/100, -16/75, -1/20, 6/25];

    kx = zeros([length(x0), length(CH)]);
    ky = zeros([length(x0), length(CH)]);

    % Compute integration steps
    deltaT = gradient(tVec);

    % Loop through integration times
    for index = 1:length(tVec) - 1

        % index intergations step
        h = deltaT(index);

        % Compute integration constants
        kx(:, 1) = h * uFunction(tVec(index)*ones(size(x0))+ ...
            A(1)*h, ...
            x0, ...
            y0);
        ky(:, 1) = h * vFunction(tVec(index)*ones(size(x0))+ ...
            A(1)*h, ...
            x0, ...
            y0);

        kx(:, 2) = h * uFunction(tVec(index)*ones(size(x0))+A(2)*h, ...
            x0+kx(1, :)*B(2, :)', ...
            y0+ky(1, :)*B(2, :)');
        ky(:, 2) = h * vFunction(tVec(index)*ones(size(x0))+A(2)*h, ...
            x0+kx(1, :)*B(2, :)', ...
            y0+ky(1, :)*B(2, :)');

        kx(:, 3) = h * uFunction(tVec(index)*ones(size(x0))+A(3)*h, ...
            x0+kx(2, :)*B(3, :)', ...
            y0+ky(2, :)*B(3, :)');
        ky(:, 3) = h * vFunction(tVec(index)*ones(size(x0))+A(3)*h, ...
            x0+kx*B(3, :)', ...
            y0+ky*B(3, :)');

        kx(:, 4) = h * uFunction(tVec(index)*ones(size(x0))+A(4)*h, ...
            x0+kx(3, :)*B(4, :)', ...
            y0+ky(3, :)*B(4, :)');
        ky(:, 4) = h * vFunction(tVec(index)*ones(size(x0))+A(4)*h, ...
            x0+kx(3, :)*B(4, :)', ...
            y0+ky(3, :)*B(4, :)');

        kx(:, 5) = h * uFunction(tVec(index)*ones(size(x0))+A(5)*h, ...
            x0+kx(4, :)*B(5, :)', ...
            y0+ky(4, :)*B(5, :)');
        ky(:, 5) = h * vFunction(tVec(index)*ones(size(x0))+A(5)*h, ...
            x0+kx(4, :)*B(5, :)', ...
            y0+ky(4, :)*B(5, :)');

        kx(:, 6) = h * uFunction(tVec(index)*ones(size(x0))+A(6)*h, ...
            x0+kx(5, :)*B(6, :)', ...
            y0+ky(5, :)*B(6, :)');
        ky(:, 6) = h * vFunction(tVec(index)*ones(size(x0))+A(6)*h, ...
            x0+kx(5, :)*B(6, :)', ...
            y0+ky(5, :)*B(6, :)');

        % Compute integration
        x0 = x0 + (CH * kx')';
        y0 = y0 + (CH * ky')';

        % Compute error estimate
        %         if errorOption == true
        %         truncationErrorX(:, index) = abs(CT(1)*k1x + ...
        %             CT(2)*k2x + CT(3)*k3x...
        %             + CT(4)*k4x + CT(5)*k5x + CT(6)*k6x);
        %         truncationErrorY(:, index) = abs(CT(1)*k1y + ...
        %             CT(2)*k2y + CT(3)*k3y...
        %             + CT(4)*k4y + CT(5)*k5y + CT(6)*k6y);
        %         end
    end

    % Assign final positions to output
    X = reshape(x0, m, n);
    Y = reshape(y0, m, n);

end
end


% Function calculateFTLE

function out = calculateFTLE(ix, iy, flowmap_x, flowmap_y, tSpan, ...
    sxcoor, sycoor)
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
[nx, ny] = size(flowmap_x);

% If the location is not on the edge of the flow map
if (ix - 1) * (ix - nx) < 0 && (iy - 1) * (iy - ny) < 0
    % Compute and assemble jacobian
    dPhixdX = (flowmap_x(ix+1, iy) - flowmap_x(ix-1, iy)) ...
        / (sxcoor(ix+1) - sxcoor(ix-1));
    dPhixdY = (flowmap_x(ix, iy+1) - flowmap_x(ix, iy-1)) ...
        / (sycoor(iy+1) - sycoor(iy-1));
    dPhiydX = (flowmap_y(ix+1, iy) - flowmap_y(ix-1, iy)) ...
        / (sxcoor(ix+1) - sxcoor(ix-1));
    dPhiydY = (flowmap_y(ix, iy+1) - flowmap_y(ix, iy-1)) ...
        / (sycoor(iy+1) - sycoor(iy-1));
    A = [dPhixdX, dPhixdY; dPhiydX, dPhiydY];
    % Determine max eigenvalue of streching and compute FTLE
    delta = A' * A;
    out = log(max(eig(delta))) / abs(tSpan);
    % If the location is on the edge of the flow map
else
    out = 0;
end


end

function NOP(varargin)
%NOP Do nothing
%
% NOP( ... )
%
% A do-nothing function for use as a placeholder when working with
%  callbacks or function handles.

% Intentionally does nothing
end