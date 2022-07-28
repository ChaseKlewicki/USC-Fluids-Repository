function [xw, yw, uw, vw, xgr, ygr, ugr, vgr] = coordinateTransform(xwall, ywall, x, ...
    y, u, v, height, varargin)
%
% Interpolates the velocity field provided by [x,y,u,v] on a grid generated
% as an offset of a 2D curve [xwall,ywall] and transforms the interpolated
% field into wall coordinates, rotating each vector according to the wall
% angle.
% Typical example: plot of the boundary layer of an airfoil, provided the
% velocity field and the airfoil points.
%
% [xwall,ywall] are the coordinates of the wall
%
% [x,y,u,v] are the points and the components of the velocity field (can be
% both structured or unstructured)
%
% [height] is the extent of the profile in the direction normal to the
% wall
%
%  Additional input arguments are:
%
% [Nh] is number of points in which the grid is discretised vertically.
% Default value is 20
%
% [Nw] is number of points in which the grid is discretised horizontally.
% By default, the points provided for the wall are used to generate the
% grid. If Nw is specified, wall points are interpolated in Nw new
% locations
%
% [interp_type]: type of interpolant ('nearest', 'linear' or 'natural'),see
% scatteredInterpolant for additional info on the interpolation method.
% Values outside the velocity field will always be linarly extrapolated
%
% [graphic_mode]: true or false, create two figures with the velocity data
% to show the result of the conversion
%
%  Output arguments
% [xw,yw,uw,vw]: are the points and the components of the velocity data
% rotated into wall coordinates. The first elements of uw and vw are the
% velocities for the first point of the provided wall. Elements for the
% successive columns are the velocities on the successive elements of the
% wall.
%
% [xgr,ygr,ugr,vgr]: are the points and the components of the velocity data
% interpolated on the wall generated grid, before the coordinate
% transformation
%
%  Additional info
% The direction in which the mesh is extended can be changed flipping the
% wall points direction. See fliplr
%
% The interpolation is performed using the scatteredInterpolant function,
% which uses a Delaunay triangulation. If no vectors are provided on the
% wall, those values will be interpolated (or extrapolated) using a
% Delaunay triangulation of the entire flow field. See help
% scatteredInterpolant for more information.
%
% Run the code without any input argument to see an example
%
% Author: Alessandro Masullo 2016 from Matlab Exchange (velocityProfile).


% Check input arguments
if nargin == 0
    % Demo mode
    % Flow around a cylinder
    [x, y] = meshgrid(-3:0.2:3);
    rh = sqrt(x.^2+y.^2);
    th = atan2(y, x);
    del = rh < 1;
    x(del) = [];
    y(del) = [];
    rh(del) = [];
    th(del) = [];
    vw = 0.5 * (1 - 1 ./ rh.^2) .* cos(th);
    vt = -0.5 * (1 + 1 ./ rh.^2) .* sin(th);
    u = vw .* cos(th) - vt .* sin(th);
    v = vw .* sin(th) + vt .* cos(th);
    % Cylinder wall
    bf = linspace(2*pi, 0);
    xwall = cos(bf);
    ywall = sin(bf);
    % Default parameters
    Nh = 20;
    Nw = 0;
    height = 1;
    interp_type = 'linear';
    graphic = true;
elseif nargin == 7
    Nh = 20;
    Nw = 0;
    interp_type = 'linear';
    graphic = false;
elseif nargin == 8
    Nh = varargin{1};
    Nw = 0;
    interp_type = 'linear';
    graphic = false;
elseif nargin == 9
    Nh = varargin{1};
    Nw = varargin{2};
    interp_type = 'linear';
    graphic = false;
elseif nargin == 10
    Nh = varargin{1};
    Nw = varargin{2};
    interp_type = varargin{3};
    graphic = false;
elseif nargin == 11
    Nh = varargin{1};
    Nw = varargin{2};
    interp_type = varargin{3};
    graphic = varargin{4};
end

% Validate arguments
validateattributes(Nh, {'numeric'}, {'numel', 1});
validateattributes(Nw, {'numeric'}, {'numel', 1});
validateattributes(interp_type, {'char'}, {'nonempty'});
validateattributes(graphic, {'logical'}, {'nonempty'});

% Prepare velocity interpolant
uint = scatteredInterpolant(x(:), y(:), u(:), interp_type, 'linear');
vint = scatteredInterpolant(x(:), y(:), v(:), interp_type, 'linear');

% Evaluate tangent and normal versors for the wall
tgv = [diff(xwall(:)), diff(ywall(:))];
module = sqrt(sum(tgv.^2, 2));
tgv = tgv ./ [module, module];
% The last versor is linearly extrapolated
tgv(end+1, :) = tgv(end-1, :) + (tgv(end, :) - tgv(end-1, :)) * 2;
nmv = [-tgv(:, 2), tgv(:, 1)];

% The versors are evalauted with the original wall points. Now consider
% the provided number of points:
if Nw == 0
    xwall_sub = xwall;
    ywall_sub = ywall;
    nmv_sub = nmv;
    tgv_sub = tgv;
else
    xwall_sub = interp1(1:numel(xwall), ...
        xwall, linspace(1, numel(xwall), Nw), '*cubic');
    ywall_sub = interp1(1:numel(ywall), ...
        ywall, linspace(1, numel(ywall), Nw), '*cubic');
    bf1 = interp1(1:numel(xwall), ...
        nmv(:, 1), linspace(1, numel(ywall), Nw), '*cubic');
    bf2 = interp1(1:numel(xwall), ...
        nmv(:, 2), linspace(1, numel(ywall), Nw), '*cubic');
    nmv_sub = [bf1(:), bf2(:)];
    bf1 = interp1(1:numel(xwall), ...
        tgv(:, 1), linspace(1, numel(ywall), Nw), '*cubic');
    bf2 = interp1(1:numel(xwall), ...
        tgv(:, 2), linspace(1, numel(ywall), Nw), '*cubic');
    tgv_sub = [bf1(:), bf2(:)];
end

% Generate the grid based on the wall. The grid is a multiple wall
% offset in the specified direction
xgr = repmat(xwall_sub(:)', Nh, 1) + ...
    kron(nmv_sub(:, 1)', linspace(0, height, Nh)');
ygr = repmat(ywall_sub(:)', Nh, 1) + ...
    kron(nmv_sub(:, 2)', linspace(0, height, Nh)');

% Interpolate velocity on the new grid points
ugr = uint(xgr, ygr);
vgr = vint(xgr, ygr);

% Rotate velocity vectors by the opposite direction of the tangent vector
mcos = repmat(tgv_sub(:, 1)', Nh, 1);
msin = -repmat(tgv_sub(:, 2)', Nh, 1);
uw = ugr .* mcos - vgr .* msin;
vw = ugr .* msin + vgr .* mcos;

% Curvilinear abscissa and structured grid
curvabs = cumsum(sqrt(diff(xwall).^2+diff(ywall).^2));
curvabs = [0; curvabs(:)];
if Nw == 0
    curvabs_sub = curvabs;
else
    curvabs_sub = interp1(1:numel(curvabs), ...
        curvabs, linspace(1, numel(curvabs), Nw), '*cubic');
end
[xw, yw] = meshgrid(curvabs_sub, linspace(0, height, Nh));

% graphic mode
if graphic
    figure, hold on
    title('Initial velocity field')
    mesh(xgr, ygr, xgr*0, 'EdgeColor', [.8, .8, .8]);
    plot(xwall, ywall, 'k', 'linewidth', 2)
    quiver(x, y, u, v, 'r')
    legend('Generated grid', 'Wall', 'Velocity field')
    axis equal
    scatter(xwall, ywall, 30, 1:numel(xwall), 'filled')

    figure, hold on
    title('Transformed velocity field')
    mesh(xw, yw, xw*0, 'EdgeColor', [.8, .8, .8]);
    quiver(xw, yw, uw, vw, 'r')
    plot(xw(:, 11:10:end)+uw(:, 11:10:end), yw(:, 11:10:end), ...
        '-k', 'linewidth', 2)
    legend('Generated grid', 'Transformed velocity field', ...
        'Velocity profile')
    axis equal
    plot([xw(1), xw(end)], [yw(1), yw(1)], 'k', 'linewidth', 2)
    scatter((xw(1, :)), yw(1, :), 30, 1:size(xw, 2), 'filled')
    xlabel('Curvilinear abscissa')
    ylabel('Normal extent')
end
