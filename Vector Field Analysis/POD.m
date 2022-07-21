function [Zeta, Sigma, Xi] = POD(uMat, vMat, dx, dy, dt)
% Function which performs a proper orthognal decomposition of a time
% series of 2-d vector fields. The function assumes a uniform grid. Nans
% in the vector field data are converted to zeros

% Inputs:
% uMat: Matrix of x-component vectors with dimensions [ny, nx, nt]
% vMat: Matrix of y-component vectors with dimensions [ny, nx, nt]
% dx: x direction grid spacing [scalar]
% dy: y direction grid spacing [scalar]
% dt: time between vector fields [scalar] (s)

% Outputs:
% Zeta: Martix of spatial Basis Functions:  [ny*nx, nt]
% Sigma: Martix of singular values representing energy of basis functions
% if nt > ny*nx: [nt, nt]          if nt < ny*nx: [ny*nx, ny*nx]
% Xi: Martix of temporal coefficients
% if nt > ny*nx: [nt, nt]          if nt < ny*nx: [ny*nx, ny*nx]

% Author: Mitual Luhar & Chase Klewicki

% Determine dimensions
[ny, nx, nt] = size(uMat);

% Reshape u and v matrices from 3 to 2 dimensions
u = reshape(uMat,[nx*ny,nt]);
v = reshape(vMat,[nx*ny,nt]);

% Combine u and v to process simultaneously
U = [u; v];
U(isnan(U)) = 0;

% Scaling matrices -- in this case, we have uniform grids
Sx = sqrt(dx*dy)*speye(2*nx*ny);
St = sqrt(dt)*speye(nt);

% Rescale U
Us = Sx*U*St;

% Compute POD modes using SVD
[Psi,Sigma,Phi] = svd(Us,"econ");
Zeta = Sx\Psi;
Xi = St\Phi;

% % Compute POD modes using Eig
% K = Us'*Us;
% [Phi,Lambda]=eig(K);
% [Lambda,inds]=sort(diag(Lambda),'descend');
% Lambda = diag(Lambda);
% Phi = Phi(:,inds);
% Psi = Us*Phi*diag(1./sqrt(diag(Lambda)));

end