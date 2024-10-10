function [eigenvalueMat, modeMat, uModeMat, vModeMat] = DMD(uMat, vMat, dx, dy, dt, r)
%DMD is a function that performs exact dynamic mode decomposition of a time
%series of two-dimensional vector fields. 
%
% INPUTS:
% uMat: Matrix of x-component velocity vectors with dimensions [ny, nx, nt]
% vMat: Matrix of y-component velocity vectors with dimensions [ny, nx, nt]
% dx: x direction grid spacing [scalar]
% dy: y direction grid spacing [scalar] 
% dt: time between vector fields [scalar] 
% r: The number of DMD modes we wish to compute [scalar]
%
% OUTPUTS:
% eigenvalueVec: Matrix with dimensions [r, r] containing the DMD eigenvalues 
%                on the diagonal. 
% modeMat: Matrix with dimensions [2*nx*ny, r] containing DMD modes in
%          columns.
% uModeMat: Matrix with dimensions [ny, nx, r] that formulates our r
%           modes as r different fields of u values. 
% vModeMat: Matrix with dimensions [ny, nx, r] that formulates our r
%           modes as r different fields of v values. 

% Authors: Oliver Kahn (main)


%Determine dimensions of uMat and vMat (which should have the same
%dimensions)
[ny, nx, nt] = size(uMat);

%Reshape u and v matrices from 3 to 2 dimensions. 
xVelocities = reshape(uMat, [nx * ny, nt]);
yVelocities = reshape(vMat, [nx * ny, nt]);

%Combine u and v to process simultaneously and set all NaN's to zeros. 
combinedVelocities = [xVelocities; yVelocities];
combinedVelocities(isnan(combinedVelocities)) = 0; 

%Let's also re-scale our combinedVelocities matrix. 
Sx = sqrt(dx*dy) * speye(2*nx*ny);
St = sqrt(dt) * speye(nt);

scaledVelocities = Sx * combinedVelocities * St; 

%We are now ready to compute the DMD modes their asscoiated eigenvalues.
X = scaledVelocities(:, 1:end-1);
Y = scaledVelocities(:, 2:end); %Same as "X-prime" in the literature.

%Compute a Singular Value Decomposition of the matrix X.
[U, Sigma, V] = svd(X, 'econ');
 
U= U(:, 1:r);
Sigma = Sigma(1:r, 1:r);  
V = V(:, 1:r);

Atilde = U' * Y * V * inv(Sigma);
[W, eigenvalueMat] = eig(Atilde); %W is a matrix of eigenvectors of Atilde.

%Initialize modeMat. 
modeMat = zeros(2*nx*ny, r); 

for index = 1:r
    modeMat(:,index) = (1/eigenvalueMat(index,index)) * Y * V * inv(Sigma) * W(:,index);
end 


%We shall now explicitly reformulate the modes as vector fields.
uModes = modeMat(1:nx*ny, :);
uModeMat = reshape(uModes, [ny, nx, r]);

vModes = modeMat((nx*ny + 1):end, :);
vModeMat = reshape(vModes, [ny, nx, r]);


end 

