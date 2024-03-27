function [ xy ] = sigmaEllipse2D( mu, Sigma, level, npoints )
%SIGMAELLIPSE2D generates x,y-points which lie on the ellipse describing
% a sigma level in the Gaussian density defined by mean and covariance.
%
%Input:
%   MU          [2 x 1] Mean of the Gaussian density
%   SIGMA       [2 x 2] Covariance matrix of the Gaussian density
%   LEVEL       Which sigma level curve to plot. Can take any positive value, 
%               but common choices are 1, 2 or 3. Default = 3.
%   NPOINTS     Number of points on the ellipse to generate. Default = 32.
%
%Output:
%   XY          [2 x npoints] matrix. First row holds x-coordinates, second
%               row holds the y-coordinates. First and last columns should 
%               be the same point, to create a closed curve.


%Setting default values, in case only mu and Sigma are specified.
if nargin < 3
    level = 3;
end
if nargin < 4
    npoints = 32;
end

% Compute A as the square root of the elements in the covariance matrix Sigma
A = sqrtm(Sigma);

% Set the angle range and generate npoints angles as described above
angles = linspace(0, 2*pi, npoints);

% Compute the points on the ellipse for each angle using the level input to control the radius condition
points = level.* [cos(angles); sin(angles)];

% Rotating and translating the calculated "points" to the ellipse to fit the covariance matrix A and mean mu
xy = A * points + mu;

end