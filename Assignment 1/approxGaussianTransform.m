function [mu_y, Sigma_y, y_s] = approxGaussianTransform(mu_x, Sigma_x, f, N)
%approxGaussianTransform takes a Gaussian density and a transformation 
%function and calculates the mean and covariance of the transformed density.
%
%Inputs
%   MU_X        [m x 1] Expected value of x.
%   SIGMA_X     [m x m] Covariance of x.
%   F           [Function handle] Function which maps a [m x 1] dimensional
%               vector into another vector of size [n x 1].
%   N           Number of samples to draw. Default = 5000.
%
%Output
%   MU_Y        [n x 1] Approximated mean of y.
%   SIGMA_Y     [n x n] Approximated covariance of y.
%   ys          [n x N] Samples propagated through f


if nargin < 4
    N = 5000;
end

%Your code here

% Define the length of the non-linear function with the given mu_x
n = length(f(mu_x));

% Draw N samples from the Gaussian density with parameters given in the input parameters and with the hint function: mvnrnd() 
x_s = mvnrnd(mu_x', Sigma_x, N)';

% Applying the non-linear function to each x in x_s, and concatenate the resulting y vectors into y_s
y_s = zeros(n, N);
for i = 1:N
    y_s(:, i) = f(x_s(:, i));
end

% Calculating the mean of each row and covariance of y_s
mu_y = mean(y_s, 2);
Sigma_y = cov(y_s');


end