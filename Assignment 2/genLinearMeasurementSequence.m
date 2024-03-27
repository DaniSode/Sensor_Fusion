function Y = genLinearMeasurementSequence(X, H, R)
n = size(X,1); % State dimension
N = size(X,2)-1; % Number of time steps
m = size(H,1); % Measurement dimension
Y = zeros(m,N);
for k = 1:N
    x = X(:,k+1); % Current state vector
    v = mvnrnd(zeros(m,1),R)'; % Measurement noise
    Y(:,k) = H*x + v; % Linear measurement model
end
end