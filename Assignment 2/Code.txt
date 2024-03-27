clear all; close all; clc

%% Scenario 1


%% First task

% Given 
sigmas = 3;
x0 = 2; 
P0 = 8;
Q = 1.5;  
R = 3;  
A = 1;  
H = 1;  
q_k_1 = normrnd(0, sqrt(Q));
r_k = normrnd(0, sqrt(R));  
N = 35; 

% Calculate the state and measurement sequence
X = genLinearStateSequence(x0,P0,A,Q,N); 
Y = genLinearMeasurementSequence(X,H,R);

% Plot results
figure;
plot(0:N, X, 'k-.', 'DisplayName', 'State','LineWidth', 2);
hold on;
grid on;
plot(1:N, Y, 'r*', 'DisplayName', 'Measurement');
xlabel('Time step');
ylabel('Value');
legend('show');



%% Second task a

% Using Kalman to estimate
[xhat, Phat] = kalmanFilter(Y, x0, P0, A, Q, H, R);

% Plot results
plot(0:N, [x0 xhat], 'b-', 'DisplayName', 'Estimate', 'LineWidth', 2);
sigma(1,:) = sqrt(Phat);
plot(0:N, [(x0 + sigmas*sqrt(P0)) (xhat + sigmas*sigma)], 'b--', 'DisplayName', 'Sigma 3+','LineWidth', 0.5);
plot(0:N, [(x0 - sigmas*sqrt(P0)) (xhat - sigmas*sigma)], 'b--', 'DisplayName', 'Sigma 3-','LineWidth', 0.5);



%% Second task b

% Define the different timesteps
k_indices = [1, 2, 4, 30];
error = zeros(length(k_indices),1);

% Plot with the error density at step 0
figure;
range = linspace(0-8, 0+8, 100);
density0 = normpdf(range, 0, sqrt(P0));
plot(range, density0, 'LineWidth', 1.5);
hold on;
grid on;

% Plot error densities at the different timesteps
for i = 1:length(k_indices)
    k = k_indices(i);
    error(i,1) = abs(X(k) - xhat(k));
    density = normpdf(range, 0, sigma(k));
    plot(range, density, 'LineWidth', 1.5);
end
labels = ["k = 1", "k = 2", "k = 4", "k = 30"];
legend("k=0", "k = 1", "k = 2", "k = 4", "k = 30");
xlabel('Error');
ylabel('Density');
title('Error densities at selected time steps');

% Write the error to command window
fprintf('The error at the selected time step:\n\n')
for j = 1:length(labels)
    fprintf('%s gives the error: %.5f\n', labels(j),error(j))
end



%% Third task

% Plot results
figure;
plot(0:N, X, 'k-.', 'DisplayName', 'State','LineWidth', 2);
hold on;
grid on;
plot(1:N, Y, 'r*', 'DisplayName', 'Measurement');
xlabel('Time step');
ylabel('Value');
legend('show');

% Using Kalman to estimate
plot(0:N, [x0 xhat], '-', 'color',[0 0.3 0.8], 'DisplayName', 'Correct Estimate', 'LineWidth', 2);
x0 = 12;
[xhat, Phat] = kalmanFilter(Y, x0, P0, A, Q, H, R);

% Plot results
plot(0:N, [x0 xhat], 'b-', 'DisplayName', 'Wrong Estimate', 'LineWidth', 2);
sigma(1,:) = sqrt(Phat);
plot(0:N, [(x0 + sigmas*sqrt(P0)) (xhat + sigmas*sigma)], 'b--', 'DisplayName', 'Sigma 3+','LineWidth', 0.5);
plot(0:N, [(x0 - sigmas*sqrt(P0)) (xhat - sigmas*sigma)], 'b--', 'DisplayName', 'Sigma 3-','LineWidth', 0.5);



%% Fourth task

% Get states and measurements
x0 = 2;
X = genLinearStateSequence(x0, P0, A, Q, N);
Y = genLinearMeasurementSequence(X, H, R);
N = size(Y,2);
xhat(:,1) = x0;
Phat(:,:,1) = P0 ;
for i=1:N
    % Prediction step
    [x_pred, P_pred] = linearPrediction(xhat(:,i), Phat(:,:,i), A, Q);
    
    % Update step
    [x_update, P_update] = linearUpdate(x_pred, P_pred, Y(:,i), H, R); 
    xhat(:,i+1) = x_update; % Save results
    Phat(:,:,i+1) = P_update; 
end

% Plot results
figure;
range = linspace(X(N)-6, X(N)+6, 100);
densityprior = normpdf(range, xhat(:,N), sqrt(Phat(:,:,N)));
densitypredict = normpdf(range, x_pred, sqrt(P_pred));
densityposterior = normpdf(range, xhat(:,N+1), sqrt(Phat(:,:,N+1)));
plot(range, densityprior, 'LineWidth', 1);
hold on;
grid on;
plot(range, densitypredict, 'LineWidth', 1);
plot(range, densityposterior, 'LineWidth', 1);
plot(ones(1,100)*Y(N),linspace(0,max(densityprior)*1.2), 'b--','linewidth',1.5)
plot(ones(1,100)*X(N),linspace(0,max(densityprior)*1.2), 'g--','linewidth',1.5)
legend('Prior', 'Predicted', 'Posterior', 'Measurement', 'State');
xlabel('States');
ylabel('Density');




%% Fifth task

% Generate state and measurement sequences
N = 1000;
amolag = 20;
X = genLinearStateSequence(x0, P0, A, Q, N);
Y = genLinearMeasurementSequence(X, H, R);

% Run Kalman filter
[xhat, Phat] = kalmanFilter(Y, x0, P0, A, Q, H, R);

% Plot results
screensize = get(0,'ScreenSize');
figure('Position',  [screensize(3)*0.35, screensize(4)*0.1, screensize(3)*0.3, screensize(4)*0.75]);
subplot(2,1,1);
histogram(X(1,2:end) - xhat(1,:), 'Normalization', 'pdf');
hold on;
grid on;
x = linspace(-6, 6, 100);
plot(x, normpdf(x, 0, Phat(end)), 'r', 'LineWidth', 2);
legend('Estimation error', 'Predicted distribution');
xlabel('x_k - xhat_k|k');
ylabel('pdf');

% Plot and calculate the autocorrelation of the innovation function along with the
subplot(2,1,2);
innovation = Y - H*xhat(:,1:end);
[acf, lags] = autocorr(innovation(1,:), amolag);
stem(lags, acf, 'filled');
hold on;
grid on;
plot(1:amolag, ones(1,amolag)*max(acf(1,2:end)),'r')
plot(1:amolag, ones(1,amolag)*min(acf(1,2:end)),'r')
ylim([min(acf(1,2:end))*3 1.05])
xlabel('Lag');
ylabel('Autocorrelation');



%% Scenario 2




%% First task

% Load data
load('SensorMeasurements.mat');
data1 = CalibrationSequenceVelocity_v0;
data2 = CalibrationSequenceVelocity_v10;
data3 = CalibrationSequenceVelocity_v20;

% Calculate the mean
mean1 = mean(data1);
mean2 = mean(data2);
mean3 = mean(data3);

% Calucalte the variance
var1 = var(data1);
var2 = var(data2);
var3 = var(data3);

% Calucalte total variance
Totvar = mean([var1 var2 var3]);

% Calculating the scaling factor C
Totmean = mean([mean2/10 mean3/20]);

% Print to command window
fprintf('\n\nThe total variance = %.5f\n\n',Totvar)
fprintf('The constant C = %.5f\n\n',Totmean)

% Plot result with the scaling factor
figure
plot(data1)
hold on
grid on
plot(data2)
plot(data3)
xlabel('Measurements')
ylabel('Speed')
figure
plot(data1)
hold on
grid on
plot(data2-10*Totmean)
plot(data3-20*Totmean)
xlabel('Measurements')
ylabel('Speed')
Chechvar = mean([var1 var(data2-10*Totmean) var(data3-20*Totmean)]);

% Write the difference between variance in both cases to command window
fprintf('The difference between the total variance before and after applying the scaling factor C is equal to: %.g\n\n',Chechvar-Totvar)



%% Second task

%% Constant velocity

% Define system parameters
x_0 = [0; 0]; 
P_0 = diag([Totvar, Totvar]); 
T = 0.1;
A = [1 T;   0 1]; 
Q = [T^4/4 T^3/2;
     T^3/2 T^2]; 
H = [1 0 ;0 1]; 
R = 1; 

% Generate data sequence from sensors
Y = Generate_y_seq();

% Use previous data for all NaN values
for i = 1:size(Y,2)
    if isnan(Y(1,i)) == 0
        val1 = Y(1,i);
    else 
        Y(1,i) = val1;
    end
    if isnan(Y(2,i)) == 0
        val2 = Y(2,i);
    else 
        Y(2,i) = val2;
    end
end

% Initialize Kalman filter
[xhat, Phat] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

% Plot results
screensize = get(0,'ScreenSize');
figure('Position',  [screensize(3)*0.35, screensize(4)*0.1, screensize(3)*0.3, screensize(4)*0.75]);
subplot(2,1,1)
plot(Y(1,:), 'r', 'MarkerSize', 2); 
hold on
grid on
plot(xhat(1,:), 'b', 'LineWidth', 2); 
legend('Noisy position measurements', 'Kalman filter estimate');
xlabel('Time');
ylabel('Position');
subplot(2,1,2)
plot(Y(2,:), 'r', 'MarkerSize', 2); 
hold on
grid on
plot(xhat(2,:), 'b', 'LineWidth', 2); 
legend('Noisy position measurements', 'Kalman filter estimate');
xlabel('Time');
ylabel('Velocity');
sgtitle('Constant velocity model')

% % Plot and calculate the autocorrelation of the innovation function along with the
% % Position
% N = 1000;
% amolag = 20;
% screensize = get(0,'ScreenSize');
% figure('Position',  [screensize(3)*0.35, screensize(4)*0.1, screensize(3)*0.3, screensize(4)*0.75]);
% subplot(2,1,1);
% innovation = Y(1,:) - H*xhat(1,:);
% [acf, lags] = autocorr(innovation(1,:), amolag);
% stem(lags, acf, 'filled');
% hold on;
% grid on;
% plot(1:amolag, ones(1,amolag)*0.05,'r')
% plot(1:amolag, ones(1,amolag)*-0.05,'r')
% ylim([-0.1 1.05])
% xlabel('Lag');
% 
% % Velocity
% subplot(2,1,2);
% innovation = Y(2,:) - H*xhat(2,:);
% [acf, lags] = autocorr(innovation(1,:), amolag);
% stem(lags, acf, 'filled');
% hold on;
% grid on;
% plot(1:amolag, ones(1,amolag)*0.05,'r')
% plot(1:amolag, ones(1,amolag)*-0.05,'r')
% ylim([-0.1 1.05])
% xlabel('Lag');
% ylabel('Autocorrelation');



%% Constant acceleration

% Define system parameters
x_0 = [0; 0; 0]; 
P_0 = diag([Totvar, Totvar, Totvar]); 
T = 0.1;
A = [1 T T^2/2;
     0 1 T;
     0 0 1]; 
Q = [T^5/20 T^4/8 T^3/6;
     T^4/8 T^3/3 T^2/2;
     T^3/6 T^2/2 T]; 
H = [1 1 0;
     0 1 1]; 
R = 1; 

% Generate data sequence from sensors
Y = Generate_y_seq();

for i = 1:size(Y,2)
    if isnan(Y(1,i)) == 0
        val1 = Y(1,i);
    else 
        Y(1,i) = val1;
    end
    if isnan(Y(2,i)) == 0
        val2 = Y(2,i);
    else 
        Y(2,i) = val2;
    end
end

% Initialize Kalman filter
[xhat, Phat] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

% Plot results
screensize = get(0,'ScreenSize');
figure('Position',  [screensize(3)*0.35, screensize(4)*0.1, screensize(3)*0.3, screensize(4)*0.75]);
subplot(2,1,1)
plot(Y(1,:), 'r', 'MarkerSize', 2); % noisy position measurements
hold on
grid on
plot(xhat(1,:), 'b', 'LineWidth', 2); % Kalman filter estimate
legend('Noisy position measurements', 'Kalman filter estimate');
xlabel('Time');
ylabel('Position');
subplot(2,1,2)
plot(Y(2,:), 'r', 'MarkerSize', 2); % noisy position measurements
hold on
grid on
plot(xhat(2,:), 'b', 'LineWidth', 2); % Kalman filter estimate
legend('Noisy position measurements', 'Kalman filter estimate');
xlabel('Time');
ylabel('Velocity');
sgtitle('Constant acceleration model')

% % Plot and calculate the autocorrelation of the innovation function along with the
% % Position
% N = 1000;
% amolag = 20;
% screensize = get(0,'ScreenSize');
% figure('Position',  [screensize(3)*0.35, screensize(4)*0.1, screensize(3)*0.3, screensize(4)*0.75]);
% subplot(2,1,1);
% innovation = Y(1,:) - H(:,1)*xhat(1,:);
% [acf, lags] = autocorr(innovation(1,:), amolag);
% stem(lags, acf, 'filled');
% hold on;
% grid on;
% plot(1:amolag, ones(1,amolag)*0.05,'r')
% plot(1:amolag, ones(1,amolag)*-0.05,'r')
% ylim([-0.1 1.05])
% xlabel('Lag');
% 
% % Velocity
% subplot(2,1,2);
% innovation = Y(2,:) - H(:,2)*xhat(2,:);
% [acf, lags] = autocorr(innovation(1,:), amolag);
% stem(lags, acf, 'filled');
% hold on;
% grid on;
% plot(1:amolag, ones(1,amolag)*0.05,'r')
% plot(1:amolag, ones(1,amolag)*-0.05,'r')
% ylim([min(acf(1,2:end))*3 1.05])
% xlabel('Lag');
% ylabel('Autocorrelation');


%% Functions

% Get linear state function
function X = genLinearStateSequence(x_0, P_0, A, Q, N)
d = length(x_0); % dimension of the state
X = zeros(d, N+1); % matrix to store the state vectors
X(:,1) = mvnrnd(x_0, P_0)'; % transpose to obtain a row vector
for k = 2:N+1
    X(:,k) = A*X(:,k-1) + mvnrnd(zeros(d,1), Q)'; % transpose to obtain a row vector
end
end

% Get measurement
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

% Prediction function
function [x, P] = linearPrediction(x, P, A, Q)
x = A * x;
P = A * P * A' + Q;
end

% Update state function
function [x, P] = linearUpdate(x, P, y, H, R)
K = P * H' / (H * P * H' + R);
x = x + K * (y - H * x);
P = (eye(size(P)) - K * H) * P;
end

% Kalmanfilter
function [X, P] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
N = size(Y,2);
n = length(x_0);
m = size(Y,1);
x = zeros(n,N);
P = zeros(n,n,N);
X(:,1) = x_0;
P(:,:,1) = P_0 ;
for i=1:N
    % Prediction step
    [x_pred, P_pred] = linearPrediction(X(:,i), P(:,:,i), A, Q);
    
    % Update step
    [x_update, P_update] = linearUpdate(x_pred, P_pred, Y(:,i), H, R); 
    X(:,i+1) = x_update; % Save results
    P(:,:,i+1) = P_update; 
end
X = X(:,2:end); 
P = P(:,:,2:end);
end