clear all; close all; clc


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
fprintf('The total variance = %.5f\n\n',Totvar)
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