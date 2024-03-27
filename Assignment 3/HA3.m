%%%%%%%%%%%%
%% Task 1 %%
%%%%%%%%%%%%



%% Part a) %%

% Clean slate
clear all; close all; clc

% Define number of samples
N = 1000;

% Define state densities with 5 states to be able to use the motionmodel
x0_1 = [125; 125];
P0_1 = [10^2 0;0 5^2];
x0_2 = [-25; 125]; 
P0_2 = [10^2 0;0 5^2];
x0_3 = [60; 60];
P0_3 = [10^2 0;0 5^2];

% Define sensor positions and measurement noise
s1 = [0; 100];
s2 = [100; 0];
sigma_phi = 0.1*pi/180; % To get radians
R = [sigma_phi^2 0;0 sigma_phi^2]; % Just 2x2 since we just care about the 2 first states

% The corresponding process noise covariance matrix is the same for all the
% cases here since the velocity v and angular velocity omega are 0. 
Gamma = [0 0;0 0;1 0;0 0;0 1];
P_v = 0;
P_w = 0;
% Q = Gamma*[P_v^2 0;0 P_w^2]*Gamma';
Q = 0; % In our case its easier to set the process noise to 0

% Define timestep (not given therefore just put it to 1)
T = 1;

% Defining function handles for needed functions
f = @(x) coordinatedTurnMotion(x, T);
h = @(x) dualBearingMeasurement(x, s1, s2);
y = @(x) genNonLinearMeasurementSequence(x, h, R);

%% Calculate the mean and covariance and printing results
[samp_mean1, samp_cov1, Y1] = approxGaussianTransform(x0_1, P0_1, y, N);
[samp_mean2, samp_cov2, Y2] = approxGaussianTransform(x0_2, P0_2, y, N);
[samp_mean3, samp_cov3, Y3] = approxGaussianTransform(x0_3, P0_3, y, N);
fprintf('The mean of the samples for case 1:\n\n'); disp(samp_mean1)
fprintf('The mean of the samples for case 2:\n\n'); disp(samp_mean2)
fprintf('The mean of the samples for case 3:\n\n'); disp(samp_mean3)
fprintf('The covariance of the samples for case 1:\n\n'); disp(samp_cov1)
fprintf('The covariance of the samples for case 2:\n\n'); disp(samp_cov2)
fprintf('The covariance of the samples for case 3:\n\n'); disp(samp_cov3)



%% Part b) %%

%% EKF simulation and printing results
[EKF_mean1, EKF_cov1] = nonLinKFprediction(x0_1, P0_1, h, Q, 'EKF');
[EKF_mean2, EKF_cov2] = nonLinKFprediction(x0_2, P0_2, h, Q, 'EKF');
[EKF_mean3, EKF_cov3] = nonLinKFprediction(x0_3, P0_3, h, Q, 'EKF');
fprintf('The mean using EKF for case 1:\n\n'); disp(EKF_mean1)
fprintf('The mean using EKF for case 2:\n\n'); disp(EKF_mean2)
fprintf('The mean using EKF for case 3:\n\n'); disp(EKF_mean3)
fprintf('The covariance using EKF for case 1:\n\n'); disp(EKF_cov1)
fprintf('The covariance using EKF for case 2:\n\n'); disp(EKF_cov2)
fprintf('The covariance using EKF for case 3:\n\n'); disp(EKF_cov3)

%% UKF simulation and printing results
[UKF_mean1, UKF_cov1] = nonLinKFprediction(x0_1, P0_1, h, Q, 'UKF');
[UKF_mean2, UKF_cov2] = nonLinKFprediction(x0_2, P0_2, h, Q, 'UKF');
[UKF_mean3, UKF_cov3] = nonLinKFprediction(x0_3, P0_3, h, Q, 'UKF');
fprintf('The mean using UKF for case 1:\n\n'); disp(UKF_mean1)
fprintf('The mean using UKF for case 2:\n\n'); disp(UKF_mean2)
fprintf('The mean using UKF for case 3:\n\n'); disp(UKF_mean3)
fprintf('The covariance using UKF for case 1:\n\n'); disp(UKF_cov1)
fprintf('The covariance using UKF for case 2:\n\n'); disp(UKF_cov2)
fprintf('The covariance using UKF for case 3:\n\n'); disp(UKF_cov3)

%% CKF simulation and printing results
[CKF_mean1, CKF_cov1] = nonLinKFprediction(x0_1, P0_1, h, Q, 'CKF');
[CKF_mean2, CKF_cov2] = nonLinKFprediction(x0_2, P0_2, h, Q, 'CKF');
[CKF_mean3, CKF_cov3] = nonLinKFprediction(x0_3, P0_3, h, Q, 'CKF');
fprintf('The mean using CKF for case 1:\n\n'); disp(CKF_mean1)
fprintf('The mean using CKF for case 2:\n\n'); disp(CKF_mean2)
fprintf('The mean using CKF for case 3:\n\n'); disp(CKF_mean3)
fprintf('The covariance using CKF for case 1:\n\n'); disp(CKF_cov1)
fprintf('The covariance using CKF for case 2:\n\n'); disp(CKF_cov2)
fprintf('The covariance using CKF for case 3:\n\n'); disp(CKF_cov3)



%% Part c) %%

% Define level
level = 3;

% Case 1
[xy1] = sigmaEllipse2D(samp_mean1, samp_cov1, level, N);

% Sample
size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.85]); 
subplot(3,3,1)
scatter(Y1(1,:),Y1(2,:),10,'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880])
title('Sample vs EKF (CASE 1)'); xlabel('x'); ylabel('y')
hold on; grid on;
scatter(samp_mean1(1),samp_mean1(2),'b*')
plot(xy1(1,:),xy1(2,:),'b')

% EKF
[xyEKF1] = sigmaEllipse2D(EKF_mean1, EKF_cov1, level, N);
scatter(EKF_mean1(1), EKF_mean1(2),'m*')
plot(xyEKF1(1,:),xyEKF1(2,:),'m')
legend('Samples','Mean','Covariance','EKF mean','EKF covariance','FontSize',5,'location','best')

% Sample
subplot(3,3,2)
scatter(Y1(1,:),Y1(2,:),10,'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880])
title('Sample vs UKF (CASE 1)'); xlabel('x'); ylabel('y')
hold on; grid on;
scatter(samp_mean1(1),samp_mean1(2),'b*')
plot(xy1(1,:),xy1(2,:),'b')

% UKF
[xyUKF1] = sigmaEllipse2D(UKF_mean1, UKF_cov1, level, N);
scatter(UKF_mean1(1), UKF_mean1(2),'m*')
plot(xyUKF1(1,:),xyUKF1(2,:),'m')
legend('Samples','Mean','Covariance','UKF mean','UKF covariance','FontSize',5,'location','best')

% Sample
subplot(3,3,3)
scatter(Y1(1,:),Y1(2,:),10,'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880])
title('Sample vs CKF (CASE 1)'); xlabel('x'); ylabel('y')
hold on; grid on;
scatter(samp_mean1(1),samp_mean1(2),'b*')
plot(xy1(1,:),xy1(2,:),'b')

% CKF
[xyCKF1] = sigmaEllipse2D(CKF_mean1, CKF_cov1, level, N);
scatter(CKF_mean1(1), CKF_mean1(2),'m*')
plot(xyCKF1(1,:),xyCKF1(2,:),'m')
legend('Samples','Mean','Covariance','CKF mean','CKF covariance','FontSize',5,'location','best')


% Case 2
[xy2] = sigmaEllipse2D(samp_mean2, samp_cov2, level, N);

% Sample
subplot(3,3,4)
scatter(Y2(1,:),Y2(2,:),10,'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880])
title('Sample vs EKF (CASE 2)'); xlabel('x'); ylabel('y')
hold on; grid on;
scatter(samp_mean2(1),samp_mean2(2),'b*')
plot(xy2(1,:),xy2(2,:),'b')

% EKF
[xyEKF2] = sigmaEllipse2D(EKF_mean2, EKF_cov2, level, N);
scatter(EKF_mean2(1), EKF_mean2(2),'m*')
plot(xyEKF2(1,:),xyEKF2(2,:),'m')
legend('Samples','Mean','Covariance','EKF mean','EKF covariance','FontSize',5,'location','best')

% Sample
subplot(3,3,5)
scatter(Y2(1,:),Y2(2,:),10,'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880])
title('Sample vs UKF (CASE 2)'); xlabel('x'); ylabel('y')
hold on; grid on;
scatter(samp_mean2(1),samp_mean2(2),'b*')
plot(xy2(1,:),xy2(2,:),'b')

% UKF
[xyUKF2] = sigmaEllipse2D(UKF_mean2, UKF_cov2, level, N);
scatter(UKF_mean2(1), UKF_mean2(2),'m*')
plot(xyUKF2(1,:),xyUKF2(2,:),'m')
legend('Samples','Mean','Covariance','UKF mean','UKF covariance','FontSize',5,'location','best')

% Sample
subplot(3,3,6)
scatter(Y2(1,:),Y2(2,:),10,'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880])
title('Sample vs CKF (CASE 2)'); xlabel('x'); ylabel('y')
hold on; grid on;
scatter(samp_mean2(1),samp_mean2(2),'b*')
plot(xy2(1,:),xy2(2,:),'b')

% CKF
[xyCKF2] = sigmaEllipse2D(CKF_mean2, CKF_cov2, level, N);
scatter(CKF_mean2(1), CKF_mean2(2),'m*')
plot(xyCKF2(1,:),xyCKF2(2,:),'m')
legend('Samples','Mean','Covariance','CKF mean','CKF covariance','FontSize',5,'location','best')


% Case 3
[xy3] = sigmaEllipse2D(samp_mean3, samp_cov3, level, N);

% Sample
subplot(3,3,7)
scatter(Y3(1,:),Y3(2,:),10,'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880])
title('Sample vs EKF (CASE 3)'); xlabel('x'); ylabel('y')
hold on; grid on;
scatter(samp_mean3(1),samp_mean3(2),'b*')
plot(xy3(1,:),xy3(2,:),'b')

%EKF
[xyEKF3] = sigmaEllipse2D(EKF_mean3, EKF_cov3, level, N);
scatter(EKF_mean3(1), EKF_mean3(2),'m*')
plot(xyEKF3(1,:),xyEKF3(2,:),'m')
legend('Samples','Mean','Covariance','EKF mean','EKF covariance','FontSize',5,'location','best')

% Sample
subplot(3,3,8)
scatter(Y3(1,:),Y3(2,:),10,'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880])
title('Sample vs UKF (CASE 3)'); xlabel('x'); ylabel('y')
hold on; grid on;
scatter(samp_mean3(1),samp_mean3(2),'b*')
plot(xy3(1,:),xy3(2,:),'b')

% UKF
[xyUKF3] = sigmaEllipse2D(UKF_mean3, UKF_cov3, level, N);
scatter(UKF_mean3(1), UKF_mean3(2),'m*')
plot(xyUKF3(1,:),xyUKF3(2,:),'m')
legend('Samples','Mean','Covariance','UKF mean','UKF covariance','FontSize',5,'location','best')

% Sample
subplot(3,3,9)
scatter(Y3(1,:),Y3(2,:),10,'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880])
title('Sample vs CKF (CASE 3)'); xlabel('x'); ylabel('y')
hold on; grid on;
scatter(samp_mean3(1),samp_mean3(2),'b*')
plot(xy3(1,:),xy3(2,:),'b')

% CKF
[xyCKF3] = sigmaEllipse2D(CKF_mean3, CKF_cov3, level, N);
scatter(CKF_mean3(1), CKF_mean3(2),'m*')
plot(xyCKF3(1,:),xyCKF3(2,:),'m')
legend('Samples','Mean','Covariance','CKF mean','CKF covariance','FontSize',5,'location','best')



%%%%%%%%%%%%
%% Task 2 %%
%%%%%%%%%%%%

%% Part a) %%
clear all; close all; clc

% Define level
level = 3;

% Define prior
x0 = [0; 0; 20; 0; 5*pi/180];
P0 = diag([10^2 10^2 2^2 (pi/180)^2 (pi/180)^2]);

% Sample time and amount of samples
T = 1;
N = 100;

% The corresponding process noise covariance matrix is the same for all the
% cases here since the velocity v and angular velocity omega are 0. 
Gamma = [0 0;0 0;1 0;0 0;0 1];
P_v = 1;
P_w = pi/180;
Q = Gamma*[T*P_v^2 0;0 T*P_w^2]*Gamma';

% Define sensors
s1 = [-200; 100];
s2 = [-200; -100];

% Different noise depending on different cases
case1_phi1 = 2*pi/180;
case1_phi2 = 2*pi/180;
case2_phi1 = 2*pi/180;
case2_phi2 = 0.1*pi/180;
case3_phi1 = 0.1*pi/180;
case3_phi2 = 0.1*pi/180;
R1 = diag([case1_phi1 case1_phi2]);
R2 = diag([case2_phi1 case2_phi2]);
R3 = diag([case3_phi1 case3_phi2]);

% Defining function handles for needed functions
f = @(x) coordinatedTurnMotion(x, T);
h = @(x) dualBearingMeasurement(x, s1, s2);

%% Genereate state
X = genNonLinearStateSequence(x0, P0, f, Q, N);

% Generate measurements
Y1 = genNonLinearMeasurementSequence(X, h, R1);
Y2 = genNonLinearMeasurementSequence(X, h, R2);
Y3 = genNonLinearMeasurementSequence(X, h, R3);

x1 = [];
y1 = [];
x2 = [];
y2 = [];
x3 = [];
y3 = [];

for i = 1:length(Y1)
     [x1(i), y1(i)] = getPosFromMeasurement(Y1(1,i), Y1(2,i), s1, s2);
     [x2(i), y2(i)] = getPosFromMeasurement(Y2(1,i), Y2(2,i), s1, s2);
     [x3(i), y3(i)] = getPosFromMeasurement(Y3(1,i), Y3(2,i), s1, s2);
end



%% Part b) %%

%% EKF
[xf_EKF_1, Pf_EKF_1, xp_EKF_1, Pp_EKF_1] = nonLinearKalmanFilter(Y1, x0, P0, f, Q, h, R1, 'EKF');
[xf_EKF_2, Pf_EKF_2, xp_EKF_2, Pp_EKF_2] = nonLinearKalmanFilter(Y2, x0, P0, f, Q, h, R2, 'EKF');
[xf_EKF_3, Pf_EKF_3, xp_EKF_3, Pp_EKF_3] = nonLinearKalmanFilter(Y3, x0, P0, f, Q, h, R3, 'EKF');

%% UKF
[xf_UKF_1, Pf_UKF_1, xp_UKF_1, Pp_UKF_1] = nonLinearKalmanFilter(Y1, x0, P0, f, Q, h, R1, 'UKF');
[xf_UKF_2, Pf_UKF_2, xp_UKF_2, Pp_UKF_2] = nonLinearKalmanFilter(Y2, x0, P0, f, Q, h, R2, 'UKF');
[xf_UKF_3, Pf_UKF_3, xp_UKF_3, Pp_UKF_3] = nonLinearKalmanFilter(Y3, x0, P0, f, Q, h, R3, 'UKF');

%% CKF
[xf_CKF_1, Pf_CKF_1, xp_CKF_1, Pp_CKF_1] = nonLinearKalmanFilter(Y1, x0, P0, f, Q, h, R1, 'CKF');
[xf_CKF_2, Pf_CKF_2, xp_CKF_2, Pp_CKF_2] = nonLinearKalmanFilter(Y2, x0, P0, f, Q, h, R2, 'CKF');
[xf_CKF_3, Pf_CKF_3, xp_CKF_3, Pp_CKF_3] = nonLinearKalmanFilter(Y3, x0, P0, f, Q, h, R3, 'CKF');


%% Plot results
axlim1 = [-350 550 -150 1000];
axlim2 = [-200 400 -150 650];
axlim3 = [-200 380 -150 450];

% Plot true trajectory and sensor positions
size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.85]); 
subplot(3,3,1)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('True vs EKF (CASE 1)'); xlabel('x'); ylabel('y')
hold on;
scatter(x1, y1, 'r.');
plot(xf_EKF_1(1,:), xf_EKF_1(2,:), 'm-')
for i = 1:5:N
    xy = sigmaEllipse2D(xf_EKF_1(1:2, i), Pf_EKF_1(1:2,1:2, i), level, N);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim1);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'EKF preduction', 'EKF covariance', 'location', 'best', 'FontSize', 5);

subplot(3,3,2)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('True vs UKF (CASE 1)'); xlabel('x'); ylabel('y')
hold on;
scatter(x1, y1, 'r.');
plot(xf_UKF_1(1,:), xf_UKF_1(2,:), 'm-')
for i = 1:5:N
    xy = sigmaEllipse2D(xf_UKF_1(1:2, i), Pf_UKF_1(1:2,1:2, i), level, N);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim1);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'UKF preduction', 'UKF covariance', 'location', 'best', 'FontSize', 5);

subplot(3,3,3)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('True vs CKF (CASE 1)'); xlabel('x'); ylabel('y')
hold on;
scatter(x1, y1, 'r.');
plot(xf_CKF_1(1,:), xf_CKF_1(2,:), 'm-')
for i = 1:5:N
    xy = sigmaEllipse2D(xf_CKF_1(1:2, i), Pf_CKF_1(1:2,1:2, i), level, N);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim1);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'CKF preduction', 'CKF covariance', 'location', 'best', 'FontSize', 5);

subplot(3,3,4)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('True vs EKF (CASE 2)'); xlabel('x'); ylabel('y')
hold on;
scatter(x2, y2, 'r.')
plot(xf_EKF_2(1,:), xf_EKF_2(2,:), 'm-')
for i = 1:5:N
    xy = sigmaEllipse2D(xf_EKF_2(1:2, i), Pf_EKF_2(1:2,1:2, i), level, N);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim2);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'EKF preduction', 'EKF covariance', 'location', 'best', 'FontSize', 5);

subplot(3,3,5)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('True vs UKF (CASE 2)'); xlabel('x'); ylabel('y')
hold on;
scatter(x2, y2, 'r.')
plot(xf_UKF_2(1,:), xf_UKF_2(2,:), 'm-')
for i = 1:5:N
    xy = sigmaEllipse2D(xf_UKF_2(1:2, i), Pf_UKF_2(1:2,1:2, i), level, N);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim2);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'UKF preduction', 'UKF covariance', 'location', 'best', 'FontSize', 5);

subplot(3,3,6)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('True vs CKF (CASE 2)'); xlabel('x'); ylabel('y')
hold on;
scatter(x2, y2, 'r.')
plot(xf_CKF_2(1,:), xf_CKF_2(2,:), 'm-')
for i = 1:5:N
    xy = sigmaEllipse2D(xf_CKF_2(1:2, i), Pf_CKF_2(1:2,1:2, i), level, N);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim2);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'CKF preduction', 'CKF covariance', 'location', 'best', 'FontSize', 5);

subplot(3,3,7)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('True vs EKF (CASE 3)'); xlabel('x'); ylabel('y')
hold on;
scatter(x3, y3, 'r.')
plot(xf_EKF_3(1,:), xf_EKF_3(2,:), 'm-')
for i = 1:5:N
    xy = sigmaEllipse2D(xf_EKF_3(1:2, i), Pf_EKF_3(1:2,1:2, i), level, N);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim3);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'EKF preduction', 'EKF covariance', 'location', 'best', 'FontSize', 5);

subplot(3,3,8)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('True vs UKF (CASE 3)'); xlabel('x'); ylabel('y')
hold on;
scatter(x3, y3, 'r.')
plot(xf_UKF_3(1,:), xf_UKF_3(2,:), 'm-')
for i = 1:5:N
    xy = sigmaEllipse2D(xf_UKF_3(1:2, i), Pf_UKF_3(1:2,1:2, i), level, N);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim3);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'UKF preduction', 'UKF covariance', 'location', 'best', 'FontSize', 5);

subplot(3,3,9)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('True vs CKF (CASE 3)'); xlabel('x'); ylabel('y')
hold on;
scatter(x3, y3, 'r.')
plot(xf_CKF_3(1,:), xf_CKF_3(2,:), 'm-')
for i = 1:5:N
    xy = sigmaEllipse2D(xf_CKF_3(1:2, i), Pf_CKF_3(1:2,1:2, i), level, N);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim3);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'CKF preduction', 'CKF covariance', 'location', 'best', 'FontSize', 5);



%% Part c %% 

size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.85]); 

N = 10;
MC = 10;
error = [];
for imc = 1:MC
% Simulate state sequence
X = genNonLinearStateSequence(x0, P0, f, Q, N);
% Simulate measurements
Y = genNonLinearMeasurementSequence(X, h, R1);
% Run Kalman filter (you need to run all three, for comparison)
[xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x0,P0,f,Q,h,R1,'EKF');
% Save the estimation errors and the prediction errors!
error = [error, [xf(1,:)-X(1,2:end); xf(2,:)-X(2,2:end)]];
end
subplot(6,3,1)
histogram(error(1,:),'Normalization','pdf')
title('delta x EKF case 1') 
subplot(6,3,4)
histogram(error(2,:),'Normalization','pdf')
title('delta y EKF case 1')

error = [];
for imc = 1:MC
% Simulate state sequence
X = genNonLinearStateSequence(x0, P0, f, Q, N);
% Simulate measurements
Y = genNonLinearMeasurementSequence(X, h, R1);
% Run Kalman filter (you need to run all three, for comparison)
[xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x0,P0,f,Q,h,R1,'UKF');
% Save the estimation errors and the prediction errors!
error = [error, [xf(1,:)-X(1,2:end); xf(2,:)-X(2,2:end)]];
end
subplot(6,3,2)
histogram(error(1,:),'Normalization','pdf')
title('delta x UKF case 1') 
subplot(6,3,5)
histogram(error(2,:),'Normalization','pdf')
title('delta y UKF case 1')

error = [];
for imc = 1:MC
% Simulate state sequence
X = genNonLinearStateSequence(x0, P0, f, Q, N);
% Simulate measurements
Y = genNonLinearMeasurementSequence(X, h, R1);
% Run Kalman filter (you need to run all three, for comparison)
[xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x0,P0,f,Q,h,R1,'CKF');
% Save the estimation errors and the prediction errors!
error = [error, [xf(1,:)-X(1,2:end); xf(2,:)-X(2,2:end)]];
end
subplot(6,3,3)
histogram(error(1,:),'Normalization','pdf')
title('delta x CKF case 1') 
subplot(6,3,6)
histogram(error(2,:),'Normalization','pdf')
title('delta y CKF case 1')

error = [];
for imc = 1:MC
% Simulate state sequence
X = genNonLinearStateSequence(x0, P0, f, Q, N);
% Simulate measurements
Y = genNonLinearMeasurementSequence(X, h, R2);
% Run Kalman filter (you need to run all three, for comparison)
[xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x0,P0,f,Q,h,R2,'EKF');
% Save the estimation errors and the prediction errors!
error = [error, [xf(1,:)-X(1,2:end); xf(2,:)-X(2,2:end)]];
end
subplot(6,3,7)
histogram(error(1,:),'Normalization','pdf')
title('delta x EKF case 2') 
subplot(6,3,10)
histogram(error(2,:),'Normalization','pdf')
title('delta y EKF case 2')

error = [];
for imc = 1:MC
% Simulate state sequence
X = genNonLinearStateSequence(x0, P0, f, Q, N);
% Simulate measurements
Y = genNonLinearMeasurementSequence(X, h, R2);
% Run Kalman filter (you need to run all three, for comparison)
[xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x0,P0,f,Q,h,R2,'UKF');
% Save the estimation errors and the prediction errors!
error = [error, [xf(1,:)-X(1,2:end); xf(2,:)-X(2,2:end)]];
end
subplot(6,3,8)
histogram(error(1,:),'Normalization','pdf')
title('delta x UKF case 2') 
subplot(6,3,11)
histogram(error(2,:),'Normalization','pdf')
title('delta y UKF case 2')

error = [];
for imc = 1:MC
% Simulate state sequence
X = genNonLinearStateSequence(x0, P0, f, Q, N);
% Simulate measurements
Y = genNonLinearMeasurementSequence(X, h, R2);
% Run Kalman filter (you need to run all three, for comparison)
[xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x0,P0,f,Q,h,R2,'CKF');
% Save the estimation errors and the prediction errors!
error = [error, [xf(1,:)-X(1,2:end); xf(2,:)-X(2,2:end)]];
end
subplot(6,3,9)
histogram(error(1,:),'Normalization','pdf')
title('delta x CKF case 2') 
subplot(6,3,12)
histogram(error(2,:),'Normalization','pdf')
title('delta y CKF case 2')

error = [];
for imc = 1:MC
% Simulate state sequence
X = genNonLinearStateSequence(x0, P0, f, Q, N);
% Simulate measurements
Y = genNonLinearMeasurementSequence(X, h, R3);
% Run Kalman filter (you need to run all three, for comparison)
[xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x0,P0,f,Q,h,R3,'EKF');
% Save the estimation errors and the prediction errors!
error = [error, [xf(1,:)-X(1,2:end); xf(2,:)-X(2,2:end)]];
end
subplot(6,3,13)
histogram(error(1,:),'Normalization','pdf')
title('delta x EKF case 3') 
subplot(6,3,16)
histogram(error(2,:),'Normalization','pdf')
title('delta y EKF case 3')

error = [];
for imc = 1:MC
% Simulate state sequence
X = genNonLinearStateSequence(x0, P0, f, Q, N);
% Simulate measurements
Y = genNonLinearMeasurementSequence(X, h, R3);
% Run Kalman filter (you need to run all three, for comparison)
[xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x0,P0,f,Q,h,R3,'UKF');
% Save the estimation errors and the prediction errors!
error = [error, [xf(1,:)-X(1,2:end); xf(2,:)-X(2,2:end)]];
end
subplot(6,3,14)
histogram(error(1,:),'Normalization','pdf')
title('delta x UKF case 3') 
subplot(6,3,17)
histogram(error(2,:),'Normalization','pdf')
title('delta y UKF case 3')

error = [];
for imc = 1:MC
% Simulate state sequence
X = genNonLinearStateSequence(x0, P0, f, Q, N);
% Simulate measurements
Y = genNonLinearMeasurementSequence(X, h, R3);
% Run Kalman filter (you need to run all three, for comparison)
[xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x0,P0,f,Q,h,R3,'CKF');
% Save the estimation errors and the prediction errors!
error = [error, [xf(1,:)-X(1,2:end); xf(2,:)-X(2,2:end)]];
end
subplot(6,3,15)
histogram(error(1,:),'Normalization','pdf')
title('delta x CKF case 3') 
subplot(6,3,18)
histogram(error(2,:),'Normalization','pdf')
title('delta y CKF case 3')



%%%%%%%%%%%%
%% Task 3 %%
%%%%%%%%%%%%

%% Part a) %% 
clear all; close all; clc 

% Level
level = 3;

% True track
% Sampling period
T = 0.1;
% Length of time sequence
K = 60;
% Allocate memory
omega = zeros(1,K+1);
% Turn rate
omega(15:45) = -pi/301/T;
% Initial state
x0 = [0 0 20 0 omega(1)]';
% Allocate memory
X = zeros(length(x0),K+1);
X(:,1) = x0;
% Create true track
for i=2:K+1
% Simulate
X(:,i) = coordinatedTurnMotion(X(:,i-1), T);
% Set turn−rate
X(5,i) = omega(i);
end

% Define prior
x0 = [0; 0; 20; 0; 5*pi/180];
P0 = diag([10^2 10^2 2^2 (5*pi/180)^2 (pi/180)^2]);

% The corresponding process noise covariance matrix is the same for all the
% cases here since the velocity v and angular velocity omega are 0. 
Gamma = [0 0;0 0;1 0;0 0;0 1];
P_v = [1 0.001 1000];
P_w = [pi/180 0.001*pi/180 1000*pi/180];
P_phi1 = pi/180;
P_phi2 = pi/180;
R = diag([P_phi1 P_phi2]);
Q1 = Gamma*[T*P_v(1)^2 0;0 T*P_w(1)^2]*Gamma'; % Start also well tuned
Q2 = Gamma*[T*P_v(3)^2 0;0 T*P_w(3)^2]*Gamma'; % both large
Q3 = Gamma*[T*P_v(2)^2 0;0 T*P_w(2)^2]*Gamma'; % both small

% Define sensors
s1 = [300; -100];
s2 = [300; -300];

% Defining function handles for needed functions
f = @(x) coordinatedTurnMotion(x, T);
h = @(x) dualBearingMeasurement(x, s1, s2);

%% Generate measurements
Y = genNonLinearMeasurementSequence(X, h, R);

% Calc cart cord
for i = 1:length(Y)
     [x1(i), y1(i)] = getPosFromMeasurement(Y(1,i), Y(2,i), s1, s2);
end

% Generate filtered sequences for all cases
[xf1, Pf1, xp1, Pp1] = nonLinearKalmanFilter(Y, x0, P0, f, Q1, h, R, 'CKF');
[xf2, Pf2, xp2, Pp2] = nonLinearKalmanFilter(Y, x0, P0, f, Q2, h, R, 'CKF');
[xf3, Pf3, xp3, Pp3] = nonLinearKalmanFilter(Y, x0, P0, f, Q3, h, R, 'CKF');


%% Calculate the relative error

error1 = mean([xf1(1,:)-X(1,2:end); xf1(2,:)-X(2,2:end)]);
error2 = mean([xf2(1,:)-X(1,2:end); xf2(2,:)-X(2,2:end)]);
error3 = mean([xf3(1,:)-X(1,2:end); xf3(2,:)-X(2,2:end)]);

%% Plot true trajectory and sensor positions of all cases
axlim1 = [-50 150 -50 70]; 
axlim2 = [-300 350 -300 600]; 
axlim3 = [-50 150 -50 70]; 

size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.85]); 
subplot(3,2,1)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('Q good tuning'); xlabel('x'); ylabel('y')
hold on;
scatter(x1, y1, 'r.');
plot(xf1(1,:), xf1(2,:), 'm-')
for i = 1:5:K
    xy = sigmaEllipse2D(xf1(1:2, i), Pf1(1:2,1:2, i), level, 100);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim1);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'CKF preduction', 'CKF covariance', 'location', 'best', 'FontSize', 5);

subplot(3,2,3)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('Q with large noise, 1000 times start value'); xlabel('x'); ylabel('y')
hold on;
scatter(x1, y1, 'r.');
plot(xf2(1,:), xf2(2,:), 'm-')
for i = 1:5:K
    xy = sigmaEllipse2D(xf2(1:2, i), Pf2(1:2,1:2, i), level, 100);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim2);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'CKF preduction', 'CKF covariance', 'location', 'best', 'FontSize', 5);

subplot(3,2,5)
plot(X(1,:), X(2,:), 'k--', s1(1), s1(2), 'r*', s2(1), s2(2), 'b*');
title('Q with small noise, 1/1000 times start value'); xlabel('x'); ylabel('y')
hold on;
scatter(x1, y1, 'r.');
plot(xf3(1,:), xf3(3,:), 'm-')
for i = 1:5:K
    xy = sigmaEllipse2D(xf3(1:2, i), Pf3(1:2,1:2, i), level, 100);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(axlim3);
legend('True trajectory', 'Sensor 1', 'Sensor 2', 'Measurements', 'CKF preduction', 'CKF covariance', 'location', 'best', 'FontSize', 5);

% Plot error
subplot(3,2,2)
plot(linspace(0,K/T,length(error1)), error1, 'k-');
title('Relative position error when having good tuning'); xlabel('t');
subplot(3,2,4)
plot(linspace(0,K/T,length(error2)), error2, 'k-');
title('Relative position error when having large noise covariance'); xlabel('t');
subplot(3,2,6)
plot(linspace(0,K/T,length(error3)), error3, 'k-');
title('Relative position error when having small noise covariance'); xlabel('t');

%% Functions %%

%% Motionmodel

function [fx, Fx] = coordinatedTurnMotion(x, T)
%COORDINATEDTURNMOTION calculates the predicted state using a coordinated
%turn motion model, and also calculated the motion model Jacobian
%
%Input:
%   x           [5 x 1] state vector
%   T           [1 x 1] Sampling time
%
%Output:
%   fx          [5 x 1] motion model evaluated at state x
%   Fx          [5 x 5] motion model Jacobian evaluated at state x
%
% NOTE: the motion model assumes that the state vector x consist of the
% following states:
%   px          X-position
%   py          Y-position
%   v           velocity
%   phi         heading
%   omega       turn-rate


% Defining the given variables
p_x = x(1);
p_y = x(2);
v = x(3);
phi = x(4);
omega = x(5);

% Your code for the motion model here
fx = [p_x + T*v*cos(phi);
      p_y + T*v*sin(phi);
      v;
      phi + T*omega;
      omega];

%Check if the Jacobian is requested by the calling function
if nargout > 1
    % Your code for the motion model Jacobian here
    Fx = [1, 0, T*cos(phi), -T*v*sin(phi), 0;
          0, 1, T*sin(phi), T*v*cos(phi), 0;
          0, 0, 1, 0, 0;
          0, 0, 0, 1, T;
          0, 0, 0, 0, 1];
end

end



%% Measurementmodel

function [hx, Hx] = dualBearingMeasurement(x, s1, s2)
%DUOBEARINGMEASUREMENT calculates the bearings from two sensors, located in 
%s1 and s2, to the position given by the state vector x. Also returns the
%Jacobian of the model at x.
%
%Input:
%   x           [n x 1] State vector, the two first element are 2D position
%   s1          [2 x 1] Sensor position (2D) for sensor 1
%   s2          [2 x 1] Sensor position (2D) for sensor 2
%
%Output:
%   hx          [2 x 1] measurement vector
%   Hx          [2 x n] measurement model Jacobian
%
% NOTE: the measurement model assumes that in the state vector x, the first
% two states are X-position and Y-position.

% Define syms to calc the jacobian using the built in command
syms X Y real

% Calc the measurement vector
angle1 = atan2(Y - s1(2), X - s1(1));
angle2 = atan2(Y - s2(2), X - s2(1));
hx = [angle1; angle2];

% Calc the measurement model Jacobian
Hx = jacobian(hx, [X;Y]);

% Convert to double from syms
hx = double(subs(hx,[X;Y],[x(1);x(2)])); 
Hx = double(subs(Hx,[X;Y],[x(1);x(2)])); 

% If size is wrong we add zeros to the end for both sensors
if size(Hx, 2) < size(x, 1)
    Hx = [Hx, zeros(2, size(x, 1) - size(Hx, 2))];
end
end



%% State sequence

function X = genNonLinearStateSequence(x_0, P_0, f, Q, N)
%GENNONLINEARSTATESEQUENCE generates an N+1-long sequence of states using a 
%    Gaussian prior and a nonlinear Gaussian process model
%
%Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   N           [1 x 1] Number of states to generate
%
%Output:
%   X           [n x N+1] State vector sequence
%

% Initialize output array
n = length(x_0);
X = zeros(n, N+1);

% Set initial state
X(:,1) = mvnrnd(x_0, P_0)';

% Generate state sequence
for i = 2:N+1
    [fx, Fx] = f(X(:,i-1));
    X(:,i) = mvnrnd(fx', Q)';
end

end



%% Measurement sequence

function Y = genNonLinearMeasurementSequence(X, h, R)
%GENNONLINEARMEASUREMENTSEQUENCE generates ovservations of the states 
% sequence X using a non-linear measurement model.
%
%Input:
%   X           [n x N+1] State vector sequence
%   h           Measurement model function handle
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state) 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% Predefining the measurement sequence
n = size(X, 1);  % state dimension
T = size(X, 2) - 1;  % number of time steps
m = size(R, 1);  % measurement dimension
Y = zeros(m, T);

for t = 1:T
    x = X(:, t+1);
    Y(:, t) = h(x) + mvnrnd(zeros(m,1), R)';
end

end



%% Sigma points 

function [SP, W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%

n = length(x);

switch type
    
    case 'UKF'
        W0 = 1 - n/3;
        SP = [x, x + sqrt(n/(1-W0))*sqrtm(P), x - sqrt(n/(1-W0))*sqrtm(P)];
        W = [W0, (1-W0)/(2*n) * ones(1, 2*n)]; 
        
    case 'CKF'      
        SP = [x + sqrt(n) * sqrtm(P), x - sqrt(n) * sqrtm(P)];
        W = 1/(2*n) * ones(1, 2*n);   
        
    otherwise
        error('Incorrect type of sigma point')
end
end



%% Prediction step

function [x, P] = nonLinKFprediction(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%

n = length(x);

    switch type
        
        case 'EKF'           
            [fx, Fx] = f(x);
            x = fx;
            P = Fx*P*Fx' + Q;     
            
        case 'UKF'            
            [SP, W] = sigmaPoints(x, P, 'UKF');
            x = 0;
            for i = 1:2*n+1
                x = x + f(SP(:,i)) * W(i);
            end
            P = Q;
            for i = 1:2*n+1
                P = P + (f(SP(:,i)) - x) * (f(SP(:,i)) - x)' * W(i);
            end
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end 
            
        case 'CKF'          
            [SP, W] = sigmaPoints(x, P, 'CKF');
            x = 0;
            for i = 1:2*n
                x = x + f(SP(:,i)) * W(i);
            end
            P = Q;
            for i = 1:2*n
                P = P + (f(SP(:,i)) - x) * (f(SP(:,i)) - x)' * W(i);
            end    
      
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end
end



%% Update step

function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement vector
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state), 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%               Function must include all model parameters for the particular model, 
%               such as sensor position for some models.
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%

n = length(x);

    switch type
        
        case 'EKF'
            [hx,Hx] = h(x);
            S = Hx*P*Hx' + R;
            K = P*Hx'/S;
            x = x + K*(y - hx);
            P = (eye(size(K,1)) - K*Hx) * P; 
            
        case 'UKF'
            [SP, W] = sigmaPoints(x, P, type);
            y_hat = 0;
            for i = 1:2*n+1
                y_hat = y_hat + h(SP(:,i)) * W(i);
            end
            Pxy = 0;
            S = R;
            for i = 1:2*n+1
                Pxy = Pxy + (SP(:,i) - x) * (h(SP(:,i)) - y_hat)' * W(i);
                S = S + (h(SP(:,i)) - y_hat) * (h(SP(:,i)) - y_hat)' * W(i);
            end
            x = x + Pxy/S*(y - y_hat);
            P = P - Pxy/S*Pxy';
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end   
            
        case 'CKF'
            [SP, W] = sigmaPoints(x, P, type);
            y_hat = 0;
            for i = 1:2*n
                y_hat = y_hat + h(SP(:,i)) * W(i);
            end
            Pxy = 0;
            S = R;
            for i = 1:2*n
                Pxy = Pxy + (SP(:,i) - x) * (h(SP(:,i)) - y_hat)' * W(i);
                S = S + (h(SP(:,i)) - y_hat) * (h(SP(:,i)) - y_hat)' * W(i);
            end
            x = x + Pxy/S*(y - y_hat);
            P = P - Pxy/S*Pxy';
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end
end



%% Kalman filter

function [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, type)
%NONLINEARKALMANFILTER Filters measurement sequence Y using a 
% non-linear Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence for times 1,...,N
%   x_0         [n x 1] Prior mean for time 0
%   P_0         [n x n] Prior covariance
%   f                   Motion model function handle
%                       [fx,Fx]=f(x) 
%                       Takes as input x (state) 
%                       Returns fx and Fx, motion model and Jacobian evaluated at x
%   Q           [n x n] Process noise covariance
%   h                   Measurement model function handle
%                       [hx,Hx]=h(x,T) 
%                       Takes as input x (state), 
%                       Returns hx and Hx, measurement model and Jacobian evaluated at x
%   R           [m x m] Measurement noise covariance
%
%Output:
%   xf          [n x N]     Filtered estimates for times 1,...,N
%   Pf          [n x n x N] Filter error convariance
%   xp          [n x N]     Predicted estimates for times 1,...,N
%   Pp          [n x n x N] Filter error convariance
%

N = size(Y, 2); % number of time steps
n = length(x_0); % dimension of state vector

% Initialize output variables
xf = zeros(n, N);
Pf = zeros(n, n, N);
xp = zeros(n, N);
Pp = zeros(n, n, N);

% Set initial predicted state and covariance
xf(:, 1) = x_0;
Pf(:, :, 1) = P_0;

% Filter the measurements using the specified non-linear Kalman filter
for i = 1:N
    % Prediction step
    [xpred, Ppred] = nonLinKFprediction(xf(:, i), Pf(:, :, i), f, Q, type);
    
    % Update step
    [xfup, Pfup] = nonLinKFupdate(xpred, Ppred, Y(:, i), h, R, type);
    
    % Save results
    xp(:,i) = xpred;
    Pp(:,:,i) = Ppred;
    xf(:,i+1) = xfup;
    Pf(:,:,i+1) = Pfup;
end

% Use all but first value
xf = xf(:, 2:end);
Pf = Pf(:, :, 2:end);

end



%% Approximate mean and covariance

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

% Draw N samples from the Gaussian density with parameters given in the input parameters and with the hint function: mvnrnd() 
x_s = mvnrnd(mu_x, Sigma_x, N)';

% Applying the non-linear function to each x in x_s, and concatenate the resulting y vectors into y_s
y_s = f(x_s);

% Calculating the mean of each row and covariance of y_s
mu_y = mean(y_s, 2);
Sigma_y = (y_s - mu_y)*(y_s - mu_y)'/(N-1);

end


%% To plot the mean and covariance
function [xy] = sigmaEllipse2D(mu, Sigma, level, npoints)
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

%% Given function to get pos
function [x, y] = getPosFromMeasurement(y1, y2, s1, s2)
%GETPOSFROMMEASUREMENT computes the intersection point
%(transformed 2D measurement in Cartesian coordinate
%system) given two sensor locations and two bearing
%measurements, one from each sensor.
%INPUT: y1: bearing measurement from sensor 1
% y2: bearing measurement from sensor 2
% s1: location of sensor 1 in 2D Cartesian
% s2: location of sensor 2 in 2D Cartesian
%OUTPUT: x: coordinate of intersection point on x axis
% y: coordinate of intersection point on y axis
%This problem can be formulated as solving a set of two
%linear equations with two unknowns. Specifically, one
%would like to obtain (x,y) by solving
%(y−s1(2))=(x−s1(1))tan(y1) and (y−s2(2))=(x−s2(1))tan(y2).

x = (s2(2)-s1(2)+tan(y1)*s1(1)-tan(y2)*s2(1))/(tan(y1)-tan(y2));
y = s1(2)+tan(y1)*(x-s1(1));

end
