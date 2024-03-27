%% HA4 Daniel Söderqvist

clear all; close all; clc

%% Problem 1

%% a)

%% True track, given from task
% Sampling period
T = 0.1;
% Length of time sequence
K = 600;
% Allocate memory
omega = zeros(1,K+1);
% Turn rate
omega(150:450) = -pi/301/T;
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


%% Definitions

x_0 = [0; 0; 0; 0; 0];
P_0 = diag([10^2, 10^2, 10^2, (5*pi/180)^2, (pi/180)^2]);
s1 = [300; -100];
s2 = [300; -300];
sigmaphi = pi/180;
R = diag([sigmaphi^2, sigmaphi^2]);

% The corresponding process noise covariance matrix is the same for all the
% cases here since the velocity v and angular velocity omega are 0. 
Gamma = [0 0;0 0;1 0;0 0;0 1];
P_v = 1;
P_w = pi/180;
Q = Gamma*[T*P_v^2 0;0 T*P_w^2]*Gamma';


%% Define functions needed as handles
f = @(x) coordinatedTurnMotion(x, T);
h = @(x) dualBearingMeasurement(x, s1, s2);
S = @sigmaPoints;

%% Generate measurements and states
Y = genNonLinearMeasurementSequence(X, h, R);
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x_0, P_0, f, T, Q, s1, s2, h, R, S, 'CKF');


%% Plot results

size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.05, size(2)*0.06, size(1)*0.85, size(2)*0.85]); 

subplot(2,2,1)
ax1 = [-20 520 -450 50];
plot(X(1,:),X(2,:),'k') % True state
grid on; hold on;
plot(xf(1,:), xf(2,:),'b') % Filtered state
plot(xs(1,:), xs(2,:),'m') % Smoothed state
plot(s1(1), s1(2), '*','color',[0.8500 0.3250 0.0980],'linewidth',15) % Sensor position
plot(s2(1), s2(2), '*','color',[0.9290 0.6940 0.1250],'linewidth',15) % Sensor position
for i = 1:length(Y)
     [x1(i), y1(i)] = getPosFromMeasurement(Y(1,i), Y(2,i), s1, s2);
end
scatter(x1,y1, 'r.')
level = 3;
for i = 1:5:K
    xy = sigmaEllipse2D(xf(1:2, i), Pf(1:2,1:2, i), level, 100);
    plot(xy(1, :),xy(2, :),'g--','linewidth',0.5);
end
axis(ax1);
legend('True state', 'Filtered state', 'Smoothed state','Sensor 1 position','Sensor 2 position', 'Measurements', 'Covariances', 'location','best')
title('Without an outlier'); xlabel('x'); ylabel('y')

%% Error 
error_kf = X(:,2:end)-xf;
error_xs = X(:,2:end)-xs;

%% Plot error
subplot(2,2,2)
plot(error_kf(1,:), 'b-', 'LineWidth', 1.5);
hold on; grid on;
plot(error_xs(1,:), 'm-', 'LineWidth', 1.5);
title('Error between True Trajectory and Estimates');
xlabel('Timestep'); ylabel('Error')
legend('Filtered', 'Smoothered');

%% b)

%% Add outlier and then generate states
kk = 300;
outlier = 5;
Y_outlier = Y;
Y_outlier(:, kk) = Y_outlier(:, kk) + outlier;
[xs_outlier, Ps_outlier, xf_outlier, Pf_outlier, xp_outlier, Pp_outlier] = nonLinRTSsmoother(Y_outlier, x_0, P_0, f, T, Q, s1, s2, h, R, S, 'CKF');

%% plot
subplot(2,2,3)
ax2 = [-20 600 -450 100];
plot(X(1,:),X(2,:),'k') % True state
grid on; hold on;
plot(xf_outlier(1,:), xf_outlier(2,:),'b') % Filtered state
plot(xs_outlier(1,:), xs_outlier(2,:),'m') % Smoothed state
plot(s1(1), s1(2), '*','color',[0.8500 0.3250 0.0980],'linewidth',15) % Sensor position
plot(s2(1), s2(2), '*','color',[0.9290 0.6940 0.1250],'linewidth',15) % Sensor position

for i = 1:length(Y_outlier)
     [x1_outlier(i), y1_outlier(i)] = getPosFromMeasurement(Y_outlier(1,i), Y_outlier(2,i), s1, s2);
end
scatter(x1_outlier, y1_outlier, 'r.')
level = 3;
for i = 1:5:K
    xy_outlier = sigmaEllipse2D(xf_outlier(1:2, i), Pf_outlier(1:2,1:2, i), level, 100);
    plot(xy_outlier(1, :),xy_outlier(2, :),'g--','linewidth',0.5);
end
axis(ax2);
legend('True state', 'Filtered state', 'Smoothed state','Sensor 1 position','Sensor 2 position', 'Measurements', 'Covariances', 'location','best')
title('With an outlier'); xlabel('x'); ylabel('y')

%% Error outlier
error_kf_outlier = X(:,2:end)-xf_outlier;
error_xs_outlier = X(:,2:end)-xs_outlier;

%% Plot error
subplot(2,2,4)
plot(error_kf_outlier(1,:), 'b-', 'LineWidth', 1.5);
hold on; grid on;
plot(error_xs_outlier(1,:), 'm-', 'LineWidth', 1.5);

% xlabel('Time Step');
% ylabel('Error');
title('Error between True Trajectory and Estimates');
xlabel('Timestep'); ylabel('Error')
legend('Filtered', 'Smoothered');

%% Problem 2
close all; clear all; clc

%% a)
%% Definition of variables
T = 0.01;
Q = 1.5;
R = 3;
x_0 = 2;
P_0 = 8;
n = 30; % Timesteps
N = 300; % Number of particles
A = 1;
H = 1;

%% Define a linear function
f = @(x) x;

%% Generate linear states and measurements
X = genLinearStateSequence(x_0, P_0, A, Q, n);
Y = genLinearMeasurementSequence(X, H, R);

%% Kalman estimates
[xf, Pf] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

%% Particle filter estimates
% Without resampling
bResample = false;
sigma = 1;
plotFunc_handle0 = @(k, Xk, Xkmin1, Wk, j) plotPostPdf(k, Xk, Wk, xf, Pf, bResample, sigma, gca);
[xfp1, Pfp1, Xp1, Wp1] = pfFilter(x_0, P_0, Y, f, Q, f, R, N, bResample, []);

% With resampling
bResample = true;
plotFunc_handle1 = @(k, Xk, Xkmin1, Wk, j) plotPostPdf(k, Xk, Wk, xf, Pf, bResample, sigma, gca);
[xfp2, Pfp2, Xp2, Wp2] = pfFilter(x_0, P_0, Y, f, Q, f, R, N, bResample, []);

%% Write and calculate MSE to command window
MSE_Kf = mse(X(2:end)-xf);
MSE_Pf0 = mse(X(2:end)-xfp1);
MSE_Pf1 = mse(X(2:end)-xfp2);

fprintf('MSE for the Kalman Filter is:\n\n'); disp(MSE_Kf)
fprintf('MSE for the Particle Filter without resampling is:\n\n'); disp(MSE_Pf0)
fprintf('MSE for the Particle Filter with resampling is:\n\n'); disp(MSE_Pf1)

%% Calculate the error
error_kf = X(2:end)-xf;
error_pf0 = X(2:end)-xfp1;
error_pf1 = X(2:end)-xfp2;

%% Plot 
size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.85]); 

plot(X,'k')
grid on; hold on;
errorbar(0:length(xf(1,:))-1 ,xf(1,:), error_kf, 'm-')
errorbar(0:length(xfp1(1,:))-1 ,xfp1(1,:), error_pf0, 'g-')
errorbar(0:length(xfp2(1,:))-1 ,xfp2(1,:), error_pf1, 'b-')
scatter(1:length(Y), Y, 100, 'r.')
legend('True state', 'Kalman Filter', 'Particle Filter without resample', 'Particle Filter with resample', 'Measurements', 'location','best')
title('Comparisson between PF and KF with right prior')
xlabel('Time-step'); ylabel('State value')

%% Plot densities
size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.85]); 

ax1 = [-10 3 0 0.4];
ax2 = [-8 5 0 0.4];
ax3 = [-12 1 0 0.4];
ax4 = [-10 1 0 0.4];
ax5 = [-7 4 0 0.4];
ax6 = [-10 0 0 0.4];

timesteps = [5 15 25];
subplot(2,3,1)
plotFunc_handle0(timesteps(1), Xp1(:,:, timesteps(1)), Xp1(:, :, timesteps(1)-1), Wp1(:, timesteps(1))', 0);
hold on; grid on;
plot(X(timesteps(1)+1)*ones(100),linspace(0, 0.4, 100),'k--')
legend('Particle filter without resampling', 'Kalman filter','True state', 'Location', 'best')
axis(ax1)

subplot(2,3,2)
plotFunc_handle0(timesteps(2), Xp1(:,:, timesteps(2)), Xp1(:, :, timesteps(2)-1), Wp1(:, timesteps(2))', 0);
hold on; grid on;
plot(X(timesteps(2)+1)*ones(100),linspace(0, 0.4, 100),'k--')
legend('Particle filter without resampling', 'Kalman filter','True state', 'Location', 'best')
axis(ax2)

subplot(2,3,3)
plotFunc_handle0(timesteps(3), Xp1(:,:, timesteps(3)), Xp1(:, :, timesteps(3)-1), Wp1(:, timesteps(3))', 0);
hold on; grid on;
plot(X(timesteps(2)+1)*ones(100),linspace(0, 0.4, 100),'k--')
legend('Particle filter without resampling', 'Kalman filter','True state', 'Location', 'best')
axis(ax3)

subplot(2,3,4)
plotFunc_handle1(timesteps(1), Xp2(:,:, timesteps(1)), Xp2(:, :, timesteps(1)-1), Wp2(:, timesteps(1))', 0);
hold on; grid on;
plot(X(timesteps(1)+1)*ones(100),linspace(0, 0.4, 100),'k--')
legend('Particle filter with resampling', 'Kalman filter','True state', 'Location', 'best')
axis(ax4)

subplot(2,3,5)
plotFunc_handle1(timesteps(2), Xp2(:,:, timesteps(2)), Xp2(:, :, timesteps(2)-1), Wp2(:, timesteps(2))', 0);
hold on; grid on;
plot(X(timesteps(2)+1)*ones(100),linspace(0, 0.4, 100),'k--')
legend('Particle filter with resampling', 'Kalman filter','True state', 'Location', 'best')
axis(ax5)

subplot(2,3,6)
plotFunc_handle1(timesteps(3), Xp2(:,:, timesteps(3)), Xp2(:, :, timesteps(3)-1), Wp2(:, timesteps(3))', 0);
hold on; grid on;
plot(X(timesteps(2)+1)*ones(100),linspace(0, 0.4, 100),'k--')
legend('Particle filter with resampling', 'Kalman filter','True state', 'Location', 'best')
axis(ax6)


%% b)
%% Wrong prior
x00 = -20;
P00 = 2;

%% Kalman estimates
[xff, Pff] = kalmanFilter(Y, x00, P00, A, Q, H, R);

%% Particle filter estimates
bResample = false;
[xffp1, Pffp1, Xpp1, Wpp1] = pfFilter(x00, P00, Y, f, Q, f, R, N, bResample, []);
bResample = true;
[xffp2, Pffp2, Xpp2, Wpp2] = pfFilter(x00, P00, Y, f, Q, f, R, N, bResample, []);

%% Calculate errors
errorr_kf = X(2:end)-xff;
errorr_pf0 = X(2:end)-xffp1;
errorr_pf1 = X(2:end)-xffp2;

%% Plot
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.8]); 
subplot(2,1,1)
plot(X,'k')
grid on; hold on;
scatter(1:length(Y), Y, 100, 'r.')
plot(xff, 'm--')
plot(xffp1(1,:), 'g--')
plot(xffp2(1,:), 'b--')
legend('True state', 'Measurements', 'Kalmanfilter', 'Particlefilter without resample', 'Particlefilter with resample','location','best')
title('Comparisson between PF and KF with wrong prior')
xlabel('Time-step'); ylabel('State value')

subplot(2,1,2)
plot(X,'k')
grid on; hold on;
errorbar(0:length(xff(1,:))-1 ,xff(1,:), errorr_kf, 'm-')
errorbar(0:length(xffp1(1,:))-1 ,xffp1(1,:), errorr_pf0, 'g-')
errorbar(0:length(xffp2(1,:))-1 ,xffp2(1,:), errorr_pf1, 'b-')
scatter(1:length(Y), Y, 100, 'r.')
legend('True state', 'Kalman Filter', 'Particle Filter without resample', 'Particle Filter with resample', 'Measurements', 'location','best')
title('Comparisson between PF and KF with right prior')
xlabel('Time-step'); ylabel('State value')


%% c) and d)
%% Plot particle trajectories for with and without resampling
N = 100; % Number of particles
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.85]); 
plot_trajs = @(k, Xk, Xkmin1, j) plotPartTrajs(k, Xk, Xkmin1, j);
subplot(2,1,1)
plot(0:length(X)-1, X(1:end), 'r-','linewidth',2)
grid on; hold on;
bResample = false;
pfFilter(x_0, P_0, Y, f, Q, f, R, N, bResample, plot_trajs);
plot(0:length(X)-1, X(1:end), 'r-','linewidth',2)
xlabel('Time-step'); ylabel('State value')
legend('True trajectory', 'Particles')

subplot(2,1,2)
plot(0:length(X)-1, X(1:end), 'r-','linewidth',2)
hold on; grid on;
bResample = true;
pfFilter(x_0, P_0, Y, f, Q, f, R, N, bResample, plot_trajs);
plot(0:length(X)-1, X(1:end), 'r-','linewidth',2)
xlabel('Time-step'); ylabel('State value')
legend('True trajectory', 'Particles')


%% Problem 3
clear all; close all; clc

%% a)
% Load data
load("Xk.mat");

% Calculate the position sequence
position_sequence = Xk;

% Calculate the velocity measurement sequence
velocity_sequence_x = diff(X);
velocity_sequence_y = diff(Y);
velocity_sequence = [velocity_sequence_x'; velocity_sequence_y'];

% Add Gaussian noise to the velocity measurements
r = 0.001;
sigma_r = diag([r r]); % Adjust the value of sigma_r as needed
noise = mvnrnd([0 0], sigma_r, length(velocity_sequence_x));
velocity_measurement_sequence = velocity_sequence + noise';

%% Plot
size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.4, size(2)*0.85]); 

subplot(2,1,1)
hold on
plot([1+i 1+9*i 5+9*i],'b')
plot([7+9*i 11+9*i 11+i 7+i],'b');plot([5+i 1+i],'b')
plot([2+5.2*i 2+8.3*i 4+8.3*i 4+5.2*i 2+5.2*i],'b')%House 1
plot([2+3.7*i 2+4.4*i 4+4.4*i 4+3.7*i 2+3.7*i],'b')%House 2
plot([2+2*i 2+3.2*i 4+3.2*i 4+2*i 2+2*i],'b')%House 3
plot([5+i 5+2.2*i 7+2.2*i 7+i],'b')%House 4
plot([5+2.8*i 5+5.5*i 7+5.5*i 7+2.8*i 5+2.8*i],'b')%House 5
plot([5+6.2*i 5+9*i],'b');plot([7+9*i 7+6.2*i 5+6.2*i],'b')%House 6
plot([8+4.6*i 8+8.4*i 10+8.4*i 10+4.6*i 8+4.6*i],'b')%House 7
plot([8+2.4*i 8+4*i 10+4*i 10+2.4*i 8+2.4*i],'b')%House 8
plot([8+1.7*i 8+1.8*i 10+1.8*i 10+1.7*i 8+1.7*i],'b')%House 9


title('A map of the village','FontSize',20); xlabel('x'); ylabel('y')
plot(Xk(1,:),Xk(2,:), 'r-*')
subplot(2,1,2)
plot(velocity_measurement_sequence(1,:),'b-')
hold on; grid on;
plot(velocity_measurement_sequence(2,:),'m-')
title('Velocity profile'); xlabel('Timestep'); ylabel('Velocity')
legend('Velocity in x coordinate', 'Velocity in Y coordinate')

%% Needed functions

%% Motion model
function [fx, Fx] = coordinatedTurnMotion(x, T)

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


%% Measurement model
function [hx, Hx] = dualBearingMeasurement(x, s1, s2)

% Define syms to calc the jacobiian using the built in command
syms X Y real

% Calc the measurement vector
angle1 = atan2(Y-s1(2), X-s1(1));
angle2 = atan2(Y-s2(2), X-s2(1));
hx = [angle1; angle2];

% Calc the measurement model Jacobian
Hx=jacobian(hx,[X;Y]);

% Convert to double from syms
hx = double(subs(hx,[X;Y],[x(1);x(2)])); 
Hx = double(subs(Hx,[X;Y],[x(1);x(2)])); 

% If size is wrong we add zeros to the end for both sensors
if size(Hx, 2) < size(x, 1)
    Hx = [Hx, zeros(2, size(x, 1) - size(Hx, 2))];
end
end


%% Smoother
function [xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x_0, P_0, f, T, Q, s1, s2, h, R, sigmaPoints, type)

% Kalman from HA3 we estimate forward
N = size(Y, 2); % number of time steps
n = length(x_0); % dimension of state vector

% Predefine output variables
xf = zeros(n, N);
Pf = zeros(n, n, N);
xp = zeros(n, N);
Pp = zeros(n, n, N);

% Save prior
xf(:, 1) = x_0;
Pf(:, :, 1) = P_0;

% Filter the measurements using the specified non-linear Kalman filter
for i = 1:N
    % Prediction step
    [xpred, Ppred] = nonLinKFprediction(xf(:, i), Pf(:, :, i), f, T, Q, sigmaPoints, type);
    
    % Update step
    [xfup, Pfup] = nonLinKFupdate(xpred, Ppred, Y(:, i), s1, s2, h, R, sigmaPoints, type);
    
    % Save results
    xp(:, i) = xpred;
    Pp(:, :, i) = Ppred;
    xf(:, i+1) = xfup;
    Pf(:, :, i+1) = Pfup;
end

% Use all but first value
xf = xf(:, 2:end);
Pf = Pf(:, :, 2:end);

% RTS smoother we estimate backwards to the beginning already calculated all the estimated steps
xs(:, N) = xf(:, N); 
Ps(:, :, N) = Pf(:, :, N);
for i = (N-1):-1:1
    [xs(:,i) Ps(:, :, i)] = nonLinRTSSupdate(xs(:, i+1), Ps(:, :, i+1), xf(:, i), Pf(:, :, i), xp(:, i+1), Pp(:, :, i+1), f, T, sigmaPoints, type);    
end

end


function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, Ps_kplus1, xf_k, Pf_k, xp_kplus1, Pp_kplus1, ...
                                     f, T, sigmaPoints, type)
n = length(xf_k);

    switch type
        
        case 'EKF'           
            [fx, Fx] = f(xf_k, T);
            P = Pf_k*Fx';     
            
        case 'UKF'            
            [SP, W] = sigmaPoints(xf_k, Pf_k, 'UKF');
            P = 0;
            for i = 1:2*n+1
                P = P + (SP(:,i) - xf_k) * (f(SP(:,i)) - xp_kplus1)' * W(i);
            end
            
        case 'CKF'          
            [SP, W] = sigmaPoints(xf_k, Pf_k, 'CKF');
            P = 0;
            for i = 1:2*n
                P = P + (SP(:,i) - xf_k) * (f(SP(:,i)) - xp_kplus1)' * W(i);
            end    
      
        
    end
    
    Gk =  P * Pp_kplus1^(-1);
    xs = xf_k + Gk * ( xs_kplus1 - xp_kplus1 );
    Ps = Pf_k - Gk * ( Pp_kplus1 - Ps_kplus1 ) * Gk';   
end


function [x, P] = nonLinKFprediction(x, P, f, T, Q, sigmaPoints, type)

    switch type
        case 'EKF'

            % Evaluate motion model
            [fx, Fx] = f(x,T);
            % State prediction
            x = fx;
            % Covariance prediciton
            P = Fx*P*Fx' + Q;
            % Make sure P is symmetric
            P = 0.5*(P + P');

        case 'UKF'

            % Predict
            [x, P] = predictMeanAndCovWithSigmaPoints(x, P, f, T, Q, sigmaPoints, type);

            if min(eig(P))<=0
                [v,e] = eig(P);
                emin = 1e-3;
                e = diag(max(diag(e),emin));
                P = v*e*v';
            end

        case 'CKF'

            % Predict
            [x, P] = predictMeanAndCovWithSigmaPoints(x, P, f, T, Q, sigmaPoints, type);

        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end
end


function [x, P] = nonLinKFupdate(x, P, y, s1, s2, h, R, sigmaPoints, type)

switch type
    case 'EKF'
        
        % Evaluate measurement model
        [hx, Hx] = h(x, s1, s2);
        
        % Innovation covariance
        S = Hx*P*Hx' + R;
        % Kalman gain
        K = (P*Hx')/S;
        
        % State update
        x = x + K*(y - hx);
        % Covariance update
        P = P - K*S*K';
        
        % Make sure P is symmetric
        P = 0.5*(P + P');
        
    case 'UKF'

        % Update mean and covariance
        [x, P] = updateMeanAndCovWithSigmaPoints(x, P, y, s1, s2, h, R, sigmaPoints, type);
        
        if min(eig(P))<=0
            [v,e] = eig(P);
            emin = 1e-3;
            e = diag(max(diag(e),emin));
            P = v*e*v';
        end
        
    case 'CKF'

        % Update mean and covariance
        [x, P] = updateMeanAndCovWithSigmaPoints(x, P, y, s1, s2, h, R, sigmaPoints, type);
        
    otherwise
        error('Incorrect type of non-linear Kalman filter')
end

end


function [x, P] = predictMeanAndCovWithSigmaPoints(x, P, f, T, Q, sigmaPoints, type)

    % Compute sigma points
    [SP,W] = sigmaPoints(x, P, type);

    % Dimension of state and number of sigma points
    [n, N] = size(SP);

    % Allocate memory
    fSP = zeros(n,N);

    % Predict sigma points
    for i = 1:N
        [fSP(:,i), ~] = f(SP(:, i));
    end

    % Compute the predicted mean
    x = sum(fSP.*repmat(W,[n, 1]),2);

    % Compute predicted covariance
    P = Q;
    for i = 1:N
        P = P + W(i)*(fSP(:,i)-x)*(fSP(:,i)-x)';
    end

    % Make sure P is symmetric
    P = 0.5*(P + P');

end


function [x, P] = updateMeanAndCovWithSigmaPoints(x, P, y, s1, s2, h, R, sigmaPoints, type)

    % Compute sigma points
    [SP,W] = sigmaPoints(x, P, type);

    % Dimension of measurement
    m = size(R,1);

    % Dimension of state and number of sigma points
    [n, N] = size(SP);

    % Predicted measurement
    yhat = zeros(m,1);
    hSP = zeros(m,N);
    for i = 1:N
        [hSP(:,i),~] = h(SP(:,i));
        yhat = yhat + W(i)*hSP(:,i);
    end

    % Cross covariance and innovation covariance
    Pxy = zeros(n,m);
    S = R;
    for i=1:N
        Pxy = Pxy + W(i)*(SP(:,i)-x)*(hSP(:,i)-yhat)';
        S = S + W(i)*(hSP(:,i)-yhat)*(hSP(:,i)-yhat)';
    end

    % Ensure symmetry
    S = 0.5*(S+S');

    % Updated mean
    x = x+Pxy*(S\(y-yhat));
    P = P - Pxy*(S\(Pxy'));

    % Ensure symmetry
    P = 0.5*(P+P');
end


%% Sigmapoints
function [SP, W] = sigmaPoints(x, P, type)

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

%% Generate measurements
function Y = genNonLinearMeasurementSequence(X, h, R)

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

%% Covariance ellipse
function [ xy ] = sigmaEllipse2D( mu, Sigma, level, npoints )

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

x = (s2(2)-s1(2)+tan(y1)*s1(1)-tan(y2)*s2(1))/(tan(y1)-tan(y2));
y = s1(2)+tan(y1)*(x-s1(1));

end

%% Generation of lienar measurements
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

%% Generation of lienar state
function X = genLinearStateSequence(x_0, P_0, A, Q, N)

d = length(x_0); % dimension of the state
X = zeros(d, N+1); % matrix to store the state vectors

% Generate the initial state vector
X(:,1) = mvnrnd(x_0, P_0)'; % transpose to obtain a row vector

% Generate the subsequent state vectors
for k = 2:N+1
    X(:,k) = A*X(:,k-1) + mvnrnd(zeros(d,1), Q)'; % transpose to obtain a row vector
end
end

%% Kalman filter
function [X, P] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)

%% Parameters
N = size(Y,2);
n = length(x_0);
m = size(Y,1);

%% Data allocation
x = zeros(n,N);
P = zeros(n,n,N);

% Initialize state estimate
X(:,1) = x_0;
P(:,:,1) = P_0 ;

% Kalman filter loop
for i=1:N
    % Prediction step
    [x_pred, P_pred] = linearPrediction(X(:,i), P(:,:,i), A, Q);
    
    % Update step
    [x_update, P_update] = linearUpdate(x_pred, P_pred, Y(:,i), H, R);
    
    % Save results
    X(:,i+1) = x_update;
    P(:,:,i+1) = P_update;
    
end
X = X(:,2:end);
P = P(:,:,2:end);
end

function [x, P] = linearPrediction(x, P, A, Q)
x = A * x;
P = A * P * A' + Q;
end

function [x, P] = linearUpdate(x, P, y, H, R)
K = P * H' / (H * P * H' + R);
x = x + K * (y - H * x);
P = (eye(size(P)) - K * H) * P;
end

function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, N, bResample, plotFunc)

% Define sizes and predefining output
n = length(x_0);
K = size(Y, 2);
xfp = zeros(n, K);
Pfp = zeros(n, n, K);
Xp = zeros(n, N, K);
Xp_plot = zeros(n, N, K); Xp_plot(:, :, 1) = x_0;
Wp = zeros(N, K);

% Define start
X = mvnrnd(x_0, P_0, N)';
W = 1/N * ones(1, N);
   
for i = 1:K

% Calculate a filter step
[X, W] = pfFilterStep(X, W, Y(:, i), proc_f, proc_Q, meas_h, meas_R);

j = 1:N;

% If true resample
if bResample
    [X, W, j] = resampl(X, W);        
end

% Save results
Xp(:, :, i) = X;
Xp_plot(:, :, i+1) = X;
Wp(:, i) = W;
xfp(:, i) = sum(X.*W, 2);
Pfp(:, :, i) = W.* (X - xfp(:, i)) * (X - xfp(:, i))';

% Plot if wanted
if ~isempty(plotFunc)
    plotFunc(i, Xp_plot(:,:,i+1), Xp_plot(:,:,i), j);
end

end
end

function [Xr, Wr, j] = resampl(X, W)
    n = length(W);
    j = discretize(rand(1, n), cumsum([0 W]));
    Xr = X(:, j);
    Wr = 1/n * ones(1, n);
end

function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)
    X_k = mvnrnd(proc_f(X_kmin1)', proc_Q)';  
    Y_k = mvnpdf(yk', meas_h(X_k)', meas_R)';
    W_k = W_kmin1.*Y_k;
    W_k = W_k/sum(W_k);
end

%% Plot function
function plotPostPdf(k, Xk, Wk, xf, Pf, bResample, sigma, ax)

    N = size(Xk,2);

    % Let us first determine the x-interval of interest:
    xmin =    min(Xk(1,:)); %ax(1);
    xmax =    max(Xk(1,:)); %ax(2); 
    X    =    linspace(xmin-(xmax-xmin)/3, xmax+(xmax-xmin)/3, 800);

    % We can now construct a continuous approximation to the posterior
    % density by placing a Gaussian kernel around each particle
    pApprox = zeros(size(X));   % A vector that will contain the pdf values

    if bResample
        sigma=(xmax-xmin)/sqrt(N);
    end
    
    for i = 1 : N
        pApprox = pApprox + Wk(1,i)*normpdf(Xk(1,i), X, sigma);
    end

    % We are now ready to plot the densities
    
    % figure;
    set(gcf, 'Name', ['p_',num2str(k), '_', 'SIR']);
    % clf
    
    plot(X, pApprox, 'LineWidth', 2)   % This is the PF approximation
    hold on
    plot(X, normpdf(xf(1,k), X, sqrt(Pf(1,1,k))), 'r-.', 'LineWidth', 2) % KF posterior density
%     legend('Particle filter approximation', 'Kalman filter', 'Location', 'southwest')
    title(['p(x_k |  y_{1:k}), k=', num2str(k)])
    hold off;
    
end

%% Plot trajectories
function plotPartTrajs(k, Xk, Xkmin1, j)

    if (size(Xk,2) <= 100) % At most 50 particles may be plotted
        for i = 1:size(Xk,2) % loop through all particles
            plot([k-1 k], [Xkmin1(1,j(i)) Xk(1,i)],'color',[0.4 0.4 0.4]);
            hold on 
        end
        title(['Particle trajectories up to time k=', num2str(k)]);
        pause(0.05);
    else
        disp('Too many particles to plot!'); 
    end
end
