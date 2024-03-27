%% Clear all
clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Main 2 code for project            %%
%%     Daniel Söderqvist and Swadesh Gandhi     %%
%% SSY345 Sensor Fusion and Nonlinear Filtering %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run startup file
run('startup.m')

%% Switch to whatever filter wanted
% [xhat, meas] = filterTemplate;
% [xhat, meas] = gyro;
% [xhat, meas] = acc;
% [xhat, meas] = mag;
% [xhat, meas] = gyro_acc;
% [xhat, meas] = gyro_mag;
% [xhat, meas] = acc_mag;
% [xhat, meas] = gyro_acc_mag;


%% Use app on phone

%% Save data
% Save measurements
gyr = meas.gyr;
acc = meas.acc;
mag = meas.mag;
orient = meas.orient;
t_meas = meas.t;

% Convert to euler angles
real_angles = q2euler(orient)*(180/pi);

% Save states and covariance
x = xhat.x;
P = xhat.P;
t_xhat = xhat.t;

% Convert to euler angles
est_angles = q2euler(x)*(180/pi);

%% Plot 
size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.5, size(2)*0.85]); 

subplot(3,1,1)
plot(t_xhat, est_angles(1,:))
hold on; grid on
plot(t_meas, real_angles(1,:))
title('Roll, real vs estimation'); xlabel('Time [s]'); ylabel('Angle [degrees]')
legend('Estimated', 'Google')
subplot(3,1,2)
plot(t_xhat, est_angles(2,:))
hold on; grid on
plot(t_meas, real_angles(2,:))
title('Pitch, real vs estimation'); xlabel('Time [s]'); ylabel('Angle [degrees]')
legend('Estimated', 'Google')
subplot(3,1,3)
plot(t_xhat, est_angles(3,:))
hold on; grid on
plot(t_meas, real_angles(3,:))
title('Yaw, real vs estimation'); xlabel('Time [s]'); ylabel('Angle [degrees]')
legend('Estimated', 'Google')
