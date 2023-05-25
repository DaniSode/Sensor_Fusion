clear all; close all; clc

run('startup.m')

%% Switch to whatever filter wanted
[xhat, meas] = filterTemplate;
% [xhat, meas] = without_update;
% [xhat, meas] = with_update;

%% Use app on phone

%% Save data
% Save measurements
gyr = meas.gyr;
acc = meas.acc;
mag = meas.mag;
orient = meas.orient;
t_meas = meas.t;

% Convert to euler angles
real_angles = q2euler(orient);

% Save states and covariance
x = xhat.x;
P = xhat.P;
t_xhat = xhat.t;

% Convert to euler angles
est_angles = q2euler(x);

%% Plot 
size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.5, size(2)*0.85]); 

subplot(3,1,1)
plot(t_xhat, est_angles(1,:))
hold on; grid on
plot(t_meas, real_angles(1,:))
title('X axis real vs estimation'); xlabel('Time-step'); ylabel('x value')
legend('Estimated', 'Real')
subplot(3,1,2)
plot(t_xhat, est_angles(2,:))
hold on; grid on
plot(t_meas, real_angles(2,:))
title('Y axis real vs estimation'); xlabel('Time-step'); ylabel('y value')
legend('Estimated', 'Real')
subplot(3,1,3)
plot(t_xhat, est_angles(3,:))
hold on; grid on
plot(t_meas, real_angles(3,:))
title('Z axis real vs estimation'); xlabel('Time-step'); ylabel('z value')
legend('Estimated', 'Real')
