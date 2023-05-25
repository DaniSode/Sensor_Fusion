clear all; close all; clc

run('startup.m')

%% Switch to whatever filter wanted
[xhat, meas] = filterTemplate;
% [xhat, meas] = without_update;
% [xhat, meas] = with_update;

%% Use app on phone

%% Save data
gyr = meas.gyr;
acc = meas.acc;
mag = meas.mag;
t_meas = meas.t;

x = xhat.x;
P = xhat.P;
t_xhat = xhat.t;

angles = q2euler(x);

figure

subplot(3,1,1)
plot(t_xhat, angles(1,:))
hold on
plot(t_meas, gyr(1,:))
subplot(3,1,2)
plot(t_xhat, angles(2,:))
hold on
plot(t_meas, gyr(2,:))
subplot(3,1,3)
plot(t_xhat, angles(3,:))
hold on
plot(t_meas, gyr(3,:))

