%% Clear all
clear all; close all ; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Main code for project             %%
%%     Daniel Söderqvist and Swadesh Gandhi     %%
%% SSY345 Sensor Fusion and Nonlinear Filtering %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Read data
run('startup.m')

[xhat, meas] = filterTemplate;
acc_data = meas.acc;
gyro_data = meas.gyr;
mag_data = meas.mag;
t = meas.t;

%% Calculate mean and covairance and print them to command window moving
mean_acc = mean(acc_data(:, ~any(isnan(acc_data), 1)), 2);
mean_gyro = mean(gyro_data(:, ~any(isnan(gyro_data), 1)), 2);
mean_mag = mean(mag_data(:, ~any(isnan(mag_data), 1)), 2);

% cov_acc = [cov(acc_data(1, ~any(isnan(acc_data), 1))), cov(acc_data(2, ~any(isnan(acc_data), 1))), cov(acc_data(3, ~any(isnan(acc_data), 1)))];
% cov_gyro = [cov(gyro_data(1, ~any(isnan(gyro_data), 1))), cov(gyro_data(2, ~any(isnan(gyro_data), 1))), cov(gyro_data(3, ~any(isnan(gyro_data), 1)))];
% cov_mag = [cov(mag_data(1, ~any(isnan(mag_data), 1))), cov(mag_data(2, ~any(isnan(mag_data), 1))), cov(mag_data(3, ~any(isnan(mag_data), 1)))];

cov_acc = cov(acc_data(:, ~any(isnan(acc_data), 1))');
cov_gyro = cov(gyro_data(:, ~any(isnan(gyro_data), 1))');
cov_mag = cov(mag_data(:, ~any(isnan(mag_data), 1))');

fprintf('Mean for the accelerometer data when stable in x, y and z respective =\n\n'); disp(mean_acc)
fprintf('Mean for the gyroscope data when stable in x, y and z respective =\n\n'); disp(mean_gyro)
fprintf('Mean for the magnetometer data when stable in x, y and z respective =\n\n'); disp(mean_mag)

fprintf('\n\n\n\nCovariance for the accelerometer data when stable in x, y and z respective =\n\n'); disp(cov_acc)
fprintf('Covariance for the gyroscope data when stable in x, y and z respective =\n\n'); disp(cov_gyro)
fprintf('Covariance for the magnetometer data when stable in x, y and z respective =\n\n'); disp(cov_mag)


%% Plot results
n = 100;
bins = 30;
size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.85]); 

subplot(3,5,4:5)
plot(t, acc_data)
legend('x','y','z')
title('Accelerometer'); xlabel('Time [s]'); ylabel('Sensor reading value')
axis('tight')
subplot(3,5,9:10)
plot(t, gyro_data)
legend('x','y','z')
title('Gyroscope'); xlabel('Time [s]'); ylabel('Sensor reading value')
axis('tight')
subplot(3,5,14:15)
plot(t, mag_data)
legend('x','y','z')
title('Magnetometer'); xlabel('Time [s]'); ylabel('Sensor reading value')
axis('tight')
sgtitle('Sensor data', 'fontsize',20)

subplot(3,5,1)
histogram(acc_data(1, :), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(acc_data(1, :)), max(acc_data(1, :)), n), mvnpdf(linspace(min(acc_data(1, :)), max(acc_data(1, :)), n)', mean_acc(1), cov_acc(1)), 'r', 'linewidth', 1.5)
title('Accelerometer data x')
subplot(3,5,2)
histogram(acc_data(2, :), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(acc_data(2, :)), max(acc_data(2, :)), n), mvnpdf(linspace(min(acc_data(2, :)), max(acc_data(2, :)), n)', mean_acc(2), cov_acc(2)), 'r', 'linewidth', 1.5)
title('Accelerometer data y')
subplot(3,5,3)
histogram(acc_data(3, :), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(acc_data(3, :)), max(acc_data(3, :)), n), mvnpdf(linspace(min(acc_data(3, :)), max(acc_data(3, :)), n)', mean_acc(3), cov_acc(3)), 'r', 'linewidth', 1.5)
title('Accelerometer data z')

subplot(3,5,6)
histogram(gyro_data(1, :), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(gyro_data(1, :)), max(gyro_data(1, :)), n), mvnpdf(linspace(min(gyro_data(1, :)), max(gyro_data(1, :)), n)', mean_gyro(1), cov_gyro(1)), 'r', 'linewidth', 1.5)
title('Gyroscope data x')
subplot(3,5,7)
histogram(gyro_data(2, :), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(gyro_data(2, :)), max(gyro_data(2, :)), n), mvnpdf(linspace(min(gyro_data(2, :)), max(gyro_data(2, :)), n)', mean_gyro(2), cov_gyro(2)), 'r', 'linewidth', 1.5)
title('Gyroscope data y')
subplot(3,5,8)
histogram(gyro_data(3, :), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(gyro_data(3, :)), max(gyro_data(3, :)), n), mvnpdf(linspace(min(gyro_data(3, :)), max(gyro_data(3, :)), n)', mean_gyro(3), cov_gyro(3)), 'r', 'linewidth', 1.5)
title('Gyroscope data z')

subplot(3,5,11)
histogram(mag_data(1, :), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(mag_data(1, :)), max(mag_data(1, :)), n), mvnpdf(linspace(min(mag_data(1, :)), max(mag_data(1, :)), n)', mean_mag(1), cov_mag(1)), 'r', 'linewidth', 1.5)
title('Magnetometer data x')
subplot(3,5,12)
histogram(mag_data(2, :), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(mag_data(2, :)), max(mag_data(2, :)), n), mvnpdf(linspace(min(mag_data(2, :)), max(mag_data(2, :)), n)', mean_mag(2), cov_mag(2)), 'r', 'linewidth', 1.5)
title('Magnetometer data y')
subplot(3,5,13)
histogram(mag_data(3, :), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(mag_data(3, :)), max(mag_data(3, :)), n), mvnpdf(linspace(min(mag_data(3, :)), max(mag_data(3, :)), n)', mean_mag(3), cov_mag(3)), 'r', 'linewidth', 1.5)
title('Magnetometer data z')
