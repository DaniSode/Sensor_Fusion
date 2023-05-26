%% Clear all
clear all; close all ; clc

%% Optional %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Main code for project             %%
%%     Daniel Söderqvist and Swadesh Gandhi     %%
%% SSY345 Sensor Fusion and Nonlinear Filtering %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Read data
[numbers_moving, strings_moving, ~] = xlsread('data_moving.xlsx');
[numbers_still, strings_still, ~] = xlsread('data_still.xlsx');


%% Handle data moving 
acc_data_moving = [];
gyro_data_moving = [];
mag_data_moving = [];

for i = 1:length(strings_moving)
     if "ACC" == strings_moving(i) 
        acc_data_moving = [acc_data_moving; numbers_moving(i, 3:5)];
     elseif "GYR" == strings_moving(i) 
        gyro_data_moving = [gyro_data_moving; numbers_moving(i, 3:5)];
     elseif "MAG" == strings_moving(i) 
        mag_data_moving = [mag_data_moving; numbers_moving(i, 3:5)];
     else
         continue
     end
end

min_length_moving = min([size(acc_data_moving, 1), size(gyro_data_moving, 1), size(mag_data_moving, 1)]);
acc_data_moving = acc_data_moving(1:min_length_moving, :);
gyro_data_moving = gyro_data_moving(1:min_length_moving, :);
mag_data_moving = mag_data_moving(1:min_length_moving, :);

timestep = 0.01;
t_moving = (1:min_length_moving).*timestep;


%% Handle data still
acc_data_still = [];
gyro_data_still = [];
mag_data_still = [];
for i = 1:length(strings_still)
     if "ACC" == strings_still(i) 
        acc_data_still = [acc_data_still; numbers_still(i,3:5)];
     elseif "GYR" == strings_still(i) 
        gyro_data_still = [gyro_data_still; numbers_still(i,3:5)];
     elseif "MAG" == strings_still(i) 
        mag_data_still = [mag_data_still; numbers_still(i,3:5)];
     else
         continue
     end
end

min_length_still = min([size(acc_data_still, 1), size(gyro_data_still, 1), size(mag_data_still, 1)]);
acc_data_still = acc_data_still(1:min_length_still, :);
gyro_data_still = gyro_data_still(1:min_length_still, :);
mag_data_still = mag_data_still(1:min_length_still, :);

t_still = (1:min_length_still).*timestep;


%% Calculate mean and covairance and print them to command window moving
mean_acc_moving = mean(acc_data_moving);
mean_gyro_moving = mean(gyro_data_moving);
mean_mag_moving = mean(mag_data_moving);

cov_acc_moving = [cov(acc_data_moving(:,1)), cov(acc_data_moving(:,2)), cov(acc_data_moving(:,3))];
cov_gyro_moving = [cov(gyro_data_moving(:,1)), cov(gyro_data_moving(:,2)), cov(gyro_data_moving(:,3))];
cov_mag_moving = [cov(mag_data_moving(:,1)), cov(mag_data_moving(:,2)), cov(mag_data_moving(:,3))];

fprintf('Mean for the accelerometer data when moving in x, y and z respective =\n\n'); disp(mean_acc_moving)
fprintf('Mean for the gyroscope data when moving in x, y and z respective =\n\n'); disp(mean_gyro_moving)
fprintf('Mean for the magnetometer data when moving in x, y and z respective =\n\n'); disp(mean_mag_moving)

fprintf('\n\n\n\nCovariance for the accelerometer data when moving in x, y and z respective =\n\n'); disp(cov_acc_moving)
fprintf('Covariance for the gyroscope data when moving in x, y and z respective =\n\n'); disp(cov_gyro_moving)
fprintf('Covariance for the magnetometer data when moving in x, y and z respective =\n\n'); disp(cov_mag_moving)


%% Calculate mean and covairance and print them to command window still
mean_acc_still = mean(acc_data_still);
mean_gyro_still = mean(gyro_data_still);
mean_mag_still = mean(mag_data_still);

% cov_acc_still = [cov(acc_data_still(:,1)), cov(acc_data_still(:,2)), cov(acc_data_still(:,3))];
% cov_gyro_still = [cov(gyro_data_still(:,1)), cov(gyro_data_still(:,2)), cov(gyro_data_still(:,3))];
% cov_mag_still = [cov(mag_data_still(:,1)), cov(mag_data_still(:,2)), cov(mag_data_still(:,3))];

cov_acc_still = cov(acc_data_still);
cov_gyro_still = cov(gyro_data_still);
cov_mag_still = cov(mag_data_still);

fprintf('\n\nMean for the accelerometer data when still in x, y and z respective =\n\n'); disp(mean_acc_still)
fprintf('Mean for the gyroscope data when still in x, y and z respective =\n\n'); disp(mean_gyro_still)
fprintf('Mean for the magnetometer data when still in x, y and z respective =\n\n'); disp(mean_mag_still)

fprintf('\n\nCovariance for the accelerometer data when still in x, y and z respective =\n\n'); disp(cov_acc_still)
fprintf('Covariance for the gyroscope data when still in x, y and z respective =\n\n'); disp(cov_gyro_still)
fprintf('Covariance for the magnetometer data when still in x, y and z respective =\n\n'); disp(cov_mag_still)


%% Plot results moving
n = 100;
bins = 30;
size = get(0,'screensize'); size = size(1,end-1:end);
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.85]); 

subplot(3,5,4:5)
plot(t_moving, acc_data_moving)
legend('x','y','z')
title('Accelerometer'); xlabel('Time [s]'); ylabel('Sensor reading value')
axis('tight')
subplot(3,5,9:10)
plot(t_moving, gyro_data_moving)
legend('x','y','z')
title('Gyroscope'); xlabel('Time [s]'); ylabel('Sensor reading value')
axis('tight')
subplot(3,5,14:15)
plot(t_moving, mag_data_moving)
legend('x','y','z')
title('Magnetometer'); xlabel('Time [s]'); ylabel('Sensor reading value')
axis('tight')
sgtitle('When the phone is moving', 'fontsize',20)

subplot(3,5,1)
histogram(acc_data_moving(:, 1), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(acc_data_moving(:, 1)), max(acc_data_moving(:, 1)), n), mvnpdf(linspace(min(acc_data_moving(:, 1)), max(acc_data_moving(:, 1)), n)', mean_acc_moving(1), cov_acc_moving(1)), 'r', 'linewidth', 1.5)
title('Accelerometer data x')
subplot(3,5,2)
histogram(acc_data_moving(:, 2), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(acc_data_moving(:, 2)), max(acc_data_moving(:, 2)), n), mvnpdf(linspace(min(acc_data_moving(:, 2)), max(acc_data_moving(:, 2)), n)', mean_acc_moving(2), cov_acc_moving(2)), 'r', 'linewidth', 1.5)
title('Accelerometer data y')
subplot(3,5,3)
histogram(acc_data_moving(:, 3), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(acc_data_moving(:, 3)), max(acc_data_moving(:, 3)), n), mvnpdf(linspace(min(acc_data_moving(:, 3)), max(acc_data_moving(:, 3)), n)', mean_acc_moving(3), cov_acc_moving(3)), 'r', 'linewidth', 1.5)
title('Accelerometer data z')

subplot(3,5,6)
histogram(gyro_data_moving(:, 1), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(gyro_data_moving(:, 1)), max(gyro_data_moving(:, 1)), n), mvnpdf(linspace(min(gyro_data_moving(:, 1)), max(gyro_data_moving(:, 1)), n)', mean_gyro_moving(1), cov_gyro_moving(1)), 'r', 'linewidth', 1.5)
title('Gyroscope data x')
subplot(3,5,7)
histogram(gyro_data_moving(:, 2), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(gyro_data_moving(:, 2)), max(gyro_data_moving(:, 2)), n), mvnpdf(linspace(min(gyro_data_moving(:, 2)), max(gyro_data_moving(:, 2)), n)', mean_gyro_moving(2), cov_gyro_moving(2)), 'r', 'linewidth', 1.5)
title('Gyroscope data y')
subplot(3,5,8)
histogram(gyro_data_moving(:, 3), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(gyro_data_moving(:, 3)), max(gyro_data_moving(:, 3)), n), mvnpdf(linspace(min(gyro_data_moving(:, 3)), max(gyro_data_moving(:, 3)), n)', mean_gyro_moving(3), cov_gyro_moving(3)), 'r', 'linewidth', 1.5)
title('Gyroscope data z')

subplot(3,5,11)
histogram(mag_data_moving(:, 1), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(mag_data_moving(:, 1)), max(mag_data_moving(:, 1)), n), mvnpdf(linspace(min(mag_data_moving(:, 1)), max(mag_data_moving(:, 1)), n)', mean_mag_moving(1), cov_mag_moving(1)), 'r', 'linewidth', 1.5)
title('Magnetometer data x')
subplot(3,5,12)
histogram(mag_data_moving(:, 2), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(mag_data_moving(:, 2)), max(mag_data_moving(:, 2)), n), mvnpdf(linspace(min(mag_data_moving(:, 2)), max(mag_data_moving(:, 2)), n)', mean_mag_moving(2), cov_mag_moving(2)), 'r', 'linewidth', 1.5)
title('Magnetometer data y')
subplot(3,5,13)
histogram(mag_data_moving(:, 3), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(mag_data_moving(:, 3)), max(mag_data_moving(:, 3)), n), mvnpdf(linspace(min(mag_data_moving(:, 3)), max(mag_data_moving(:, 3)), n)', mean_mag_moving(3), cov_mag_moving(3)), 'r', 'linewidth', 1.5)
title('Magnetometer data z')


%% Plot results still
figure('Position', [size(1)*0.1, size(2)*0.06, size(1)*0.7, size(2)*0.85]); 

subplot(3,5,4:5)
plot(t_still, acc_data_still)
legend('x','y','z')
title('Accelerometer'); xlabel('Time [s]'); ylabel('Sensor reading value')
axis('tight')
subplot(3,5,9:10)
plot(t_still, gyro_data_still)
legend('x','y','z')
title('Gyroscope'); xlabel('Time [s]'); ylabel('Sensor reading value')
axis('tight')
subplot(3,5,14:15)
plot(t_still, mag_data_still)
legend('x','y','z')
title('Magnetometer'); xlabel('Time [s]'); ylabel('Sensor reading value')
axis('tight')
titleObj = sgtitle('When the phone is still', 'fontsize', 20);

subplot(3,5,1)
histogram(acc_data_still(:, 1), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(acc_data_still(:, 1)), max(acc_data_still(:, 1)), n), mvnpdf(linspace(min(acc_data_still(:, 1)), max(acc_data_still(:, 1)), n)', mean_acc_still(1), cov_acc_still(1)), 'r', 'linewidth', 1.5)
title('Accelerometer data x')
subplot(3,5,2)
histogram(acc_data_still(:, 2), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(acc_data_still(:, 2)), max(acc_data_still(:, 2)), n), mvnpdf(linspace(min(acc_data_still(:, 2)), max(acc_data_still(:, 2)), n)', mean_acc_still(2), cov_acc_still(2)), 'r', 'linewidth', 1.5)
title('Accelerometer data y')
subplot(3,5,3)
histogram(acc_data_still(:, 3), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(acc_data_still(:, 3)), max(acc_data_still(:, 3)), n), mvnpdf(linspace(min(acc_data_still(:, 3)), max(acc_data_still(:, 3)), n)', mean_acc_still(3), cov_acc_still(3)), 'r', 'linewidth', 1.5)
title('Accelerometer data z')

subplot(3,5,6)
histogram(gyro_data_still(:, 1), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(gyro_data_still(:, 1)), max(gyro_data_still(:, 1)), n), mvnpdf(linspace(min(gyro_data_still(:, 1)), max(gyro_data_still(:, 1)), n)', mean_gyro_still(1), cov_gyro_still(1)), 'r', 'linewidth', 1.5)
title('Gyroscope data x')
subplot(3,5,7)
histogram(gyro_data_still(:, 2), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(gyro_data_still(:, 2)), max(gyro_data_still(:, 2)), n), mvnpdf(linspace(min(gyro_data_still(:, 2)), max(gyro_data_still(:, 2)), n)', mean_gyro_still(2), cov_gyro_still(2)), 'r', 'linewidth', 1.5)
title('Gyroscope data y')
subplot(3,5,8)
histogram(gyro_data_still(:, 3), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(gyro_data_still(:, 3)), max(gyro_data_still(:, 3)), n), mvnpdf(linspace(min(gyro_data_still(:, 3)), max(gyro_data_still(:, 3)), n)', mean_gyro_still(1), cov_gyro_still(3)), 'r', 'linewidth', 1.5)
title('Gyroscope data z')

subplot(3,5,11)
histogram(mag_data_still(:, 1), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(mag_data_still(:, 1)), max(mag_data_still(:, 1)), n), mvnpdf(linspace(min(mag_data_still(:, 1)), max(mag_data_still(:, 1)), n)', mean_mag_still(1), cov_mag_still(1)), 'r', 'linewidth', 1.5)
title('Magnetometer data x')
subplot(3,5,12)
histogram(mag_data_still(:, 2), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(mag_data_still(:, 2)), max(mag_data_still(:, 2)), n), mvnpdf(linspace(min(mag_data_still(:, 2)), max(mag_data_still(:, 2)), n)', mean_mag_still(2), cov_mag_still(2)), 'r', 'linewidth', 1.5)
title('Magnetometer data y')
subplot(3,5,13)
histogram(mag_data_still(:, 3), bins, 'Normalization', 'pdf')
hold on; grid on;
plot(linspace(min(mag_data_still(:, 3)), max(mag_data_still(:, 3)), n), mvnpdf(linspace(min(mag_data_still(:, 3)), max(mag_data_still(:, 3)), n)', mean_mag_still(3), cov_mag_still(3)), 'r', 'linewidth', 1.5)
title('Magnetometer data z')
