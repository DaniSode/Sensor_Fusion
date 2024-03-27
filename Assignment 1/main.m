clear all; close all; clc

%% Task 1


% Given values
mu_q = [0; 10];
Sigma_q = [0.3 0; 0 8];
A = [1 0.5; 0 1];

% Compute the mean and covariance of z using Result 1b)
mu_z = A * mu_q;
Sigma_z = A * Sigma_q * A';

% Generate points for the sigma ellipses of q and z
q_points = sigmaEllipse2D(mu_q, Sigma_q);
z_points = sigmaEllipse2D(mu_z, Sigma_z);

% Plot the sigma ellipses for q and z in the same figure
figure
hold on
plot(q_points(1,:), q_points(2,:), 'b')
plot(z_points(1,:), z_points(2,:), 'r')
axis equal
xlabel('q_1')
ylabel('q_2')
legend('q', 'z')


%% Task 2


% Define x as a Gaussian random variable
mu_x = 0;
sigma_x = sqrt(2);

% Define the transformation function
f = @(x) 3*x;

% Analytical calculation of z distribution
mu_z = 0;
sigma_z = sqrt(18);
z_samples_analytical = normrnd(mu_z, sigma_z, [1, 5000]);

% Approximate calculation of z distribution using approxGaussianTransform
[mu_z_approx, sigma_z_approx, z_samples_approx] = approxGaussianTransform(mu_x, sigma_x, f,500);

% Plot the histograms of the analytical and approximated z distributions
figure
subplot(1,2,1)
histogram(z_samples_analytical,'Normalization', 'pdf')
hold on
x = linspace(-20, 20, 5000);
analytical_pdf = normpdf(x, mu_z, sigma_z);
plot(x, analytical_pdf, 'LineWidth', 2)
title('Analytical z Distribution')
legend('Histogram', 'Gaussian PDF')
hold off

subplot(1,2,2)
histogram(z_samples_approx, 'Normalization', 'pdf')
hold on
x = linspace(-30, 30, 5000);
approx_pdf = normpdf(x, mu_z_approx, sigma_z_approx);
plot(x, approx_pdf, 'LineWidth', 2)
title('Approximated z Distribution')
legend('Histogram', 'Gaussian PDF')



% Define x as a Gaussian random variable
mu_x = 0;
sigma_x = sqrt(2);

% Define the transformation function
f = @(x) x^3;

% Analytical calculation of z distribution
mu_z = 0;
sigma_z = sqrt(15);
z_samples_analytical = normrnd(mu_z, sigma_z, [1, 5000]);

% Approximate calculation of z distribution using approxGaussianTransform
[mu_z_approx, sigma_z_approx, z_samples_approx] = approxGaussianTransform(mu_x, sigma_x, f,500);

% Plot the histograms of the analytical and approximated z distributions
figure
subplot(1,2,1)
histogram(z_samples_analytical,'Normalization', 'pdf')
hold on
x = linspace(-20, 20, 5000);
analytical_pdf = normpdf(x, mu_z, sigma_z);
plot(x, analytical_pdf, 'LineWidth', 2)
title('Analytical z Distribution')
legend('Histogram', 'Gaussian PDF')
hold off

subplot(1,2,2)
histogram(z_samples_approx, 'Normalization', 'pdf')
hold on
x = linspace(-30, 30, 5000);
approx_pdf = normpdf(x, mu_z_approx, sigma_z_approx);
plot(x, approx_pdf, 'LineWidth', 2)
title('Approximated z Distribution')
legend('Histogram', 'Gaussian PDF')

%% 3

% Set parameters
n = 10000; % number of samples
sigma_r = 1; % variance of r
a = -10; % lower bound of x
b = 10; % upper bound of x

% Define x, r and h(x)
x = a+(b-a).*rand(n, 1); % Uniform distribution between -1 and 1
r = sigma_r*randn(n,1); % normally distributed with mean 0 and variance sigma_r^2
hxlin = 2*x;
hxnon = x.^2; 

% Define y
ylin = hxlin + r;
ynon = hxnon + r;

% Plot a)
figure
subplot(2, 1, 1);
histogram(ylin, 100, 'Normalization', 'pdf');
edges = linspace(min(ylin), max(ylin), 100);
pdf_values = unifpdf(edges, min(ylin), max(ylin));
hold on
plot(edges, pdf_values,'r', 'LineWidth', 1);
plot([min(edges) min(edges)], [0 min(pdf_values)],'r', 'LineWidth', 1);
plot([max(edges) max(edges)], [0 max(pdf_values)], 'r','LineWidth', 1);
title('Linear, p(y)');

subplot(2, 1, 2);
histogram(ynon, 100, 'Normalization', 'pdf');
edges = linspace(min(ynon), max(ynon), 100);
pdf_values = unifpdf(edges, min(ynon), max(ynon));
hold on
plot(edges, pdf_values, 'r','LineWidth', 1);
plot([min(edges) min(edges)], [0 min(pdf_values)],'r', 'LineWidth', 1);
plot([max(edges) max(edges)], [0 max(pdf_values)], 'r','LineWidth', 1);
title('Non-linear, p(y)');


% Redefine x, r and h(x)
x = 1; % Fixed value 3
r = sigma_r*randn(n,1); % normally distributed with mean 0 and variance sigma_r^2
hxlin = 2*x;
hxnon = x.^2; 

% Redefine y
ylin = hxlin + r;
ynon = hxnon + r;

% Plot b)

figure
subplot(2, 1, 1);
histogram(ylin, 100, 'Normalization', 'pdf');
edges = linspace(min(ylin), max(ylin), 100);
pdf_values = normpdf(edges, mean(ylin), std(ylin));
hold on
plot(edges, pdf_values,'r', 'LineWidth', 1);
title('Linear, p(y|x)');

subplot(2, 1, 2);
histogram(ynon, 100, 'Normalization', 'pdf');
edges = linspace(min(ynon), max(ynon), 100);
pdf_values = normpdf(edges, mean(ynon), std(ynon));
hold on
plot(edges, pdf_values,'r', 'LineWidth', 1);
title('Non-linear, p(y|x)');

% d)

% Set parameters
n = 10000; % number of samples
sigma_r = 1; % variance of r
mu_x = 0;
sigma_x = 10;
x = normrnd(mu_x, sigma_x, [1, n]);
r = normrnd(0, sigma_r, [1, n]);
hxlin = 2*x;
hxnon = x.^2; 

% Redefine y
ylin = hxlin + r;
ynon = hxnon + r;

% Plot result
figure
subplot(2, 1, 1);
histogram(ylin, 100, 'Normalization', 'pdf');
edges = linspace(min(ylin), max(ylin), 100);
pdf_values = normpdf(edges, mean(ylin), std(ylin));
hold on
plot(edges, pdf_values,'r', 'LineWidth', 1);
title('Linear, p(y)');

subplot(2, 1, 2);
histogram(ynon, 100, 'Normalization', 'pdf');
edges = linspace(min(ynon), max(ynon), 100);
pdf_values = normpdf(edges, mean(ynon), std(ynon));
hold on
plot(edges, pdf_values,'r', 'LineWidth', 1);
title('Non-linear, p(y)');

% d)

% Set parameters
n = 10000; % number of samples
sigma_r = 1; % variance of r
mu_x = 0;
sigma_x = 10;
x = 1;
r = normrnd(0, sigma_r, [1, n]);
hxlin = 2*x;
hxnon = x.^2; 

% Redefine y
ylin = hxlin + r;
ynon = hxnon + r;

figure
subplot(2, 1, 1);
histogram(ylin, 100, 'Normalization', 'pdf');
edges = linspace(min(ylin), max(ylin), 100);
pdf_values = normpdf(edges, mean(ylin), std(ylin));
hold on
plot(edges, pdf_values,'r', 'LineWidth', 1);
title('Linear, p(y|x)');

subplot(2, 1, 2);
histogram(ynon, 100, 'Normalization', 'pdf');
edges = linspace(min(ynon), max(ynon), 100);
pdf_values = normpdf(edges, mean(ynon), std(ynon));
hold on
plot(edges, pdf_values,'r', 'LineWidth', 1);
title('Non-linear, p(y|x)');

%% 4

% Defining given
p_theta_minus1 = 0.5;
p_theta_1 = 0.5;
sigma2 = 0.5^2;

% Defining sample
num_samples = 10000;
y = zeros(num_samples, 1);

% Looping to create a set of random varaibles creating the distribution we
% want
for i = 1:num_samples
    if rand < p_theta_minus1
        theta = -1;
    else
        theta = 1;
    end
    w = randn * sqrt(sigma2);
    y(i) = theta + w;
end

% Plot
figure
histogram(y);

% d)

% Define the parameters
sigma = 0.5;
theta_values = [-1, 1];
prior_prob = [0.5, 0.5];
y = 0.7;

% Compute the conditional probability densities
p_y_given_theta1 = (1/sqrt(2*pi*sigma^2)) * exp(-(y-1)^2/(2*sigma^2));
p_y_given_theta2 = (1/sqrt(2*pi*sigma^2)) * exp(-(y+1)^2/(2*sigma^2));

% Compute the marginal probability density of y
p_y = sum(prior_prob .* [p_y_given_theta1, p_y_given_theta2]);

% Compute the posterior probabilities of theta
posterior_prob = [p_y_given_theta1 * prior_prob(1), p_y_given_theta2 * prior_prob(2)] / p_y;

% Find the most likely value of theta
[~, idx] = max(posterior_prob);
theta_hat = theta_values(idx);

disp(['Most likely value of theta: ', num2str(theta_hat)])

% Compute the MMSE estimator
theta_hat_MMSE = -1*posterior_prob(1) + 1*posterior_prob(2);

% Display the MMSE estimator
disp(['MMSE estimator = ', num2str(theta_hat_MMSE)])