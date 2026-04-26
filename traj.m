%% Double Lane Change (DLC) Reference Path Generation
clear; close all; clc;

% Parameters from Table 1
dx1 = 25.00;    % [m]
dx2 = 21.95;    % [m]
dy1 = 4.05;     % [m]
dy2 = 5.70;     % [m]
Xs1 = 27.19;    % [m]
Xs2 = 56.46;    % [m]
alpha = 2.4;    % Shape parameter

% Define the longitudinal distance range
X = 0:0.1:100;  % [m]

% Calculate intermediate variables z1 and z2
z1 = (alpha/dx1) * (X - Xs1) - (alpha/2);
z2 = (alpha/dx2) * (X - Xs2) - (alpha/2);

% Calculate Reference Lateral Displacement (y_ref)
y_ref = (dy1/2) * (1 + tanh(z1)) - (dy2/2) * (1 + tanh(z2));

% Calculate Reference Yaw Angle (psi_ref)
% Note: Using dot operators for element-wise array operations
term1 = (dy1 * alpha) ./ (dx1 * cosh(z1).^2);
term2 = (dy2 * alpha) ./ (dx2 * cosh(z2).^2);
psi_ref = atan(term1 - term2);

%% Plotting Results Side-by-Side
figure('Name', 'DLC Reference Path', 'Color', 'w', 'Position', [100, 100, 1000, 400]);

% Plot 1: y_ref vs X
subplot(1, 2, 1);
plot(X, y_ref, 'b', 'LineWidth', 2);
xlabel('Longitudinal Distance X (m)');
ylabel('Lateral Displacement y_{ref} (m)');
title('Reference Lateral Displacement');
grid on;

% Plot 2: psi_ref vs X
subplot(1, 2, 2);
plot(X, rad2deg(psi_ref), 'b', 'LineWidth', 2); % Converted to degrees for better readability
xlabel('Longitudinal Distance X (m)');
ylabel('Yaw Angle \psi_{ref} (deg)');
title('Reference Yaw Angle');
grid on;
