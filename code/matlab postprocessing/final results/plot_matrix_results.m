%% ========================================================================
%  TU/e & ASML BEP - 30-Case Matrix Post-Processing
%  Plots Aerodynamic Drag Force and Drag Coefficient (Cd) vs Velocity
%  ========================================================================
clear; clc; close all;

% 1. Load the Production Data
filename = 'matrix_results_20260601_125031.csv';
data = readtable(filename);

% Sort data just to be safe
data = sortrows(data, {'Diameter_m', 'Velocity_ms'});

% Background force correction: fit F₀ per diameter and subtract
[unique_D, ~, ~] = unique(data.Diameter_m);
F0_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:length(unique_D)
    idx = data.Diameter_m == unique_D(i);
    v_i = data.Velocity_ms(idx);
    F_i = abs(data.Drag_Fx(idx));
    p = polyfit(v_i, F_i, 1);
    F0_map(unique_D(i)) = p(2);
end
F_corrected = zeros(height(data), 1);
for i = 1:height(data)
    F_corrected(i) = abs(data.Drag_Fx(i)) - F0_map(data.Diameter_m(i));
end
data.Drag_Fx = F_corrected;

% Filter out zero-force numerical artifacts before any processing
data(data.Drag_Fx <= 0, :) = [];

% Extract unique dimensions for the loop
diameters = unique(data.Diameter_m);
velocities = unique(data.Velocity_ms);

% Define H2 Gas Constants (100 Pa, 293.15 K) for Cd calculation
n_inf = 2.47e22;                % Number density [m^-3]
m_H2 = 3.346e-27;               % Mass of H2 molecule [kg]
rho = n_inf * m_H2;             % Gas density [kg/m^3] (~8.26e-5 kg/m^3)

% Setup visual formatting arrays
colors = lines(length(diameters));
markers = {'o', 's', '^', 'd', 'v'};

%% FIGURE 1: Absolute Drag Force vs Velocity (Log-Scale)
fig1 = figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; grid on;

for i = 1:length(diameters)
    d = diameters(i);
    % Extract the subset of data for this specific diameter
    idx = (data.Diameter_m == d);
    V_subset = data.Velocity_ms(idx);
    F_subset = data.Drag_Fx(idx);
    
    % Plot using logarithmic Y-axis because forces span multiple magnitudes
    semilogy(V_subset, F_subset, ['-', markers{i}], 'LineWidth', 2, ...
        'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 8, ...
        'DisplayName', sprintf('d = %g \\mu m', d * 1e6));
end

% Thesis-grade formatting
set(gca, 'FontSize', 12, 'YMinorGrid', 'on');
xlabel('Inflow Velocity $v$ (m/s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Time-Averaged Drag Force $|F_x|$ (N)', 'Interpreter', 'latex', 'FontSize', 14);
title('\bf{Aerodynamic Drag Force on Tin Droplets in H_2 Vacuum (100 Pa)}', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 11);


%% FIGURE 2: Drag Coefficient (Cd) vs Velocity
fig2 = figure('Color', 'w', 'Position', [150, 150, 800, 500]);
hold on; grid on;

for i = 1:length(diameters)
    d = diameters(i);
    idx = (data.Diameter_m == d);
    V_subset = data.Velocity_ms(idx);
    F_subset = data.Drag_Fx(idx);
    
    % Calculate Cross-Sectional Area
    A = pi * (d / 2)^2;
    
    % Calculate Drag Coefficient: Cd = (2 * F) / (rho * V^2 * A)
    Cd = (2 .* F_subset) ./ (rho .* (V_subset.^2) .* A);
    
    plot(V_subset, Cd, ['-', markers{i}], 'LineWidth', 2, ...
        'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 8, ...
        'DisplayName', sprintf('d = %g \\mu m', d * 1e6));
end

% Thesis-grade formatting
set(gca, 'FontSize', 12);
xlabel('Inflow Velocity $v$ (m/s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Drag Coefficient $C_d$ (-)', 'Interpreter', 'latex', 'FontSize', 14);
title('\bf{Droplet Drag Coefficient ($C_d$) Across Velocity Regimes}', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 11);