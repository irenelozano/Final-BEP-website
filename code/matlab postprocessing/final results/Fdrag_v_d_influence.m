%% ========================================================================
%  TU/e Mechanical Engineering - BEP DSMC Matrix Results Visualizer
%  Section 5.1: Parametric Trends (Drag vs. Diameter & Velocity)
%  ========================================================================
clear; clc; close all;

% 1. Load Data
data_file = 'matrix_results_20260601_125031.csv';
if ~isfile(data_file)
    error('Matrix results CSV not found! Check file name and path.');
end

data = readtable(data_file, 'VariableNamingRule', 'preserve');

% Extract columns
D_vals = data.Diameter_m;
V_vals = data.Velocity_ms;
F_raw = abs(data.Drag_Fx);

% Background force correction: fit F₀ per diameter and subtract
[unique_D, ~, ~] = unique(D_vals);
F0_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:length(unique_D)
    idx = D_vals == unique_D(i);
    v_i = V_vals(idx);
    F_i = F_raw(idx);
    p = polyfit(v_i, F_i, 1);
    F0_map(unique_D(i)) = p(2);
end
F_vals = zeros(size(D_vals));
for i = 1:length(D_vals)
    F_vals(i) = F_raw(i) - F0_map(D_vals(i));
end

% Find unique values
unique_V = unique(V_vals);
unique_D = unique(D_vals);

% 2. Setup Figure
fig = figure('Name', 'Section 5.1: Parametric Trends', 'Position', [100, 100, 1100, 500], 'Color', 'w');

% Define a high-contrast colormap (Cool to Warm)
colors_V = jet(length(unique_V)); % Colors for velocity lines
colors_D = jet(length(unique_D)); % Colors for diameter lines

%% --- PANEL A: Drag Force vs. Droplet Diameter ---
ax1 = subplot(1, 2, 1);
hold on; grid on; box on;
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 11);

% Loop through each Velocity and plot Fx vs D
for i = 1:length(unique_V)
    v_current = unique_V(i);
    
    % Find rows for this velocity
    idx = (V_vals == v_current);
    d_subset = D_vals(idx);
    f_subset = F_vals(idx);
    
    % Sort by diameter to ensure smooth lines
    [d_subset, sort_idx] = sort(d_subset);
    f_subset = f_subset(sort_idx);
    
    % Plot
    plot(d_subset, f_subset, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
        'Color', colors_V(i,:), 'MarkerFaceColor', 'w', ...
        'DisplayName', sprintf('$v = %d$ m/s', v_current));
end

xlabel('Droplet Diameter $d$ (m)', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Aerodynamic Drag Force $F_x$ (N)', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('(a) Drag Force Scaling with Droplet Size', 'FontSize', 13);
leg1 = legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10);
title(leg1, 'Flow Velocity');

%% --- PANEL B: Drag Force vs. Flow Velocity ---
ax2 = subplot(1, 2, 2);
hold on; grid on; box on;
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 11);

% Loop through each Diameter and plot Fx vs V
for i = 1:length(unique_D)
    d_current = unique_D(i);
    
    % Find rows for this diameter
    idx = (D_vals == d_current);
    v_subset = V_vals(idx);
    f_subset = F_vals(idx);
    
    % Sort by velocity to ensure smooth lines
    [v_subset, sort_idx] = sort(v_subset);
    f_subset = f_subset(sort_idx);
    
    % Format diameter for legend (e.g., 100 nm, 1 um)
    if d_current < 1e-6
        d_label = sprintf('$%d$ nm', round(d_current*1e9));
    else
        d_label = sprintf('$%d$ \\mu m', round(d_current*1e6));
    end
    
    % Plot
    plot(v_subset, f_subset, '-s', 'LineWidth', 2, 'MarkerSize', 6, ...
        'Color', colors_D(i,:), 'MarkerFaceColor', 'w', ...
        'DisplayName', d_label);
end

xlabel('Flow Velocity $v_0$ (m/s)', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Aerodynamic Drag Force $F_x$ (N)', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
title('(b) Drag Force Scaling with Flow Speed', 'FontSize', 13);

% Adjust x-axis ticks for velocity to show specific values clearly
xticks([20 100 250 400 550 700]);
xticklabels({'20', '100', '250', '400', '550', '700'});

leg2 = legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10);
title(leg2, 'Droplet Diameter');

%% --- Final Formatting ---
sgtitle('\textbf{Parametric Mapping of Tin Droplet Aerodynamic Drag in 100 Pa $H_2$}', ...
    'Interpreter', 'latex', 'FontSize', 15);