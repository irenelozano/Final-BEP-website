%% ========================================================================
%  ASML Debris Mitigation - Smooth Penetration Contour Map
%  Display-only interpolation in log(d)
%  ========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%  1. Load production data
%  ------------------------------------------------------------------------

filename = 'matrix_results_20260601_125031.csv';
data = readtable(filename);

diameters        = unique(data.Diameter_m);
velocities_data  = unique(data.Velocity_ms);

rho_sn = 6990; % Density of liquid tin [kg/m^3]

% Reference distance planes [m]
critical_distance_m  = 1.30;
reference_distance_m = 0.48;

% Plotting limit for colour scale [m]
plot_limit_m = 3.5;

%% ------------------------------------------------------------------------
%  2. Background force correction: F0 subtraction
%  ------------------------------------------------------------------------

F0_map = containers.Map('KeyType', 'double', 'ValueType', 'double');

fprintf('--- Background force (F0) correction ---\n');

for i = 1:length(diameters)

    idx = data.Diameter_m == diameters(i);

    v_i = data.Velocity_ms(idx);
    F_i = data.Drag_Fx(idx);

    % Linear fit to estimate residual zero-velocity force offset
    p = polyfit(v_i, F_i, 1);
    F0_map(diameters(i)) = p(2);

    if diameters(i) < 1e-6
        label = sprintf('%.0f nm', diameters(i)*1e9);
    else
        label = sprintf('%.0f um', diameters(i)*1e6);
    end

    fprintf('  d = %s: F0 = %.4e N\n', label, F0_map(diameters(i)));
end

fprintf('---\n\n');

% Apply F0 correction
F0_vals = zeros(size(data.Diameter_m));

for k = 1:length(diameters)
    F0_vals(data.Diameter_m == diameters(k)) = F0_map(diameters(k));
end

data.F_Corrected = data.Drag_Fx - F0_vals;

%% ------------------------------------------------------------------------
%  3. Fit drag law at simulated diameters
%  ------------------------------------------------------------------------

fprintf('Fitting drag law F = a v for each diameter...\n');

a_fit_by_d = zeros(size(diameters));
xstop_slope_by_d = zeros(size(diameters)); % x_stop / v0 = m/a

fprintf('--- Fitted drag law: F = a v ---\n');

for j = 1:length(diameters)

    d = diameters(j);

    % Droplet mass
    mass = rho_sn * (4/3) * pi * (d/2)^3;

    % Extract corrected drag data for this diameter
    idx = data.Diameter_m == d;

    V_data = data.Velocity_ms(idx);
    F_data = data.F_Corrected(idx);

    [V_data, sort_idx] = sort(V_data);
    F_data = F_data(sort_idx);

    % Fit through origin: F = a*v
    a_fit = V_data(:) \ F_data(:);

    % Avoid non-physical or zero fitted drag coefficients
    if a_fit <= 0
        warning('Non-positive fitted drag coefficient for d = %.3e m. Replacing with eps.', d);
        a_fit = eps;
    end

    a_fit_by_d(j) = a_fit;
    xstop_slope_by_d(j) = mass / a_fit;

    if d < 1e-6
        label = sprintf('%.0f nm', d*1e9);
    else
        label = sprintf('%.0f um', d*1e6);
    end

    fprintf('  d = %s: a = %.4e N/(m/s), x_stop/v0 = %.4e m/(m/s)\n', ...
        label, a_fit_by_d(j), xstop_slope_by_d(j));
end

fprintf('---\n\n');

%% ------------------------------------------------------------------------
%  4. Build smooth display grid
%  ------------------------------------------------------------------------

fprintf('Building smooth display grid...\n');

v_max = max(velocities_data);

% Dense velocity grid including v0 -> 0
velocities_plot = unique([ ...
    0, 0.1, 1, 2, 5, 10, 20, ...
    linspace(25, v_max, 200) ...
]);

% Dense diameter grid in log-space for display only
diameters_plot = logspace( ...
    log10(min(diameters)), ...
    log10(max(diameters)), ...
    300 ...
);

% Interpolate x_stop/v0 in log(d)
logd_data = log10(diameters);
logd_plot = log10(diameters_plot);

xstop_slope_plot = interp1( ...
    logd_data, ...
    xstop_slope_by_d, ...
    logd_plot, ...
    'pchip' ...
);

% Avoid small interpolation artifacts
xstop_slope_plot = max(xstop_slope_plot, 0);

% Smooth plotting grid
[D_plot_grid, V_plot_grid] = meshgrid(diameters_plot, velocities_plot);

% Analytical stopping distance on smooth grid
Travel_Distance_m = V_plot_grid .* reshape(xstop_slope_plot, 1, []);

% Cap plotted travel distance for readable colour scale
Travel_Distance_plot_m = min(Travel_Distance_m, plot_limit_m);

%% ------------------------------------------------------------------------
%  5. Plot smooth contour map
%  ------------------------------------------------------------------------

figure('Color', 'w', 'Position', [150, 150, 850, 600]);
hold on;

levels = linspace(0, plot_limit_m, 31);

contourf(D_plot_grid * 1e6, V_plot_grid, Travel_Distance_plot_m, ...
    levels, ...
    'LineColor', 'none');

set(gca, 'XScale', 'log', 'FontSize', 12);

xlim([min(diameters)*1e6, max(diameters)*1e6]);
ylim([0, v_max]);

xticks([0.1 0.5 1 5 10]);
xticklabels({'0.1', '0.5', '1', '5', '10'});

colormap(turbo);

c = colorbar;
c.Label.String = 'Total travel distance (m)';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 13;

caxis([0 plot_limit_m]);

c.Ticks = 0:0.5:plot_limit_m;
c.TickLabels = string(c.Ticks);

%% ------------------------------------------------------------------------
%  6. Reference distance contours
%  ------------------------------------------------------------------------

% 1.30 m reference distance
contour(D_plot_grid * 1e6, V_plot_grid, Travel_Distance_m, ...
    [critical_distance_m critical_distance_m], ...
    'LineColor', [1 1 1]*0.6, ...
    'LineWidth', 2.8, ...
    'LineStyle', '--');

% 0.48 m reference distance
contour(D_plot_grid * 1e6, V_plot_grid, Travel_Distance_m, ...
    [reference_distance_m reference_distance_m], ...
    'LineColor', [1 1 1]*0.6, ...
    'LineWidth', 2.0, ...
    'LineStyle', ':');

% Legend handles only, no labels printed on the contour lines
h_13 = plot(nan, nan, '--', ...
    'Color', [1 1 1]*0.6, ...
    'LineWidth', 2.8, ...
    'DisplayName', '$s = 1.30$ m');

h_048 = plot(nan, nan, ':', ...
    'Color', [1 1 1]*0.6, ...
    'LineWidth', 2.0, ...
    'DisplayName', '$s = 0.48$ m');

legend([h_13, h_048], ...
    'Location', 'northwest', ...
    'Interpreter', 'latex', ...
    'FontSize', 10);

%% ------------------------------------------------------------------------
%  7. Formatting
%  ------------------------------------------------------------------------

xlabel('Droplet diameter $d$ ($\mu$m)', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);

ylabel('Initial ejection velocity $v_0$ (m/s)', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);

title('Debris Penetration Map: Travel Distance in 100 Pa $H_2$', ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'FontWeight', 'bold');

box on;
hold off;