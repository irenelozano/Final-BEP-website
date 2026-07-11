% =========================================================================
% plot_grid_conv.m  —  Spatial Grid Convergence Dashboard
% =========================================================================
%
% PURPOSE:
%   Reads the grid convergence CSV from runtime_output and generates a
%   2-panel figure:
%     Panel 1: Drag force Fx [pN] vs Grid Resolution N
%     Panel 2: Relative Error |ΔFx/Fx_fine| [%] vs N
%
% INPUT:
%   ../runtime_output/grid_convergence_20260614_220618.csv
%
% OUTPUTS:
%   post_processing/figures/grid_conv_dashboard.png
%   post_processing/figures/grid_conv_dashboard.fig
%
% USAGE:
%   Run from the post_processing folder:
%     >> plot_grid_conv
% =========================================================================

clearvars; close all; clc;

%% ---- 0. Paths -----------------------------------------------------------

csv_file = '../runtime_output/grid_convergence_20260614_220618.csv';

if ~isfile(csv_file)
    error('CSV file not found: %s', csv_file);
end

fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% ---- 1. Load data -------------------------------------------------------

T = readtable(csv_file, 'VariableNamingRule', 'preserve');

% Expected columns in this CSV:
% N, Total_Cells, Drag_Fx_N, Drag_Fx_pN, PPC_actual, Fnum_used, Percent_Change

N_vals      = T.N;
Fx_vals_N   = T.Drag_Fx_N;
Fx_vals_pN  = T.Drag_Fx_pN;
total_cells = T.Total_Cells;
ppc_actual  = T.PPC_actual;
fnum_used   = T.Fnum_used;

% Sort by ascending grid resolution
[N_vals, idx] = sort(N_vals, 'ascend');

Fx_vals_N   = Fx_vals_N(idx);
Fx_vals_pN  = Fx_vals_pN(idx);
total_cells = total_cells(idx);
ppc_actual  = ppc_actual(idx);
fnum_used   = fnum_used(idx);

%% ---- 2. Relative error calculation -------------------------------------

% Reference = finest grid resolution
Fx_ref_N = Fx_vals_N(end);

% Relative error with respect to finest solution [%]
rel_err = abs((Fx_vals_N - Fx_ref_N) ./ Fx_ref_N) * 100;

% Finest value is the reference, so its error is zero.
% For log-scale plotting, replace it with NaN.
rel_err_plot = rel_err;
rel_err_plot(end) = NaN;

% Step-to-step percentage change, calculated from Fx values
step_change = NaN(size(Fx_vals_N));
for i = 2:length(Fx_vals_N)
    step_change(i) = ((Fx_vals_N(i) - Fx_vals_N(i-1)) / Fx_vals_N(i-1)) * 100;
end

%% ---- 3. Estimate convergence order -------------------------------------

% Approximate cell size using N only.
% Since L cancels out in the slope, using dx ~ 1/N is enough for p estimate.
dx_relative = 1 ./ N_vals;

valid = ~isnan(rel_err_plot) & rel_err_plot > 0;

if sum(valid) >= 2
    log_dx  = log(dx_relative(valid));
    log_err = log(rel_err_plot(valid));

    p_fit   = polyfit(log_dx, log_err, 1);
    p_order = p_fit(1);
else
    p_order = NaN;
end

%% ---- 4. Richardson extrapolation ---------------------------------------

if length(N_vals) >= 2 && ~isnan(p_order)
    Fx_f = Fx_vals_N(end);
    Fx_c = Fx_vals_N(end-1);

    r = N_vals(end) / N_vals(end-1);

    % Avoid unstable Richardson extrapolation if p is too small or negative
    p_used = max(abs(p_order), 1);

    Fx_RE_N  = Fx_f + (Fx_f - Fx_c) / (r^p_used - 1);
    Fx_RE_pN = Fx_RE_N * 1e12;
else
    Fx_RE_N  = Fx_ref_N;
    Fx_RE_pN = Fx_ref_N * 1e12;
end

%% ---- 5. Print summary ---------------------------------------------------

fprintf('\n=== Grid Convergence Analysis ===\n');
fprintf('CSV file        : %s\n', csv_file);
fprintf('Finest N        : %d\n', N_vals(end));
fprintf('Fx finest       : %.4e N = %.4f pN\n', Fx_ref_N, Fx_ref_N * 1e12);
fprintf('Fx Richardson   : %.4e N = %.4f pN\n', Fx_RE_N, Fx_RE_pN);

if ~isnan(p_order)
    fprintf('Estimated order : p ≈ %.2f\n', p_order);
else
    fprintf('Estimated order : not enough valid data\n');
end

fprintf('\nStep-to-step changes:\n');
for i = 1:length(N_vals)
    if i == 1
        fprintf('N = %d: reference coarse case\n', N_vals(i));
    else
        fprintf('N = %d: %+6.3f %% from previous grid\n', N_vals(i), step_change(i));
    end
end

fprintf('\nRelative error against finest grid:\n');
for i = 1:length(N_vals)
    fprintf('N = %d: %.4f %%\n', N_vals(i), rel_err(i));
end

%% ---- 6. Plot ------------------------------------------------------------

fig = figure('Name', 'Grid Convergence Dashboard', ...
             'Position', [100 80 900 720], ...
             'Color', 'w');

palette = [0.12 0.47 0.71;    % blue
           0.89 0.10 0.11;    % red
           0.20 0.60 0.20];   % green

%% Panel 1: Fx vs N

ax1 = subplot(2,1,1);
hold(ax1, 'on');

plot(ax1, N_vals, Fx_vals_pN, 'o-', ...
     'Color', palette(1,:), ...
     'LineWidth', 2, ...
     'MarkerSize', 8, ...
     'MarkerFaceColor', palette(1,:), ...
     'DisplayName', 'DSMC F_x');

yline(ax1, Fx_RE_pN, '--', ...
      'Color', palette(2,:), ...
      'LineWidth', 1.6, ...
      'Label', sprintf('Richardson extrap. %.3f pN', Fx_RE_pN), ...
      'LabelVerticalAlignment', 'bottom', ...
      'DisplayName', 'Richardson extrapolation');

xlabel(ax1, 'Grid resolution N [cells per side]', 'FontSize', 12);
ylabel(ax1, 'Steady-state drag force F_x [pN]', 'FontSize', 12);

title(ax1, 'Grid Convergence: Drag Force vs Grid Resolution', ...
      'FontSize', 13);

legend(ax1, 'Location', 'best', 'FontSize', 10);
grid(ax1, 'on');
box(ax1, 'on');

set(ax1, ...
    'FontSize', 11, ...
    'XTick', N_vals);

%% Panel 2: Relative error vs N

ax2 = subplot(2,1,2);
hold(ax2, 'on');

plot_idx = find(~isnan(rel_err_plot));

semilogy(ax2, N_vals(plot_idx), rel_err_plot(plot_idx), 's-', ...
         'Color', palette(2,:), ...
         'LineWidth', 2, ...
         'MarkerSize', 8, ...
         'MarkerFaceColor', palette(2,:), ...
         'DisplayName', '|F_x - F_{x,fine}| / |F_{x,fine}|');

yline(ax2, 1.0, 'k--', ...
      'LineWidth', 1.4, ...
      'Label', '1% threshold', ...
      'LabelVerticalAlignment', 'bottom');

if ~isnan(p_order) && length(plot_idx) >= 2
    text(ax2, N_vals(plot_idx(2)), rel_err_plot(plot_idx(2)) * 1.3, ...
         sprintf('Estimated slope p ≈ %.2f', p_order), ...
         'FontSize', 11, ...
         'Color', palette(2,:));
end

xlabel(ax2, 'Grid resolution N [cells per side]', 'FontSize', 12);
ylabel(ax2, 'Relative error against finest grid [%]', 'FontSize', 12);

title(ax2, 'Grid Convergence Error', 'FontSize', 13);

legend(ax2, 'Location', 'best', 'FontSize', 10);
grid(ax2, 'on');
box(ax2, 'on');

set(ax2, ...
    'FontSize', 11, ...
    'XTick', N_vals, ...
    'YScale', 'log');

% Annotate first N below 1% error
below1 = find(rel_err < 1.0 & rel_err > 0, 1, 'first');

if ~isempty(below1)
    xline(ax2, N_vals(below1), ':', ...
          'Color', palette(3,:), ...
          'LineWidth', 1.4);

    text(ax2, N_vals(below1), rel_err(below1) * 0.7, ...
         sprintf('Converged below 1%%\\nN = %d', N_vals(below1)), ...
         'FontSize', 10, ...
         'Color', palette(3,:), ...
         'FontWeight', 'bold');
end

%% ---- 7. Save ------------------------------------------------------------

fig_path_png = fullfile(fig_dir, 'grid_conv_dashboard.png');
fig_path_fig = fullfile(fig_dir, 'grid_conv_dashboard.fig');

saveas(fig, fig_path_fig);
exportgraphics(fig, fig_path_png, 'Resolution', 300);

fprintf('\nFigure saved to:\n%s\n', fig_path_png);