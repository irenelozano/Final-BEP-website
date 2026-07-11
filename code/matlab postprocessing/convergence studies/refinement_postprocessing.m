%% ========================================================================
%  TU/e Mechanical Engineering - BEP Grid Convergence Post-Processor
%  Two Stacked Plots with Aligned X-Axes
%  ========================================================================

clear; clc; close all;

csv_file = '../runtime_output/grid_convergence_20260614_220618.csv';
if ~isfile(csv_file)
    error('CSV file not found: %s', csv_file);
end

data_table = readtable(csv_file, 'VariableNamingRule', 'preserve');
N_values   = data_table.N;
Drag_FX_N  = data_table.Drag_Fx_N;
n_cases    = length(N_values);

% Convert to milli-Newtons (N)
% Drag_FX_N = Drag_FX_N * 1000; 

% Recalculate % change dynamically
percent_change = zeros(n_cases, 1);
percent_change(1) = NaN;

for i = 2:n_cases
    percent_change(i) = ((Drag_FX_N(i) - Drag_FX_N(i-1)) / Drag_FX_N(i-1)) * 100;
end

abs_changes = abs(percent_change);
drag_final = Drag_FX_N(end);

%% ------------------------------------------------------------------------
%  Figure generation: two stacked plots
%  ------------------------------------------------------------------------

fig1 = figure('Name', 'Grid Convergence Dashboard', ...
    'Position', [100, 100, 1000, 700], ...
    'Color', 'w');

t = tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

% -------------------------------------------------------------------------
% Top plot: Drag force
% -------------------------------------------------------------------------
ax1 = nexttile;
hold(ax1, 'on'); 
grid(ax1, 'on'); 
box(ax1, 'on');

h1 = plot(ax1, N_values, Drag_FX_N, 'o-', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'Color', [0.2 0.4 0.8], ...
    'MarkerFaceColor', [0.3010 0.7450 0.9330], ...
    'DisplayName', 'Drag force $F_x$');

ylabel(ax1, 'Drag force $|F_x|$ (N)', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

title(ax1, 'Spatial Grid Convergence: Worst Case ($d=10\,\mu$m, $v=700$ m/s)', ...
    'Interpreter', 'latex', ...
    'FontSize', 15, ...
    'FontWeight', 'bold');

legend(ax1, h1, ...
    'Location', 'best', ...
    'Interpreter', 'latex');

ax1.FontSize = 11;
ax1.XTick = N_values;

% -------------------------------------------------------------------------
% Bottom plot: Percentage change
% -------------------------------------------------------------------------
ax2 = nexttile;
hold(ax2, 'on'); 
grid(ax2, 'on'); 
box(ax2, 'on');

h2 = plot(ax2, N_values, abs_changes, 's-', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 9, ...
    'Color', [0.8500 0.3250 0.0980], ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'DisplayName', '$|\Delta F_x/F_x|$');

% Highlight points below the 2% threshold
converged_idx = abs_changes < 2.0;

h3 = plot(ax2, N_values(converged_idx), abs_changes(converged_idx), 's', ...
    'MarkerSize', 9, ...
    'MarkerFaceColor', [0.4660 0.6740 0.1880], ...
    'MarkerEdgeColor', [0.4660 0.6740 0.1880], ...
    'LineStyle', 'none', ...
    'DisplayName', 'Below 2\% threshold');

h4 = yline(ax2, 2.0, 'k:', ...
    'LineWidth', 2, ...
    'DisplayName', '2\% threshold');

ylabel(ax2, '$|\Delta F_x/F_x|$ between refinement steps [\%]', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

xlabel(ax2, 'Grid resolution $N$', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

legend(ax2, [h2, h3, h4], ...
    'Location', 'best', ...
    'Interpreter', 'latex');

ax2.FontSize = 11;
ax2.XTick = N_values;

% Make sure bottom y-axis starts from 0
valid_changes = abs_changes(~isnan(abs_changes));
if ~isempty(valid_changes)
    ylim(ax2, [0, max(2.5, max(valid_changes) * 1.2)]);
else
    ylim(ax2, [0, 2.5]);
end

% -------------------------------------------------------------------------
% Align x-axes
% -------------------------------------------------------------------------
linkaxes([ax1, ax2], 'x');
xlim(ax1, [min(N_values), max(N_values)]);