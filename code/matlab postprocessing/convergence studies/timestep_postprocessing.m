%% ========================================================================
%  TU/e Mechanical Engineering - BEP Timestep Convergence Post-Processor
%  2-Graph Thesis Format
%
%  Graph 1:
%    - Left y-axis: Drag force |Fx| [N]
%    - Right y-axis: Step-to-step percentage change [%]
%
%  Graph 2:
%    - Computational workload tracker
%
%  Important logic:
%  - Data is sorted from coarse timestep to fine timestep for convergence.
%  - Drag force is kept in Newtons.
%  - Percentage change is plotted at the refined timestep value.
%  - Selected timestep is the first point after which all later refinements
%    remain below the 2% step-to-step threshold.
%  ========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%  1. Load latest timestep convergence CSV
%  ------------------------------------------------------------------------

output_dir = '../runtime_output/';
file_pattern = fullfile(output_dir, 'time_convergence_*.csv');

dir_info = dir(file_pattern);

if isempty(dir_info)
    error('No time convergence spreadsheet found! Run the timestep refinement pipeline first.');
end

% Sort by date to always pick the most recent file
[~, idx] = sort([dir_info.datenum], 'descend');
target_file = fullfile(output_dir, dir_info(idx(1)).name);

fprintf('Loading latest timestep convergence file: %s\n', dir_info(idx(1)).name);

time_table = readtable(target_file, 'VariableNamingRule', 'preserve');

%% ------------------------------------------------------------------------
%  2. Extract data
%  ------------------------------------------------------------------------

dt_raw           = time_table.Timestep_s;
total_steps_raw  = time_table.Total_Steps;
drag_forces_raw  = abs(time_table.Drag_Fx);   % Force kept in Newtons [N]

% Filter missing or failed output points
valid = drag_forces_raw > 0;

if any(~valid)
    fprintf('WARNING: %d point(s) with zero drag excluded.\n', sum(~valid));
    dt_raw           = dt_raw(valid);
    total_steps_raw  = total_steps_raw(valid);
    drag_forces_raw  = drag_forces_raw(valid);
end

% Sort from smallest to largest timestep for normal x-axis plotting
[dt_plot, plot_idx] = sort(dt_raw, 'ascend');
steps_plot          = total_steps_raw(plot_idx);
drag_plot_N         = drag_forces_raw(plot_idx);

% Sort from largest to smallest timestep for refinement logic
% This means refinement goes from coarse dt -> fine dt.
[dt_refine, refine_idx] = sort(dt_raw, 'descend');
steps_refine            = total_steps_raw(refine_idx);
drag_refine_N           = drag_forces_raw(refine_idx);

n_cases = length(dt_refine);

%% ------------------------------------------------------------------------
%  3. Compute convergence metrics
%  ------------------------------------------------------------------------

% Final value is the finest timestep case
dt_finest = min(dt_raw);
drag_final_N = drag_plot_N(1);   % because dt_plot is ascending

% Step-by-step percentage change in the physically correct refinement order:
% coarse dt -> fine dt
percent_change = NaN(n_cases, 1);

for i = 2:n_cases
    percent_change(i) = ((drag_refine_N(i) - drag_refine_N(i-1)) / drag_refine_N(i-1)) * 100;
end

abs_percent_change = abs(percent_change);

% For plotting percentage change:
% Each change is assigned to the coarser timestep value so the plot reaches
% the rightmost timestep on the x-axis.
change_dt = dt_refine(1:end-1);
change_values = abs_percent_change(2:end);

% Sort change points for normal numerical plotting
[change_dt_plot, change_sort_idx] = sort(change_dt, 'ascend');
change_values_plot = change_values(change_sort_idx);

% Find selected timestep:
% first timestep in the coarse-to-fine sequence after which all following
% refinement changes are below 2%.
idx_selected = n_cases;   % default to finest timestep

for i = 2:n_cases
    if all(abs(percent_change(i:end)) < 2.0)
        idx_selected = i;
        break;
    end
end

dt_selected = dt_refine(idx_selected);
drag_selected_N = drag_refine_N(idx_selected);

% For report summary
drag_min_N = min(drag_plot_N);
drag_max_N = max(drag_plot_N);
drag_spread_percent = ((drag_max_N - drag_min_N) / drag_final_N) * 100;

%% ------------------------------------------------------------------------
%  4. Figure generation: two stacked plots with aligned x-axes
%  ------------------------------------------------------------------------

fig1 = figure( ...
    'Name', 'Timestep Convergence Dashboard', ...
    'NumberTitle', 'off', ...
    'Position', [100, 100, 1000, 700], ...
    'Color', 'w');

tiledlayout(2, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

dt_CFL = 4.8e-9;

%% ------------------------------------------------------------------------
%  Top plot: Drag force
%  ------------------------------------------------------------------------

ax1 = nexttile;
hold(ax1, 'on'); 
grid(ax1, 'on'); 
box(ax1, 'on');

h1 = semilogx(ax1, dt_plot, drag_plot_N, '-o', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 9, ...
    'Color', [0.2 0.4 0.8], ...
    'MarkerFaceColor', [0.3010 0.7450 0.9330], ...
    'DisplayName', 'Drag force $|F_x|$');

h2 = xline(ax1, dt_CFL, 'k-.', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'CFL limit');

ylabel(ax1, 'Drag force $|F_x|$ [N]', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

title(ax1, '(a) Drag Force', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);

legend(ax1, [h1, h2], ...
    'Location', 'best', ...
    'Interpreter', 'latex', ...
    'FontSize', 9);

ax1.FontSize = 11;

%% ------------------------------------------------------------------------
%  Bottom plot: Step-to-step percentage change
%  ------------------------------------------------------------------------

ax2 = nexttile;
hold(ax2, 'on'); 
grid(ax2, 'on'); 
box(ax2, 'on');

h3 = semilogx(ax2, change_dt_plot, change_values_plot, 's-', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 9, ...
    'Color', [0.8500 0.3250 0.0980], ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'DisplayName', '$|\Delta F_x/F_x|$');

% Highlight points below 2%
below_threshold = change_values_plot < 2.0;

h4 = semilogx(ax2, change_dt_plot(below_threshold), change_values_plot(below_threshold), 's', ...
    'MarkerSize', 9, ...
    'MarkerFaceColor', [0.4660 0.6740 0.1880], ...
    'MarkerEdgeColor', [0.4660 0.6740 0.1880], ...
    'LineStyle', 'none', ...
    'DisplayName', 'Below 2\% threshold');

h5 = yline(ax2, 2.0, 'k:', ...
    'LineWidth', 2, ...
    'DisplayName', '2\% threshold');

h6 = xline(ax2, dt_CFL, 'k-.', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'CFL limit');

ylabel(ax2, '$|\Delta F_x/F_x|$ between refinement steps [\%]', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

xlabel(ax2, 'Timestep $\Delta t$ [s]', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

title(ax2, '(b) Step-to-Step Convergence Criterion', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);

legend(ax2, [h3, h4, h5, h6], ...
    'Location', 'best', ...
    'Interpreter', 'latex', ...
    'FontSize', 9);

ax2.FontSize = 11;

if isempty(change_values_plot)
    ylim(ax2, [0, 2.5]);
else
    ylim(ax2, [0, max(2.5, max(change_values_plot) * 1.2)]);
end

%% ------------------------------------------------------------------------
%  Align x-axes
%  ------------------------------------------------------------------------

linkaxes([ax1, ax2], 'x');
xlim(ax1, [min(dt_plot), max(dt_plot)]);

%% ------------------------------------------------------------------------
%  Overall title
%  ------------------------------------------------------------------------

sgtitle('Temporal Convergence Study: Worst Case ($N=40$, $d=10~\mu$m, $v=700~\mathrm{m/s}$)', ...
    'Interpreter', 'latex', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');%% ------------------------------------------------------------------------
%  5. Export figure
%  ------------------------------------------------------------------------

exportgraphics(fig1, '../post_processing/timestep_convergence_2graphs.png', ...
    'Resolution', 300);

%% ------------------------------------------------------------------------
%  6. Console summary
%  ------------------------------------------------------------------------

fprintf('\n==========================================================\n');
fprintf(' TIMESTEP CONVERGENCE ANALYSIS SUMMARY\n');
fprintf('==========================================================\n\n');

fprintf('Input file: %s\n', dir_info(idx(1)).name);
fprintf('Number of timestep levels: %d\n', n_cases);

fprintf('\nDrag force range:\n');
fprintf('  Minimum:       %.4e N\n', drag_min_N);
fprintf('  Maximum:       %.4e N\n', drag_max_N);
fprintf('  Finest dt:     %.2e s\n', dt_finest);
fprintf('  Finest value:  %.4e N\n', drag_final_N);
fprintf('  Selected dt:   %.2e s\n', dt_selected);
fprintf('  Selected drag: %.4e N\n', drag_selected_N);
fprintf('  Total spread:  %.4f %% of finest-dt value\n\n', drag_spread_percent);

fprintf('Step-by-step refinement changes, coarse dt -> fine dt:\n');

for i = 2:n_cases
    fprintf('  dt %.2e -> %.2e: %.4f %% change\n', ...
        dt_refine(i-1), dt_refine(i), percent_change(i));
end

fprintf('\nCONVERGENCE STATUS:\n');

if n_cases >= 2
    last_change = abs(percent_change(end));

    fprintf('  Last refinement: dt %.2e -> %.2e gives %.4f %% change\n', ...
        dt_refine(end-1), dt_refine(end), last_change);

    if all(abs(percent_change(idx_selected:end)) < 2.0)
        fprintf('  CONVERGED: from dt = %.2e s onward, all finer refinements are below 2%%.\n', ...
            dt_selected);
    else
        fprintf('  NOT CONVERGED: final refinement sequence does not remain below 2%%.\n');
    end
else
    fprintf('  Not enough timestep levels to assess convergence.\n');
end

fprintf('\n==========================================================\n\n');