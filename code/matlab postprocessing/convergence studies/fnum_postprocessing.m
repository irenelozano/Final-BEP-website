%% ========================================================================
%  TU/e Mechanical Engineering - BEP Particle Convergence Post-Processor
%  Single-Panel Thesis Format with Two Y-Axes
%  ========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%  1. Load latest fnum convergence CSV
%  ------------------------------------------------------------------------

output_dir = '../runtime_output/';
file_pattern = fullfile(output_dir, 'fnum_convergence_*.csv');

dir_info = dir(file_pattern);

if isempty(dir_info)
    error('No fnum convergence spreadsheet found! Run the fnum refinement pipeline first.');
end

% Sort by date to always pick the most recent file
[~, idx] = sort([dir_info.datenum], 'descend');
target_file = fullfile(output_dir, dir_info(idx(1)).name);

fprintf('Loading latest fnum convergence file: %s\n', dir_info(idx(1)).name);

fnum_table = readtable(target_file, 'VariableNamingRule', 'preserve');

%% ------------------------------------------------------------------------
%  2. Extract and sort data
%  ------------------------------------------------------------------------

fnum_raw = fnum_table.Fnum;
total_particles_raw = fnum_table.Total_Simulated_Particles;
drag_forces_raw = abs(fnum_table.Drag_Fx);   % Force kept in Newtons [N]

% Grid is fixed at N = 40 for the fnum convergence study
grid_res = 40;
grid_cells = grid_res^3;

% Calculate particles per cell
PPC_raw = ceil(total_particles_raw ./ grid_cells);

% Sort from lowest to highest PPC
[PPC_values, sort_idx] = sort(PPC_raw, 'ascend');
fnum_values = fnum_raw(sort_idx);
total_particles = total_particles_raw(sort_idx);
drag_forces_N = drag_forces_raw(sort_idx);

n_cases = length(PPC_values);

%% ------------------------------------------------------------------------
%  3. Compute convergence metrics
%  ------------------------------------------------------------------------

% Final value is the highest particle-density case
drag_final = drag_forces_N(end);

% Step-by-step percentage change
percent_change = zeros(n_cases, 1);
percent_change(1) = NaN;

for i = 2:n_cases
    percent_change(i) = ((drag_forces_N(i) - drag_forces_N(i-1)) / drag_forces_N(i-1)) * 100;
end

abs_changes = abs(percent_change);

% Find first case within 2% of final value
percent_deviation_from_final = ((drag_forces_N - drag_final) / drag_final) * 100;
idx_converged = find(abs(percent_deviation_from_final) <= 2.0, 1, 'first');

if isempty(idx_converged)
    idx_converged = n_cases;
end

PPC_converged = PPC_values(idx_converged);
drag_converged = drag_forces_N(idx_converged);

%% ------------------------------------------------------------------------
%  4. Figure generation: two stacked plots with aligned x-axes
%  ------------------------------------------------------------------------

fig1 = figure( ...
    'Name', 'Particle Convergence Dashboard', ...
    'NumberTitle', 'off', ...
    'Position', [100, 100, 1000, 700], ...
    'Color', 'w');

t = tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Top plot: Drag force

ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');

h1 = plot(ax1, PPC_values, drag_forces_N, 'o-', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'Color', [0.2 0.4 0.8], ...
    'MarkerFaceColor', [0.3010 0.7450 0.9330], ...
    'DisplayName', 'Drag force $|F_x|$');

h2 = xline(ax1, 20, 'k--', ...
    'LineWidth', 1.5, ...
    'DisplayName', '20 PPC reference');

ylabel(ax1, 'Drag force $|F_x|$ [N]', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

title(ax1, 'Statistical Particle Convergence: Worst Case ($N=40$, $d=10~\mu$m, $v=700~\mathrm{m/s}$)', ...
    'Interpreter', 'latex', ...
    'FontSize', 15, ...
    'FontWeight', 'bold');

legend(ax1, [h1, h2], ...
    'Location', 'best', ...
    'Interpreter', 'latex', ...
    'FontSize', 9);

ax1.FontSize = 11;
ax1.XTick = PPC_values;

%% Bottom plot: Percentage change

ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

h3 = plot(ax2, PPC_values, abs_changes, 's-', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 9, ...
    'Color', [0.8500 0.3250 0.0980], ...
    'MarkerFaceColor', [0.9290 0.6940 0.1250], ...
    'DisplayName', '$|\Delta F_x/F_x|$');

% Highlight points below the 2% convergence threshold
converged_idx = abs_changes < 2.0;

h4 = plot(ax2, PPC_values(converged_idx), abs_changes(converged_idx), 's', ...
    'MarkerSize', 9, ...
    'MarkerFaceColor', [0.4660 0.6740 0.1880], ...
    'MarkerEdgeColor', [0.4660 0.6740 0.1880], ...
    'LineStyle', 'none', ...
    'DisplayName', 'Below 2\% threshold');

h5 = yline(ax2, 2.0, 'k:', ...
    'LineWidth', 2, ...
    'DisplayName', '2\% threshold');

h6 = xline(ax2, 20, 'k--', ...
    'LineWidth', 1.5, ...
    'DisplayName', '20 PPC reference');

ylabel(ax2, '$|\Delta F_x/F_x|$ [\%]', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

xlabel(ax2, 'Particles per cell (PPC)', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

legend(ax2, [h3, h4, h5, h6], ...
    'Location', 'best', ...
    'Interpreter', 'latex', ...
    'FontSize', 9);

ax2.FontSize = 11;
ax2.XTick = PPC_values;

valid_changes = abs_changes(~isnan(abs_changes));
if ~isempty(valid_changes)
    ylim(ax2, [0, max(2.5, max(valid_changes) * 1.2)]);
else
    ylim(ax2, [0, 2.5]);
end

%% Align x-axes exactly
linkaxes([ax1, ax2], 'x');
xlim(ax1, [min(PPC_values), max(PPC_values)]);

%% ------------------------------------------------------------------------
%  5. Export figure
%  ------------------------------------------------------------------------

exportgraphics(fig1, '../post_processing/fnum_convergence_2yaxis.png', ...
    'Resolution', 300);

%% ------------------------------------------------------------------------
%  6. Console summary
%  ------------------------------------------------------------------------

fprintf('\n==========================================================\n');
fprintf(' PARTICLE / FNUM CONVERGENCE ANALYSIS SUMMARY\n');
fprintf('==========================================================\n\n');

fprintf('Input file: %s\n', dir_info(idx(1)).name);
fprintf('Number of refinement levels: %d\n', n_cases);
fprintf('Grid resolution: N = %d, total cells = %d\n\n', grid_res, grid_cells);

fprintf('Drag force range:\n');
fprintf('  Minimum:      %.4e N\n', min(drag_forces_N));
fprintf('  Final:        %.4e N at PPC = %.0f\n', drag_final, PPC_values(end));
fprintf('  Converged at: %.4e N at PPC = %.0f\n', drag_converged, PPC_converged);
fprintf('  Spread:       %.4f %%\n\n', ...
    ((max(drag_forces_N) - min(drag_forces_N)) / drag_final) * 100);

fprintf('Step-by-step convergence:\n');

for i = 2:n_cases
    fprintf('  PPC %.0f -> %.0f: %.4f %% change\n', ...
    PPC_values(i-1), PPC_values(i), percent_change(i));
end

fprintf('\nCONVERGENCE STATUS:\n');

last_change = abs(percent_change(end));

if ~isnan(last_change)
    fprintf('  Last step: PPC %.0f -> %.0f gives %.4f %% change\n', ...
    PPC_values(end-1), PPC_values(end), last_change);

    if last_change < 2.0
        fprintf('  CONVERGED: final step is below the 2%% threshold.\n');
    else
        fprintf('  NOT CONVERGED: final step is above the 2%% threshold.\n');
    end
end

fprintf('\n==========================================================\n\n');