%% ========================================================================
%  Mesh Convergence Post-Processing
%  Reads the CSV from run_mesh_refinement.sh and produces a publication-ready
%  plot of drag force vs. number of surface triangles.
%  ========================================================================
clear; clc; close all;

%% 1. Find the newest mesh_convergence CSV
csv_pattern = '../runtime_output/mesh_convergence_*.csv';
files = dir(csv_pattern);

if isempty(files)
    error('No mesh_convergence_*.csv found in ../runtime_output/');
end

% Sort by datenum (newest last)
[~, idx] = sort([files.datenum]);
csv_file = fullfile(files(idx(end)).folder, files(idx(end)).name);

fprintf('Reading: %s\n', csv_file);

%% 2. Load data
opts = detectImportOptions(csv_file, 'NumHeaderLines', 1);
data = readtable(csv_file, opts);

N_tri = data.Ntriangles;
Drag_N = data.Drag_Fx_N;
Drag_pN = data.Drag_Fx_pN;
PPC = data.PPC_actual;
Mesh_File = data.Mesh_File;

%% 3. Compute statistics
drag_final = Drag_pN(end);
n_levels = length(N_tri);

% Deviation from final (finest mesh) value
pct_dev = ((Drag_pN - drag_final) / drag_final) * 100;

% Stepwise change
step_change = nan(1, n_levels);
for i = 2:n_levels
    step_change(i) = ((Drag_pN(i) - Drag_pN(i-1)) / Drag_pN(i-1)) * 100;
end

fprintf('\n==========================================================\n');
fprintf(' MESH CONVERGENCE ANALYSIS SUMMARY\n');
fprintf('==========================================================\n\n');

for i = 1:n_levels
    % Extract just the filename from the path
    [~, fname, fext] = fileparts(char(Mesh_File{i}));
    short_name = [fname, fext];
    
    fprintf('  %6d triangles (%s):\n', N_tri(i), short_name);
    fprintf('    Drag = %.4f pN  (%.6e N)\n', Drag_pN(i), Drag_N(i));
    fprintf('    Deviation from finest: %+.3f%%\n', pct_dev(i));
    if ~isnan(step_change(i))
        fprintf('    Change from previous:  %+.3f%%\n', step_change(i));
    end
    fprintf('    PPC = %.1f\n', PPC(i));
    fprintf('\n');
end

% Convergence assessment
fprintf('CONVERGENCE ASSESSMENT:\n');
if n_levels >= 2
    last_change = step_change(end);
    if abs(last_change) < 2.0
        fprintf('  ✓ CONVERGED: %d→%d triangles change = %+.3f%% (< 2%%)\n', ...
            N_tri(end-1), N_tri(end), last_change);
    else
        fprintf('  ⚠ NOT CONVERGED: last change = %+.3f%% (≥ 2%%)\n', last_change);
        fprintf('    Consider a finer mesh beyond %d triangles.\n', N_tri(end));
    end
end
fprintf('==========================================================\n\n');

%% 4. Plot
fig = figure('Name', 'Mesh Convergence', ...
    'NumberTitle', 'off', 'Position', [100, 100, 900, 500]);

% Drag force vs number of triangles
subplot(1, 2, 1);
hold on; grid on;

plot(N_tri, Drag_pN, 'o-', 'LineWidth', 2.5, 'MarkerSize', 12, ...
    'MarkerFaceColor', [0.2 0.4 0.8], 'Color', [0.2 0.4 0.8]);

yline(drag_final, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Final: %.2f pN', drag_final));

xlabel('Number of surface triangles', 'FontSize', 12);
ylabel('Drag force |F_x| (pN)', 'FontSize', 12);
title('(a) Drag vs. mesh resolution', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'XScale', 'log', 'FontSize', 10);
xticks(N_tri);
xlim([min(N_tri)*0.8, max(N_tri)*1.2]);

% Add triangle count labels above markers
for i = 1:n_levels
    text(N_tri(i), Drag_pN(i)*1.05, sprintf('%d', N_tri(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

% Percent deviation from finest mesh
subplot(1, 2, 2);
hold on; grid on;

colors = repmat([0.2 0.8 0.3], n_levels, 1);
for i = 1:n_levels
    if abs(pct_dev(i)) > 2.0
        colors(i, :) = [0.9 0.2 0.2];
    end
end

bar(1:n_levels, pct_dev, 'FaceColor', 'flat', 'CData', colors, ...
    'EdgeColor', 'black', 'LineWidth', 1);

yline(0, 'k-', 'LineWidth', 1);
yline(2, 'r--', 'LineWidth', 1.5, 'DisplayName', '±2% threshold');
yline(-2, 'r--', 'LineWidth', 1.5);

xlabel('Mesh refinement level', 'FontSize', 12);
ylabel('% deviation from finest mesh', 'FontSize', 12);
title('(b) Relative error vs. finest mesh', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);

xticklabels(cellstr(num2str(N_tri)));
xlim([0.5, n_levels + 0.5]);

sgtitle(sprintf('Surface Mesh Convergence — Worst Case (d=10\\mum, v=700 m/s)'), ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save
savefig(fig, '../runtime_output/mesh_convergence_plots.fig');
exportgraphics(fig, '../runtime_output/mesh_convergence_plots.png', 'Resolution', 150);
fprintf('Saved: ../runtime_output/mesh_convergence_plots.png\n');
