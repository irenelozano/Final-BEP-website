%% ========================================================================
%  steadystate_proof_plot.m
%
%  PURPOSE:
%    Generates the steady-state proof figure for the BEP report.
%    Reads the output of run_steadystate_proof.sh and produces a
%    two-panel publication-quality figure showing:
%      (a) Full force time history (transient + steady), with the
%          transient region shaded and the steady window highlighted.
%      (b) Cumulative running mean in the steady window only,
%          with ±2% convergence band.
%
%  USAGE:
%    steadystate_proof_plot()                    % auto-detects newest run
%    steadystate_proof_plot('../runtime_output/steadystate_proof')
%    steadystate_proof_plot(timeseries_file, config_file)  % explicit files
%
%  WHAT THIS PROVES:
%    1. The force time series shows no trend after the transient —
%       the flow field has fully developed.
%    2. The cumulative running mean converges to within ±2% of the
%       final value within the steady window — the time average is stable.
%    3. The ±1σ noise band around the mean is small relative to the
%       mean — signal-to-noise is acceptable for DSMC.
%
%  OUTPUT:
%    steadystate_proof_<label>.png  (in the same folder as the data)
%    steadystate_proof_<label>.fig
% =========================================================================

function steadystate_proof_plot(varargin)

% -------------------------------------------------------------------------
% 1. FILE RESOLUTION
% -------------------------------------------------------------------------
if nargin == 0 
    % Auto-detect newest timeseries file
    search_dir = '../runtime_output/steadystate_proof';
    files = dir(fullfile(search_dir, '*_timeseries.txt'));
    if isempty(files)
        error('No *_timeseries.txt files found in: %s', search_dir);
    end
    [~, idx] = sort([files.datenum], 'descend');
    ts_file = fullfile(files(idx(1)).folder, files(idx(1)).name);
else
    ts_file = varargin{1};
end

fprintf('Processing timeseries: %s\n', ts_file);

% -------------------------------------------------------------------------
% 2. HARDCODED PARAMETERS (For d=10um, v=700m/s case)
% -------------------------------------------------------------------------
dt           = 1.72748e-08; 
trans_steps  = 25000;
steady_steps = 50000;
record_every = 200;

% -------------------------------------------------------------------------
% 3. READ TIMESERIES DATA (SPARTA surf dump format)
%    Format per block (after ITEM: TIMESTEP):
%      ITEM: NUMBER OF SURFS / 192 / ITEM: BOX BOUNDS / 3 box lines / ITEM: SURFS header
%      then 192 data rows: id Fx Fy Fz
%    Strategy: scan every line; accept it as a data row only when all 4
%    fields parse as finite numbers. This is format-agnostic.
% -------------------------------------------------------------------------
raw_text = fileread(ts_file);
blocks   = regexp(raw_text, 'ITEM: TIMESTEP', 'split');
blocks   = blocks(2:end);   % discard leading empty string before first block

if isempty(blocks)
    error('No ITEM: TIMESTEP blocks found in: %s', ts_file);
end

n_blocks  = length(blocks);
steps     = nan(n_blocks, 1);
Fx_nN     = nan(n_blocks, 1);

for k = 1:n_blocks
    lines = splitlines(strtrim(blocks{k}));

    % First non-empty line after the split is always the timestep number
    for li = 1:length(lines)
        val = str2double(strtrim(lines{li}));
        if ~isnan(val)
            steps(k) = val;
            break;
        end
    end

    % Accumulate Fx (col 2) over all data rows in this block.
    % Data rows have exactly 4 finite numeric fields: id Fx Fy Fz.
    % Skip all ITEM: headers and box-bound lines (they have 1-3 or non-numeric fields).
    fx_sum = 0;
    n_data = 0;
    for j = 1:length(lines)
        ln = strtrim(lines{j});
        if isempty(ln) || startsWith(ln, 'ITEM:'), continue; end
        parts = strsplit(ln);
        if length(parts) ~= 4, continue; end
        vals = str2double(parts);
        if any(isnan(vals)), continue; end   % skip header/bound lines
        fx_sum = fx_sum + vals(2);           % col 2 = Fx (signed)
        n_data = n_data + 1;
    end

    if n_data > 0
        Fx_nN(k) = abs(fx_sum) * 1e9;   % N -> nN (drag magnitude)
    else
        Fx_nN(k) = 0;   % no surface hits in this block -> zero drag
    end
end

% Drop any blocks where the timestep couldn't be parsed
valid = ~isnan(steps);
steps  = steps(valid);
Fx_nN  = Fx_nN(valid);

% Physical time at each recorded step
phys_time_us = steps * dt * 1e6;   % microseconds

% Split into transient and steady windows
trans_mask  = steps <= trans_steps;
steady_mask = steps >  trans_steps;

Fx_trans    = Fx_nN(trans_mask);
t_trans     = phys_time_us(trans_mask);
Fx_steady   = Fx_nN(steady_mask);
t_steady    = phys_time_us(steady_mask);
steps_steady = steps(steady_mask);

% -------------------------------------------------------------------------
% 4. STATISTICS (steady window only)
% -------------------------------------------------------------------------
mean_Fx   = mean(Fx_steady);
std_Fx    = std(Fx_steady);
cv_pct    = 0;
if mean_Fx > 0
    cv_pct = (std_Fx / mean_Fx) * 100;
end

% Cumulative running mean over ALL data from step 0
cum_mean_full = cumsum(Fx_nN) ./ (1:length(Fx_nN))';

% Cumulative running mean in steady window only (for panel b convergence check)
cum_mean  = cumsum(Fx_steady) ./ (1:length(Fx_steady))';

% Relative deviation of running mean from final value (%)
rel_dev   = zeros(size(cum_mean));
if mean_Fx > 0
    rel_dev = (cum_mean - mean_Fx) / mean_Fx * 100;
end

% When does cumulative mean enter and STAY within ±2% band?
band_pct  = 2.0;
outside   = abs(rel_dev) >= band_pct;
last_out  = find(outside, 1, 'last');

if isempty(last_out)
    conv_idx = 1; % It was always inside
elseif last_out < length(rel_dev)
    conv_idx = last_out + 1; % The point right after the last time it was outside
else
    conv_idx = []; % It never converged
end
% -------------------------------------------------------------------------
% 5. FIGURE LAYOUT
% -------------------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [80 80 950 700]);

% Color palette (matches your existing report style)
C_transient = [0.85 0.85 0.85];   % light grey shading for transient region
C_steady    = [0.90 0.96 1.00];   % light blue shading for steady region
C_data      = [0.00 0.18 0.35];   % dark navy — block data dots
C_cumul     = [0.85 0.33 0.10];   % orange-red — running mean
C_band      = [0.30 0.75 0.93];   % sky blue — ±1σ band
C_converge  = [0.47 0.67 0.19];   % green — convergence marker

% =========================================================================
% PANEL (a): Full time history — transient + steady
% =========================================================================
ax1 = subplot(2, 1, 1);
hold on; box on; grid on;

% Shaded transient region
t_end_trans = t_trans(end);
t_start_all = phys_time_us(1);
t_end_all   = phys_time_us(end);
y_lo = min(Fx_nN) * 0.90;
y_hi = max(Fx_nN) * 1.10;
if y_lo == y_hi || isnan(y_lo) || isnan(y_hi)
    y_lo = 0; y_hi = 1;
end

fill([t_start_all t_end_trans t_end_trans t_start_all], ...
     [y_lo y_lo y_hi y_hi], C_transient, ...
     'EdgeColor', 'none', 'FaceAlpha', 0.7, 'HandleVisibility', 'off');

% Shaded steady region
fill([t_end_trans t_end_all t_end_all t_end_trans], ...
     [y_lo y_lo y_hi y_hi], C_steady, ...
     'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

% ±1σ band (steady window only) — drawn as horizontal band
fill([t_steady(1) t_steady(end) t_steady(end) t_steady(1)], ...
     [mean_Fx - std_Fx, mean_Fx - std_Fx, mean_Fx + std_Fx, mean_Fx + std_Fx], ...
     C_band, 'EdgeColor', 'none', 'FaceAlpha', 0.35, 'DisplayName', '\pm1\sigma noise band');

% All data points (transient + steady) as scatter
scatter(t_trans,  Fx_trans,  28, [0.6 0.6 0.6], 'filled', 'DisplayName', 'Transient blocks');
scatter(t_steady, Fx_steady, 35, C_data, 'filled', 'DisplayName', 'Steady-state blocks');

% Running cumulative mean from t=0 over ALL data
plot(phys_time_us, cum_mean_full, '-', 'Color', C_cumul, 'LineWidth', 2, ...
     'DisplayName', 'Cumulative running mean (from t=0)');

% Final mean line
yline(mean_Fx, 'k--', 'LineWidth', 1.8, ...
      'DisplayName', sprintf('Mean: %.3f nN', mean_Fx));

% Transient boundary marker
xline(t_end_trans, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, ...
     'HandleVisibility', 'off');

ylim([y_lo y_hi]);
xlim([t_start_all t_end_all]);
xlabel('Physical time (\mus)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Drag force |F_x| (nN)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('(a) Full force time history — transient + steady-state sampling', ''), ...
      'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);

% Annotation box: noise metric
annotation_str = sprintf('CV = %.1f%%\nMean = %.3f nN\n\\sigma = %.3f nN', ...
    cv_pct, mean_Fx, std_Fx);
text(0.98, 0.05, annotation_str, 'Units', 'normalized', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
     'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], ...
     'Margin', 4);

% =========================================================================
% PANEL (b): Cumulative running mean convergence (steady window only)
% =========================================================================
ax2 = subplot(2, 1, 2);
hold on; box on; grid on;

% ±2% convergence band (shaded)
fill([t_steady(1) t_steady(end) t_steady(end) t_steady(1)], ...
     [-band_pct -band_pct band_pct band_pct], ...
     [0.85 0.95 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
     'DisplayName', sprintf('\\pm%.0f%% convergence band', band_pct));

% Zero line (final mean reference)
yline(0, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Final mean (reference)');

% ±2% band borders
yline( band_pct, '--', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(-band_pct, '--', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Relative deviation of cumulative mean
plot(t_steady, rel_dev, 's-', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.8, ...
     'MarkerSize', 5, 'MarkerFaceColor', [0.75 0.58 0.83], ...
     'DisplayName', 'Cumulative mean deviation');

% Mark convergence point
if ~isempty(conv_idx)
    scatter(t_steady(conv_idx), rel_dev(conv_idx), 120, C_converge, ...
            'pentagram', 'filled', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Converged at t = %.1f \\mus', t_steady(conv_idx)));
end

xlim([t_steady(1) t_steady(end)]);
y_range = max(abs(rel_dev)) * 1.3;
ylim([-max(y_range, band_pct*1.5), max(y_range, band_pct*1.5)]);

xlabel('Physical time (\mus)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative deviation from final mean (%)', 'FontSize', 12, 'FontWeight', 'bold');
title('(b) Convergence trajectory — cumulative mean relative to final value', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);

% -------------------------------------------------------------------------
% 6. SUPERTITLE AND FORMATTING
% -------------------------------------------------------------------------
sgtitle('Steady-State Verification: d = 10 \mum,  v_0 = 700 m/s', ...
    'FontSize', 14, 'FontWeight', 'bold');

set(findall(fig, '-property', 'FontName'), 'FontName', 'Helvetica');

% -------------------------------------------------------------------------
% 7. SAVE
% -------------------------------------------------------------------------
[ts_dir, ts_name] = fileparts(ts_file);
out_base = fullfile(ts_dir, strrep(ts_name, '_timeseries', '_proof_figure'));

savefig(fig, [out_base '.fig']);
print(fig, [out_base '.png'], '-dpng', '-r200');

fprintf('  Figure saved: %s.png\n\n', out_base);



