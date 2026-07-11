%% ========================================================================
%  DSMC Validation against Epstein, Stokes, Cunningham, Baines
%  Three figures:
%    1. Epstein-collapse (F/F_Epstein vs speed ratio S)
%    2. Continuum-model failure (Cunningham/DSMC, Stokes/DSMC vs d)
%    3. Scaling exponents (d^n and v^m)
%  ========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%  1. Load latest matrix results
%  ------------------------------------------------------------------------

output_dir = '../production_output/';
files = dir(fullfile(output_dir, 'matrix_results_*.csv'));

if isempty(files)
    error('No matrix_results_*.csv file found in ../production_output/');
end

[~, idx] = sort([files.datenum], 'descend');
target_file = fullfile(output_dir, files(idx(1)).name);

fprintf('Loading: %s\n', target_file);

data = readtable(target_file, 'VariableNamingRule', 'preserve');

diam_col = find_column(data, {'Diameter_m', 'diameter_m', 'Diameter'});
vel_col  = find_column(data, {'Velocity_ms', 'velocity_ms', 'Velocity'});
drag_col = find_column(data, {'Drag_Fx', 'Drag_Fx_N', 'DragForce', 'Fx'});

d = data.(diam_col);
v = data.(vel_col);
F_raw = abs(data.(drag_col));

valid = d > 0 & v > 0 & F_raw > 0;
d = d(valid);
v = v(valid);
F_raw = F_raw(valid);

diameters = unique(d);
velocities = unique(v);

%% ------------------------------------------------------------------------
%  2. Background force correction (containers.Map — same as other scripts)
%  ------------------------------------------------------------------------

F0_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
fprintf('\n--- Background force (F0) correction ---\n');
for i = 1:length(diameters)
    di = diameters(i);
    idx = d == di;
    v_i = v(idx);
    F_i = F_raw(idx);
    p = polyfit(v_i, F_i, 1);
    F0_map(di) = p(2);
    fprintf('  d=%s: F0 = %.4e N\n', diameter_label_plain(di), F0_map(di));
end
fprintf('---\n\n');

F_dsmc = zeros(size(F_raw));
for i = 1:length(F_raw)
    F_dsmc(i) = F_raw(i) - F0_map(d(i));
end
F_dsmc = max(F_dsmc, 0);

%% ------------------------------------------------------------------------
%  3. Gas properties and analytical models
%  ------------------------------------------------------------------------

T       = 293.15;              % K
T_wall  = T;                   % isothermal
P       = 100;                 % Pa
R_gas   = 8.314462618;         % J/(mol K)
M_H2    = 2.016e-3;            % kg/mol
mu_H2   = 8.8e-6;              % Pa s
k_B     = 1.380649e-23;        % J/K
m_H2    = 3.34e-27;            % kg
d_m_H2  = 2.72e-10;            % molecular diameter of H2 [m]
n_H2    = 2.47e22;             % number density [1/m^3]
lambda  = 1 / (sqrt(2) * pi * d_m_H2^2 * n_H2);
alpha   = 0.3;                 % accommodation coefficient

rho_gas = P * M_H2 / (R_gas * T);
v_mean  = sqrt(8 * k_B * T / (pi * m_H2));
v_beta  = sqrt(2 * k_B * T / m_H2);

r  = d ./ 2;
Kn = lambda ./ d;
S  = v ./ v_beta;

% Davies (1945) slip correction constants
A_const = 1.257;
B_const = 0.400;
C_const = 1.100;

% Epstein with accommodation correction
delta_eff = 1 + alpha * pi / 8;
F_epstein = (4/3) * pi .* r.^2 .* rho_gas .* v_mean .* v .* delta_eff;

% Stokes continuum
F_stokes = 6 * pi * mu_H2 .* r .* v;

% Cunningham slip correction (Davies 1945) — uses Knudsen based on radius
Kn_slip = 2 * lambda ./ d;     % radius-based: Kn_slip = 2*lambda/d
Cc = 1 + Kn_slip .* (A_const + B_const .* exp(-C_const ./ Kn_slip));
F_cunningham = F_stokes ./ Cc;

%% ------------------------------------------------------------------------
%  4. Baines generalized free-molecular model
%  ------------------------------------------------------------------------

v_theory = logspace(log10(5), log10(800), 500);
S_theory = v_theory ./ v_beta;

base  = (S_theory + 1./(2*S_theory)).*exp(-S_theory.^2)./sqrt(pi) ...
      + (S_theory.^2 + 1 - 1./(4*S_theory.^2)).*erf(S_theory);

% Diffuse reflection part (accommodated)
R_d = 2*pi * (diameters(:)/2).^2 * P ...
    .* ( base + (S_theory/3).*sqrt(pi).*sqrt(T_wall/T) );

% Specular reflection part
R_s = 2*pi * (diameters(:)/2).^2 * P .* base;

% Combined Baines force
F_baines = alpha .* R_d + (1 - alpha) .* R_s;
% F_baines: size = [n_diameters, n_v_theory]

%% ------------------------------------------------------------------------
%  5. Ratios
%  ------------------------------------------------------------------------

R_epstein    = F_dsmc ./ F_epstein;
R_stokes     = F_stokes ./ F_dsmc;
R_cunningham = F_cunningham ./ F_dsmc;

% Baines/Epstein ratio for the theoretical curve (per diameter)
S_data = v ./ v_beta;
F_baines_data = zeros(size(v));
for i = 1:length(diameters)
    di = diameters(i);
    idx = d == di;
    for j = find(idx)'
        [~, tidx] = min(abs(S_theory - S_data(j)));
        F_baines_data(j) = F_baines(i, tidx);
    end
end

%% ========================================================================
%  FIGURE 1: DSMC collapse against Epstein scaling (+ Baines reference)
%  ========================================================================

fig1 = figure('Color','w','Position',[100,100,850,560]);
hold on; grid on; box on;

colors = lines(length(diameters));

% Factor-of-two band
x_band = [min(S)*0.8, max(S)*1.2, max(S)*1.2, min(S)*0.8];
y_band = [0.5, 0.5, 2.0, 2.0];
patch(x_band, y_band, [0.85 0.92 1.00], ...
    'EdgeColor','none', 'FaceAlpha',0.35, 'HandleVisibility','off');

% Baines/F_Epstein reference curves per diameter (dashed)
for i = 1:length(diameters)
    di = diameters(i);
    ratio_b = F_baines(i,:) ./ ...
        ((4/3)*pi*(di/2)^2*rho_gas*v_mean.*v_theory*delta_eff);
    plot(S_theory, ratio_b, '--', 'LineWidth', 1.3, ...
        'Color', colors(i,:), 'HandleVisibility','off');
end

% DSMC data
for i = 1:length(diameters)
    idx_d = d == diameters(i);
    S_i = S(idx_d);
    R_i = R_epstein(idx_d);
    [S_i, sort_idx] = sort(S_i);
    R_i = R_i(sort_idx);
    plot(S_i, R_i, 'o-', 'LineWidth', 2, 'MarkerSize', 7, ...
        'Color', colors(i,:), 'MarkerFaceColor', 'w', ...
        'DisplayName', diameter_label(diameters(i)));
end

% Perfect-agreement line
plot(S_theory, ones(size(S_theory)), 'k--', 'LineWidth', 1.8, ...
    'DisplayName', 'Epstein ($F/F_E = 1$)');

% Baines legend entry (single style representative)
h_baines = plot(nan, nan, '--k', 'LineWidth', 1.3, ...
    'DisplayName', 'Baines (per $d$)');

set(gca, 'XScale','log', 'YScale','log', 'FontSize',12);
xlabel('Molecular speed ratio $S = v_0/v_\beta$', ...
    'Interpreter','latex', 'FontSize',14);
ylabel('$F_{\mathrm{DSMC}}/F_{\mathrm{Epstein}}$  (dashed = Baines)', ...
    'Interpreter','latex', 'FontSize',14);
title('DSMC drag vs Epstein scaling with Baines free-molecular reference', ...
    'Interpreter','latex', 'FontSize',14);
legend('Location','best', 'Interpreter','latex', 'FontSize',10);

% Annotate Kn for each diameter in the margin
yl = ylim;
for i = 1:length(diameters)
    Kn_i = lambda / diameters(i);
    s_i = min(S(idx_d)) * 0.85;
    text(s_i, yl(2)*0.92 - 0.06*(i-1)*yl(2), ...
        sprintf('$d$=%s, $Kn$=$%.0f$', ...
            strrep(diameter_label_plain(diameters(i)), ' ', '~'), Kn_i), ...
        'FontSize', 8, 'Color', colors(i,:), 'Interpreter','latex');
end

exportgraphics(fig1, 'validation_epstein_collapse.png', 'Resolution', 300);
savefig(fig1, 'validation_epstein_collapse.fig');

%% ========================================================================
%  FIGURE 2: Continuum model failure (+ Baines/DSMC panel)
%  ========================================================================

fig2 = figure('Color','w','Position',[100,100,1100,720]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% --- Panel (a): Cunningham / DSMC ---
ax1 = nexttile;
hold on; grid on; box on;
for i = 1:length(velocities)
    idx_v = v == velocities(i);
    d_i = d(idx_v); R_i = R_cunningham(idx_v);
    [d_i, sort_idx] = sort(d_i); R_i = R_i(sort_idx);
    loglog(d_i*1e6, R_i, 'o-', 'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', sprintf('$v_0=%g$ m/s', velocities(i)));
end
yline(1, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('Droplet diameter $d$ ($\mu$m)', 'Interpreter','latex', 'FontSize',13);
ylabel('$F_{\mathrm{Cunningham}}/F_{\mathrm{DSMC}}$', 'Interpreter','latex', 'FontSize',13);
title('(a) Slip-corrected continuum vs DSMC', 'Interpreter','latex', 'FontSize',13);
legend('Location','best', 'Interpreter','latex', 'FontSize',7);
ax1.FontSize = 11;

% --- Panel (b): Stokes / DSMC ---
ax2 = nexttile;
hold on; grid on; box on;
for i = 1:length(velocities)
    idx_v = v == velocities(i);
    d_i = d(idx_v); R_i = R_stokes(idx_v);
    [d_i, sort_idx] = sort(d_i); R_i = R_i(sort_idx);
    loglog(d_i*1e6, R_i, 'o-', 'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', sprintf('$v_0=%g$ m/s', velocities(i)));
end
yline(1, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('Droplet diameter $d$ ($\mu$m)', 'Interpreter','latex', 'FontSize',13);
ylabel('$F_{\mathrm{Stokes}}/F_{\mathrm{DSMC}}$', 'Interpreter','latex', 'FontSize',13);
title('(b) No-slip continuum vs DSMC', 'Interpreter','latex', 'FontSize',13);
legend('Location','best', 'Interpreter','latex', 'FontSize',7);
ax2.FontSize = 11;

% --- Panel (c): Baines / DSMC ---
ax3 = nexttile;
hold on; grid on; box on;
R_baines = F_baines_data ./ F_dsmc;
for i = 1:length(velocities)
    idx_v = v == velocities(i);
    d_i = d(idx_v); R_i = R_baines(idx_v);
    [d_i, sort_idx] = sort(d_i); R_i = R_i(sort_idx);
    loglog(d_i*1e6, R_i, 'o-', 'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', sprintf('$v_0=%g$ m/s', velocities(i)));
end
yline(1, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('Droplet diameter $d$ ($\mu$m)', 'Interpreter','latex', 'FontSize',13);
ylabel('$F_{\mathrm{Baines}}/F_{\mathrm{DSMC}}$', 'Interpreter','latex', 'FontSize',13);
title('(c) Baines free-molecular vs DSMC', 'Interpreter','latex', 'FontSize',13);
legend('Location','best', 'Interpreter','latex', 'FontSize',7);
ax3.FontSize = 11;

% --- Panel (d): F/v linear fit slope comparison ---
ax4 = nexttile;
hold on; grid on; box on;
v_plot = logspace(log10(10), log10(800), 200);
for i = 1:length(diameters)
    di = diameters(i);
    ri = di / 2;
    Kn_i = lambda / di;
    % Baines F/v at low v → asymptotes to Epstein
    S_p = v_plot / v_beta;
    base_p = (S_p + 1./(2*S_p)).*exp(-S_p.^2)./sqrt(pi) ...
           + (S_p.^2 + 1 - 1./(4*S_p.^2)).*erf(S_p);
    Rd_p = 2*pi*ri^2*P.*(base_p + (S_p/3).*sqrt(pi).*sqrt(T_wall/T));
    Rs_p = 2*pi*ri^2*P.*base_p;
    Fb_p = alpha.*Rd_p + (1-alpha).*Rs_p;
    % DSMC data points
    idx_d = d == di;
    v_di = v(idx_d); F_di = F_dsmc(idx_d);
    [v_di, si] = sort(v_di); F_di = F_di(si);
    loglog(v_di, F_di, 'o-', 'LineWidth', 2, 'MarkerSize', 7, ...
        'Color', colors(i,:), 'MarkerFaceColor', 'w', ...
        'DisplayName', [diameter_label_plain(di) ' DSMC']);
    loglog(v_plot, Fb_p, '--', 'LineWidth', 1.3, ...
        'Color', colors(i,:), 'HandleVisibility','off');
end
xlabel('Velocity $v_0$ (m/s)', 'Interpreter','latex', 'FontSize',13);
ylabel('$F_x$ (N)', 'Interpreter','latex', 'FontSize',13);
title('(d) DSMC vs Baines (dashed) per diameter', 'Interpreter','latex', 'FontSize',13);
legend('Location','best', 'Interpreter','latex', 'FontSize',7);
ax4.FontSize = 11;

sgtitle('Continuum and kinetic model comparison with DSMC', ...
    'Interpreter','latex', 'FontSize',15, 'FontWeight','bold');

exportgraphics(fig2, 'validation_continuum_failure.png', 'Resolution', 300);
savefig(fig2, 'validation_continuum_failure.fig');

%% ========================================================================
%  FIGURE 3: Scaling exponent validation
%  ========================================================================

fig3 = figure('Color','w','Position',[100,100,1050,480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% --- Diameter exponent: F ~ d^n at fixed velocity ---
n_exp = zeros(length(velocities),1);
for i = 1:length(velocities)
    idx_v = v == velocities(i);
    p = polyfit(log10(d(idx_v)), log10(F_dsmc(idx_v)), 1);
    n_exp(i) = p(1);
end

ax1 = nexttile;
hold on; grid on; box on;
plot(velocities, n_exp, 'o-', 'LineWidth', 2.3, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
yline(2, 'k--', 'LineWidth', 1.6, 'DisplayName', '$d^2$ projected-area scaling');
xlabel('Flow velocity $v_0$ (m/s)', 'Interpreter','latex', 'FontSize',13);
ylabel('Fitted diameter exponent $n$ in $F_x \propto d^n$', 'Interpreter','latex', 'FontSize',13);
title('(a) Diameter scaling exponent', 'Interpreter','latex', 'FontSize',14);
legend('Location','best', 'Interpreter','latex', 'FontSize',9);
ax1.FontSize = 11;

% --- Velocity exponent: F ~ v^m at fixed diameter ---
m_exp = zeros(length(diameters),1);
for i = 1:length(diameters)
    idx_d = d == diameters(i);
    p = polyfit(log10(v(idx_d)), log10(F_dsmc(idx_d)), 1);
    m_exp(i) = p(1);
end

ax2 = nexttile;
hold on; grid on; box on;
semilogx(diameters*1e6, m_exp, 's-', 'LineWidth', 2.3, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
yline(1, 'k--', 'LineWidth', 1.6, 'DisplayName', 'Linear velocity scaling');
xlabel('Droplet diameter $d$ ($\mu$m)', 'Interpreter','latex', 'FontSize',13);
ylabel('Fitted velocity exponent $m$ in $F_x \propto v_0^m$', 'Interpreter','latex', 'FontSize',13);
title('(b) Velocity scaling exponent', 'Interpreter','latex', 'FontSize',14);
legend('Location','best', 'Interpreter','latex', 'FontSize',9);
ax2.FontSize = 11;

sgtitle('DSMC scaling checks for rarefied molecular drag', ...
    'Interpreter','latex', 'FontSize',15, 'FontWeight','bold');

exportgraphics(fig3, 'validation_scaling_exponents.png', 'Resolution', 300);
savefig(fig3, 'validation_scaling_exponents.fig');

%% ========================================================================
%  Console summary
%  ========================================================================

fprintf('\n=====================================================\n');
fprintf(' VALIDATION SUMMARY (with Baines)\n');
fprintf('=====================================================\n');

fprintf('\nEpstein-normalized DSMC ratio:\n');
fprintf('  min(F_DSMC/F_Epstein) = %.3f\n', min(R_epstein));
fprintf('  max(F_DSMC/F_Epstein) = %.3f\n', max(R_epstein));
fprintf('  mean                  = %.3f\n', mean(R_epstein));

fprintf('\nBaines/DSMC ratio:\n');
fprintf('  min(F_Baines/F_DSMC)  = %.3f\n', min(R_baines));
fprintf('  max(F_Baines/F_DSMC)  = %.3f\n', max(R_baines));
fprintf('  mean                  = %.3f\n', mean(R_baines));

fprintf('\nContinuum overprediction:\n');
fprintf('  Cunningham/DSMC range = %.2e to %.2e\n', min(R_cunningham), max(R_cunningham));
fprintf('  Stokes/DSMC range     = %.2e to %.2e\n', min(R_stokes), max(R_stokes));

fprintf('\nMean absolute deviation from analytical models:\n');
fprintf('  |DSMC - Epstein|/Epstein     = %.3f\n', mean(abs(R_epstein - 1)));
fprintf('  |DSMC - Baines|/Baines       = %.3f\n', mean(abs(R_baines - 1)));

fprintf('\nDiameter scaling exponents:\n');
for i = 1:length(velocities)
    fprintf('  v = %g m/s: n = %.3f\n', velocities(i), n_exp(i));
end

fprintf('\nVelocity scaling exponents:\n');
for i = 1:length(diameters)
    fprintf('  d = %s: m = %.3f\n', diameter_label_plain(diameters(i)), m_exp(i));
end

fprintf('\nSaved figures:\n');
fprintf('  validation_epstein_collapse.png\n');
fprintf('  validation_continuum_failure.png\n');
fprintf('  validation_scaling_exponents.png\n');
fprintf('\n=====================================================\n');

%% ========================================================================
%  Local helper functions
%  ========================================================================

function col = find_column(tbl, candidates)
    names = tbl.Properties.VariableNames;
    for i = 1:length(candidates)
        match = strcmp(names, candidates{i});
        if any(match)
            idx = find(match, 1);
            col = names{idx};
            return;
        end
    end
    error('None of the expected columns were found: %s', strjoin(candidates, ', '));
end

function label = diameter_label(d)
    if d < 1e-6
        label = sprintf('$%.0f$ nm', d*1e9);
    else
        label = sprintf('$%.0f~\\mu\\mathrm{m}$', d*1e6);
    end
end

function label = diameter_label_plain(d)
    if d < 1e-6
        label = sprintf('%.0f nm', d*1e9);
    else
        label = sprintf('%.0f um', d*1e6);
    end
end
