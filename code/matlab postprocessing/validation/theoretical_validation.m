%% ========================================================================
%  DSMC Validation: Drag Force vs. Theoretical Models
%  Log-log comparison for selected velocity cases
%  ========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%  1. User settings
%  ------------------------------------------------------------------------

data_file = 'matrix_results_20260601_125031.csv';

% Velocity cases to plot [m/s]
velocity_cases = [20, 100, 550, 700];

% Droplet diameters used in the DSMC matrix
d_case_labels = [100e-9, 500e-9, 1e-6, 5e-6, 10e-6];
d_case_text   = {'100 nm', '500 nm', '1 \mum', '5 \mum', '10 \mum'};

%% ------------------------------------------------------------------------
%  2. Load DSMC data once
%  ------------------------------------------------------------------------

data = readtable(data_file);

%% ------------------------------------------------------------------------
%  Background force correction: fit and subtract F0 per diameter
% ------------------------------------------------------------------------
[unique_D, ~, ~] = unique(data.Diameter_m);
F0_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
fprintf('--- Background force (F0) correction ---\n');

for i = 1:length(unique_D)
    idx = data.Diameter_m == unique_D(i);
    v_i = data.Velocity_ms(idx);
    F_i = abs(data.Drag_Fx(idx));
    p = polyfit(v_i, F_i, 1);
    F0_map(unique_D(i)) = p(2);

    if unique_D(i) < 1e-6
        lbl = sprintf('%.0f nm', unique_D(i)*1e9);
    else
        lbl = sprintf('%.0f um', unique_D(i)*1e6);
    end

    fprintf('  d=%s: F0 = %.4e N\n', lbl, p(2));
end

fprintf('---\n\n');

%% ------------------------------------------------------------------------
%  3. Gas properties: H2 at 100 Pa and 293 K
%  ------------------------------------------------------------------------

T = 293.15;              % Gas temperature [K]
P = 100;                 % Gas pressure [Pa]
R = 8.314;               % Universal gas constant [J/(mol K)]
M_H2 = 2.016e-3;         % H2 molar mass [kg/mol]
mu_H2 = 8.8e-6;          % Dynamic viscosity of H2 [Pa s]

k_B = 1.38e-23;          % Boltzmann constant [J/K]
m_H2 = 3.34e-27;         % H2 molecular mass [kg]

d_m_H2 = 2.72e-10;       % Molecular diameter of H2 (from hydrogen.vhs) [m]
n_H2   = 2.47e22;        % Number density (from in.drag: global nrho 2.47e22) [1/m^3]
lambda = 1 / (sqrt(2) * pi * d_m_H2^2 * n_H2); % Mean free path [m]
alpha_CLL = 0.3;         % Accommodation coefficient / diffuse weighting

T_wall = T;              % Isothermal droplet wall assumption [K]

rho_gas = (P * M_H2) / (R * T);                 % Gas density [kg/m^3]
v_mean = sqrt((8 * k_B * T) / (pi * m_H2));     % Mean molecular speed [m/s]
v_beta = sqrt((2 * k_B * T) / m_H2);            % Most probable molecular speed [m/s]

%% ------------------------------------------------------------------------
%  4. Continuous diameter range for analytical curves
%  ------------------------------------------------------------------------

d_theory = logspace(log10(100e-9), log10(10e-6), 300);
r_theory = d_theory / 2;
Kn_regime = lambda ./ d_theory;     % diameter-based Kn — regime labels only
Kn_slip   = lambda ./ r_theory;     % radius-based Kn = 2*lambda./d — Cunningham

%% ========================================================================
%  5. Loop over velocity cases
%  ========================================================================

for v_target = velocity_cases

    %% --------------------------------------------------------------------
    %  5.1 Extract DSMC data for selected velocity
    %  --------------------------------------------------------------------

    idx = data.Velocity_ms == v_target;
    d_dsmc = data.Diameter_m(idx);
    F_dsmc = abs(data.Drag_Fx(idx));

    % Subtract per-diameter F0
    for j = 1:length(d_dsmc)
        F_dsmc(j) = F_dsmc(j) - F0_map(d_dsmc(j));
    end

    % Check if your CSV contains the standard deviation.
    % If not, it defaults to your 2% bound.
    if ismember('Drag_Fx_Std', data.Properties.VariableNames)
        F_std = data.Drag_Fx_Std(idx);
    else
        F_std = F_dsmc * 0.02; % 2% conservative noise floor bound
    end

    % Sort data by diameter
    [d_dsmc, sort_idx] = sort(d_dsmc);
    F_dsmc = F_dsmc(sort_idx);
    F_std = F_std(sort_idx);

    %% --------------------------------------------------------------------
    %  5.2 Analytical drag models
    %  --------------------------------------------------------------------

    % A. Continuum Stokes drag
    F_stokes = 6 * pi * mu_H2 .* r_theory .* v_target;

    % B. Cunningham slip correction (Davies 1945 coefficients)
    % Davies coefficients require Kn on the particle radius:
    % Kn = lambda/r = 2*lambda/d
    alpha_c = 1.257;
    beta_c  = 0.400;
    gamma_c = 1.100;

    C_factor = 1 + Kn_slip .* ...
        (alpha_c + beta_c .* exp(-gamma_c ./ Kn_slip));

    F_cunningham = F_stokes ./ C_factor;

    % C. Epstein low-speed free-molecular drag
    F_epstein = (4/3) * pi * rho_gas .* (r_theory.^2) .* v_mean .* v_target .* ...
        (1 + (pi/8) * alpha_CLL);

    % D. Baines generalized free-molecular resistance
    % Accounts for finite molecular speed ratio effects.
    S = v_target / v_beta;

    R_d = (2 * pi .* r_theory.^2 .* P) .* ...
        ( ...
        (S + 1./(2*S)) .* exp(-S.^2) .* (1/sqrt(pi)) + ...
        (S.^2 + 1 - 1./(4*S.^2)) .* erf(S) + ...
        (S/3) .* sqrt(pi) .* sqrt(T_wall/T) ...
        );

    R_s = (2 * pi .* r_theory.^2 .* P) .* ...
        ( ...
        (S + 1./(2*S)) .* exp(-S.^2) .* (1/sqrt(pi)) + ...
        (S.^2 + 1 - 1./(4*S.^2)) .* erf(S) ...
        );

    % Approximate diffuse/specular mixture.
    % Note: this is not exactly equivalent to CLL, but is useful as a
    % generalized kinetic reference.
    F_baines = alpha_CLL .* R_d + (1 - alpha_CLL) .* R_s;

    %% --------------------------------------------------------------------
    %  6. Create figure
    %  --------------------------------------------------------------------

    figure('Color', 'w', 'Position', [100, 100, 900, 640]);
    hold on; grid on; box on;

    set(gca, ...
        'XScale', 'log', ...
        'YScale', 'log', ...
        'FontSize', 12, ...
        'LineWidth', 1.0);

    xlabel('Droplet Diameter $d$ (m)', ...
        'Interpreter', 'latex', 'FontSize', 14);

    ylabel('Aerodynamic Drag Force $F_d$ (N)', ...
        'Interpreter', 'latex', 'FontSize', 14);

    title(sprintf(['\\textbf{Validation of DSMC Drag against Theoretical Regimes} ', ...
        '($v_0 = %d~\\mathrm{m/s}$)'], v_target), ...
        'Interpreter', 'latex', 'FontSize', 14);

    %% --------------------------------------------------------------------
    %  7. Plot theoretical models and DSMC data
    %  --------------------------------------------------------------------

    h1 = plot(d_theory, F_stokes, '--b', ...
        'LineWidth', 2, ...
        'DisplayName', 'Continuum no-slip (Stokes)');

    h2 = plot(d_theory, F_cunningham, '-.g', ...
        'LineWidth', 2, ...
        'DisplayName', 'Slip correction (Cunningham)');

    h3 = plot(d_theory, F_epstein, '-r', ...
        'LineWidth', 2, ...
        'DisplayName', 'Low-speed free-molecular (Epstein)');

    h4 = plot(d_theory, F_baines, ':m', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Generalized free-molecular (Baines)');

    % Discrete DSMC data only: markers and error bars, no connecting line
    h5 = errorbar(d_dsmc, F_dsmc, F_std, F_std, 'o', ...
        'LineStyle', 'none', ...
        'Color', [0 0.6 0.7], ...
        'LineWidth', 2.0, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'w', ...
        'CapSize', 6, ...
        'DisplayName', 'SPARTA DSMC data (\pm 1\sigma)');

    %% --------------------------------------------------------------------
    %  8. Fix axis limits before adding shaded regions
    %  --------------------------------------------------------------------

    xlim([1e-7, 1e-5]);

    % Let MATLAB choose y-limits from data/curves first
    drawnow;
    yl = ylim;

    %% --------------------------------------------------------------------
    %  9. Background regime shading
    %  --------------------------------------------------------------------

    p1 = patch([1e-7 1e-6 1e-6 1e-7], ...
               [yl(1) yl(1) yl(2) yl(2)], ...
               [0.90 0.95 1.00], ...
               'EdgeColor', 'none', ...
               'FaceAlpha', 0.25, ...
               'HandleVisibility', 'off');

    p2 = patch([1e-6 1e-5 1e-5 1e-6], ...
               [yl(1) yl(1) yl(2) yl(2)], ...
               [1.00 0.95 0.85], ...
               'EdgeColor', 'none', ...
               'FaceAlpha', 0.20, ...
               'HandleVisibility', 'off');

    % Send patches to the back
    uistack(p1, 'bottom');
    uistack(p2, 'bottom');

    %% --------------------------------------------------------------------
    %  10. Reference diameter lines
    %  --------------------------------------------------------------------

    for i = 1:length(d_case_labels)
        xline(d_case_labels(i), ':k', d_case_text{i}, ...
            'LabelVerticalAlignment', 'bottom', ...
            'LabelHorizontalAlignment', 'center', ...
            'FontSize', 9, ...
            'HandleVisibility', 'off');
    end

    %% --------------------------------------------------------------------
    %  11. Legend
    %  --------------------------------------------------------------------

    legend([h1, h2, h3, h4, h5], ...
        'Location', 'northwest', ...
        'Interpreter', 'latex', ...
        'FontSize', 9);

    %% --------------------------------------------------------------------
    %  12. Export figure
    %  --------------------------------------------------------------------

    set(gcf, 'PaperPositionMode', 'auto');

    filename = sprintf('validation_v%d.png', v_target);
    exportgraphics(gcf, filename, 'Resolution', 300);

end