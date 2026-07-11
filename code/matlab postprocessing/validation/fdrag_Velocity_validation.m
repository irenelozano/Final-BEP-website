%% ========================================================================
%  TU/e Mechanical Engineering - BEP DSMC Validation Post-Processor
%  Section 2.4: Theoretical Drag Validation (Force vs. Velocity)
%  ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
%  1. Load DSMC Matrix Data
%  ------------------------------------------------------------------------
data_file = 'matrix_results_20260601_125031.csv';
if ~isfile(data_file)
    error('Matrix results CSV not found! Check file name and path.');
end
data = readtable(data_file, 'VariableNamingRule', 'preserve');

%% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%  Background force correction: subtract the v=0 offset per diameter
%  (from xlo outflow+inlet pressure imbalance)
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
        label = sprintf('%.0f nm', unique_D(i)*1e9);
    else
        label = sprintf('%.0f um', unique_D(i)*1e6);
    end
    fprintf('  d=%s: F0 = %.4e N\n', label, F0_map(unique_D(i)));
end
fprintf('---\n\n');

% 2. Physical Constants & Gas Properties (100 Pa H2 at 293 K)
%  ------------------------------------------------------------------------
T = 293.15;              % Gas temperature [K]
T_wall = T;              % Isothermal droplet wall [K]
P = 100;                 % Gas pressure [Pa]
R = 8.314;               % Universal gas constant [J/(mol K)]
M_H2 = 2.016e-3;         % H2 molar mass [kg/mol]
mu_H2 = 8.8e-6;          % Dynamic viscosity of H2 [Pa s]
k_B = 1.38e-23;          % Boltzmann constant [J/K]
m_H2 = 3.34e-27;         % H2 molecular mass [kg]
d_m_H2 = 2.72e-10;       % Molecular diameter of H2 (from hydrogen.vhs) [m]
n_H2   = 2.47e22;        % Number density (from in.drag: global nrho 2.47e22) [1/m^3]
lambda = 1 / (sqrt(2) * pi * d_m_H2^2 * n_H2); % Mean free path [m]
alpha_CLL = 0.3;         % Accommodation coefficient (diffuse/specular ratio)

rho_gas = (P * M_H2) / (R * T);             % Gas density [kg/m^3]
v_mean = sqrt((8 * k_B * T) / (pi * m_H2)); % Mean molecular speed [m/s]
v_beta = sqrt((2 * k_B * T) / m_H2);        % Most probable molecular speed [m/s]

% Davies (1945) Slip Correction Constants
alpha_c = 1.257;
beta_c  = 0.400;
gamma_c = 1.100;

%% ------------------------------------------------------------------------
%  3. Setup Theoretical Velocity Range
%  ------------------------------------------------------------------------
v_theory = logspace(log10(10), log10(800), 200); % Continuous velocity array

% Target diameters to highlight the extremes of the matrix
target_D = [100e-9, 10e-6]; 
panel_titles = {'(a) Free-Molecular (100 nm)', '(b) Free-Molecular (10 \mum)'};

fig = figure('Name', 'Section 2.4: Force vs Velocity Validation', 'Position', [100, 100, 1100, 550], 'Color', 'w');

%% ------------------------------------------------------------------------
%  4. Loop over the two target diameters
%  ------------------------------------------------------------------------
for i = 1:length(target_D)
    D_p = target_D(i);
    r_p = D_p / 2;
    Kn_regime = lambda / D_p;        % diameter-based: for the displayed "Kn ≈ ..." annotation (Table 1.1)
    Kn_slip   = 2 * lambda / D_p;    % radius-based: REQUIRED by Davies coefficients (Eq. 2.18)
    
    % --- A. Empirical Continuum Models ---
    % 1. Classic Stokes (No-slip)
    F_stokes = 3 * pi * mu_H2 * D_p .* v_theory;
    
    % 2. Cunningham-Knudsen Slip Correction
    C_Kn = 1 + Kn_slip * (alpha_c + beta_c * exp(-gamma_c / Kn_slip));
    F_corrected = F_stokes ./ C_Kn;
    
    % --- B. Kinetic Models ---
    % 3. Epstein Low-Speed Limit
    F_epstein = (4/3) * pi * rho_gas * (r_p^2) * v_mean .* v_theory * (1 + (pi/8) * alpha_CLL);
    
    % 4. Baines Generalized Free-Molecular
    S = v_theory ./ v_beta;
    base = (S + 1./(2*S)).*exp(-S.^2)./sqrt(pi) + (S.^2 + 1 - 1./(4*S.^2)).*erf(S);
    R_d  = 2*pi * r_p^2 * P .* ( base + (S/3).*sqrt(pi).*sqrt(T_wall/T) );
    R_s  = 2*pi * r_p^2 * P .* ( base );
    F_baines = alpha_CLL .* R_d + (1 - alpha_CLL) .* R_s;
        
  
    % --- C. Extract DSMC Data with F0 subtraction ---
    idx = (data.Diameter_m == D_p);
    v_dsmc = data.Velocity_ms(idx);
    F_dsmc = abs(data.Drag_Fx(idx)) - F0_map(D_p);
    
    % Check for standard deviation in CSV, else default to 2% noise floor
    if ismember('Drag_Fx_Std', data.Properties.VariableNames)
        F_std = data.Drag_Fx_Std(idx);
    else
        F_std = F_dsmc * 0.02; 
    end
    
    % Sort for clean plotting
    [v_dsmc, sort_idx] = sort(v_dsmc);
    F_dsmc = F_dsmc(sort_idx);
    F_std = F_std(sort_idx);

    % --- D. Plotting Subplot ---
    ax = subplot(1, 2, i);
    hold on; grid on; box on;
    set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 11);
    
    % Plot Theoretical Curves
    plot(v_theory, F_stokes, '--b', 'LineWidth', 2, 'DisplayName', 'Stokes (No-Slip)');
    plot(v_theory, F_corrected, '-.g', 'LineWidth', 2, 'DisplayName', 'Cunningham-Weber Corrected');
    plot(v_theory, F_baines, '-m', 'LineWidth', 1.5, 'DisplayName', 'Baines (generalized free-molecular)');

    plot(v_theory, F_epstein, ':r', 'LineWidth', 2.5, 'DisplayName', 'Epstein (Low-Speed Kinetic)');
    %plot(v_theory, F_baines, '-m', 'LineWidth', 1.5, 'DisplayName', 'Baines (Generalized Kinetic)');
    
    % Plot DSMC Data with Error Bars
    errorbar(v_dsmc, F_dsmc, F_std, F_std, 'ks', 'MarkerSize', 7, 'MarkerFaceColor', [0 0.6 0.7], ...
        'LineWidth', 1.5, 'DisplayName', 'SPARTA DSMC (\pm 1\sigma)');

    % Formatting
    xlabel('Flow Velocity $v_0$ (m/s)', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Aerodynamic Drag Force $F_x$ (N)', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
    title(panel_titles{i}, 'FontSize', 13);
    
    % Annotate Knudsen Number
    text(0.05, 0.90, sprintf('$Kn \\approx %.1f$', Kn_regime), 'Units', 'normalized', ...
        'Interpreter', 'latex', 'FontSize', 14, 'BackgroundColor', 'w', 'EdgeColor', 'k');
    
    % Add limits padding
    xlim([15, 800]);
    xticks([20 100 250 400 550 700]);
    xticklabels({'20', '100', '250', '400', '550', '700'});
    
    if i == 1
        legend('Location', 'northwest', 'FontSize', 9);
    end
end

sgtitle('\textbf{Validation of DSMC Drag vs. Continuum and Kinetic Theories}', ...
    'Interpreter', 'latex', 'FontSize', 16);
for i = 1:length(unique_D)
    d = unique_D(i); r = d/2;
    idx = data.Diameter_m==d;
    v = data.Velocity_ms(idx);
    F = abs(data.Drag_Fx(idx)) - F0_map(d);
    p_all  = polyfit(v, F, 1);                  a_all  = p_all(1);          % fit using all 6 points
    hi = v>=250;                                                            % high-S, reliable only
    p_hi   = polyfit(v(hi), F(hi), 1);          a_hi   = p_hi(1);
    a_E = (4/3)*pi*rho_gas*r^2*v_mean*(1+(pi/8)*alpha_CLL);
    fprintf('d=%6.0fnm:  a_all/a_E=%.2f   a_hi/a_E=%.2f\n', d*1e9, a_all/a_E, a_hi/a_E);
end