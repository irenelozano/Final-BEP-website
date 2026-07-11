%% ========================================================================
%  ASML Debris Mitigation - Full Kinematic Trajectory Simulator
%  Plots Velocity Decay vs Distance for ALL initial velocities and diameters
%  Reference distances: x_ref,1 = 480 mm, x_ref,2 = 1300 mm
%  ========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%  1. Load production data
%  ------------------------------------------------------------------------

filename = 'matrix_results_20260601_125031.csv';
data = readtable(filename);

diameters  = unique(data.Diameter_m);
velocities = unique(data.Velocity_ms);

% Tin (Sn) material properties
rho_sn = 6990; % Density of liquid tin [kg/m^3]

% Reference downstream distances [mm]
ref1_mm = 480;
ref2_mm = 1300;

%% ------------------------------------------------------------------------
%  2. Set up figure
%  ------------------------------------------------------------------------

figure('Color', 'w', 'Position', [50, 50, 1200, 800]);

sgtitle('\textbf{Tin Debris Velocity Decay over Distance in 100 Pa $H_2$}', ...
    'Interpreter', 'latex', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

% Color palette for the initial velocities
colors = turbo(length(velocities));

%% ------------------------------------------------------------------------
%  3. F0 background correction
%  ------------------------------------------------------------------------

[unique_D_f0, ~, ~] = unique(data.Diameter_m);

F0_map = containers.Map('KeyType', 'double', 'ValueType', 'double');

for di = 1:length(unique_D_f0)

    idx_f0 = data.Diameter_m == unique_D_f0(di);

    v_f0 = data.Velocity_ms(idx_f0);
    F_f0 = abs(data.Drag_Fx(idx_f0));

    % Estimate residual zero-velocity force offset
    p_f0 = polyfit(v_f0, F_f0, 1);
    F0_map(unique_D_f0(di)) = p_f0(2);
end

fprintf('Calculating kinematic trajectories for all cases...\n');

%% ------------------------------------------------------------------------
%  4. Loop through each diameter
%  ------------------------------------------------------------------------

for i = 1:length(diameters)

    d = diameters(i);

    % Create a 2x3 grid of subplots
    subplot(2, 3, i);
    hold on; grid on; box on;

    % Droplet mass
    volume = (4/3) * pi * (d/2)^3;
    mass = rho_sn * volume;

    % Extract force-velocity data for this diameter
    idx = data.Diameter_m == d;

    V_data = data.Velocity_ms(idx);
    F_data = abs(data.Drag_Fx(idx)) - F0_map(d);

    % Sort data
    [V_data, sort_idx] = sort(V_data);
    F_data = F_data(sort_idx);

    % Remove non-positive corrected-force points
    valid_mask = F_data > 0;
    V_valid = V_data(valid_mask);
    F_valid = F_data(valid_mask);

    if isempty(V_valid)
        warning('All corrected forces are non-positive for d = %.1e m. Skipping diameter.', d);
        continue;
    end

    %% --------------------------------------------------------------------
    %  Fit drag model: F = a v + b v^2
    %  --------------------------------------------------------------------

    if length(V_valid) >= 3

        % Fit F/v = a + b v
        coeffs = polyfit(V_valid, F_valid ./ V_valid, 1);

        a_fit = max(0, coeffs(2));
        b_fit = max(0, coeffs(1));

    elseif length(V_valid) == 2

        % Two points: use linear model through origin
        a_fit = F_valid(2) / V_valid(2);
        b_fit = 0;

    else

        % One point: use linear model through origin
        a_fit = F_valid(1) / V_valid(1);
        b_fit = 0;
    end

    fprintf('d = %.3g um: F0 = %.2e N, a = %.2e, b = %.2e, tau = %.1f ms\n', ...
        d*1e6, F0_map(d), a_fit, b_fit, mass/a_fit*1000);

    %% --------------------------------------------------------------------
    %  Integrate trajectories
    %  --------------------------------------------------------------------

    % Characteristic decay time
    tau = mass / max(a_fit, 1e-30);

    % Simulate long enough to capture decay
    t_span = [0, max(50e-3, 12*tau)];

    for j = 1:length(velocities)

        v0 = velocities(j);

        % Differential equation:
        % Y(1) = x, Y(2) = v
        ode_func = @(t, Y) [
            Y(2);
            -(a_fit * max(0, Y(2)) + b_fit * max(0, Y(2))^2) / mass
        ];

        % Stop once velocity drops below threshold
        options = odeset('Events', @(t,Y) stopEvent(t,Y));

        % Initial conditions
        Y0 = [0; v0];

        % Solve ODE
        [~, Y_out] = ode45(ode_func, t_span, Y0, options);

        % Extract distance and velocity
        x_mm = Y_out(:,1) * 1000;
        v_current = Y_out(:,2);

        % Plot trajectory
        plot(x_mm, v_current, ...
            'LineWidth', 2, ...
            'Color', colors(j,:), ...
            'DisplayName', sprintf('$v_0 = %g$ m/s', v0));
    end

    %% --------------------------------------------------------------------
    %  Format subplot
    %  --------------------------------------------------------------------

    title(sprintf('$d = %g$ $\\mu$m', d*1e6), ...
        'Interpreter', 'latex', ...
        'FontSize', 12);

    xlabel('Distance $x$ (mm)', ...
        'Interpreter', 'latex');

    ylabel('Velocity $v$ (m/s)', ...
        'Interpreter', 'latex');

    % Extend x-axis to include both reference distances
    xlim([0, 1400]);
    ylim([0, 750]);

    % Reference distance 1: 0.48 m
    xline(ref1_mm, '--', '$x_{\mathrm{ref},1}=0.48\,\mathrm{m}$', ...
        'Color', [0.85 0.35 0], ...
        'LineWidth', 2, ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelHorizontalAlignment', 'left', ...
        'FontSize', 10, ...
        'Interpreter', 'latex');

    % Reference distance 2: 1.30 m
    xline(ref2_mm, '--b', '$x_{\mathrm{ref},2}=1.30\,\mathrm{m}$', ...
        'LineWidth', 2, ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelHorizontalAlignment', 'left', ...
        'FontSize', 10, ...
        'Interpreter', 'latex');

    % Only show legend on first subplot
    if i == 1
        legend('Location', 'northeast', ...
            'Interpreter', 'latex', ...
            'FontSize', 9);
    end
end

%% ------------------------------------------------------------------------
%  5. ODE stop event
%  ------------------------------------------------------------------------

function [value, isterminal, direction] = stopEvent(~, Y)
    value = Y(2) - 0.1; % Stop when velocity drops below 0.1 m/s
    isterminal = 1;     % Halt integration
    direction = -1;
end