%% ========================================================================
%  Stopping Distance vs Initial Velocity — Realistic Tin Debris (≤1 µm)
%  Reference distances: ref_1 = 480 mm, ref_2 = 1300 mm
%  ========================================================================

clear; clc; close all;

filename = 'matrix_results_20260601_125031.csv';
data = readtable(filename);

diameters = unique(data.Diameter_m);
diameters = diameters(diameters <= 1e-6);
velocities = unique(data.Velocity_ms);

rho_sn = 6990; % Liquid tin density [kg/m^3]

% Reference distances [m]
ref_dist  = [0.48, 1.30];
ref_label = {'0.48 m', '1.30 m'};

% --- F0 background correction ---
[unique_D_f0, ~, ~] = unique(data.Diameter_m);
F0_map = containers.Map('KeyType', 'double', 'ValueType', 'double');

for di = 1:length(unique_D_f0)
    idx_f0 = data.Diameter_m == unique_D_f0(di);
    v_f0 = data.Velocity_ms(idx_f0);
    F_f0 = abs(data.Drag_Fx(idx_f0));

    p_f0 = polyfit(v_f0, F_f0, 1);
    F0_map(unique_D_f0(di)) = p_f0(2);
end

v_range = linspace(0, max(velocities), 200);
v_range(1) = 1e-6;

colors = lines(length(diameters));

figure('Color', 'w', 'Position', [100, 100, 800, 600]);
hold on; grid on; box on;

fit_results = cell(length(diameters), 1);

for i = 1:length(diameters)

    d = diameters(i);

    volume = (4/3) * pi * (d/2)^3;
    mass = rho_sn * volume;

    idx = (data.Diameter_m == d);

    V_data = data.Velocity_ms(idx);
    F_data = abs(data.Drag_Fx(idx)) - F0_map(d);

    [V_data, sort_idx] = sort(V_data);
    F_data = F_data(sort_idx);

    valid_mask = (F_data > 0);
    V_valid = V_data(valid_mask);
    F_valid = F_data(valid_mask);

    if length(V_valid) >= 3
        coeffs = polyfit(V_valid, F_valid ./ V_valid, 1);
        a_fit = max(0, coeffs(2));
        b_fit = max(0, coeffs(1));
    elseif length(V_valid) == 2
        a_fit = F_valid(2) / V_valid(2);
        b_fit = 0;
    else
        a_fit = F_valid(1) / V_valid(1);
        b_fit = 0;
    end

    if b_fit > 0
        x_stop = (mass / b_fit) * log(1 + b_fit * v_range / a_fit);
    else
        x_stop = (mass / a_fit) * v_range;
    end

    plot(v_range, x_stop * 1000, ...
        'LineWidth', 2.5, ...
        'Color', colors(i,:), ...
        'DisplayName', sprintf('$d = %g$ nm', d * 1e9));

    x_stop_points = interp1(v_range, x_stop, V_valid);

    plot(V_valid, x_stop_points * 1000, 'o', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', 'k', ...
        'HandleVisibility', 'off');

    fit_results{i} = struct( ...
        'd', d, ...
        'mass', mass, ...
        'a_fit', a_fit, ...
        'b_fit', b_fit);
end

%% ------------------------------------------------------------------------
%  Reference distance lines
%  ------------------------------------------------------------------------

yline(ref_dist(1)*1000, '--', ...
    'LineWidth', 2, ...
    'Color', [0.8 0.4 0], ...
    'Label', sprintf('$\\mathrm{ref}_1 = %s$', ref_label{1}), ...
    'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 11, ...
    'Interpreter', 'latex');

yline(ref_dist(2)*1000, '--r', ...
    'LineWidth', 2, ...
    'Label', sprintf('$\\mathrm{ref}_2 = %s$', ref_label{2}), ...
    'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 11, ...
    'Interpreter', 'latex');

%% ------------------------------------------------------------------------
%  Formatting
%  ------------------------------------------------------------------------

xlabel('Initial Velocity $v_0$ (m/s)', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);

ylabel('Stopping Distance $x_{\rm stop}$ (mm)', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);

title('Tin Debris Stopping Distance in 100 Pa H$_2$', ...
    'Interpreter', 'latex', ...
    'FontSize', 16);

% Shaded region below ref_2 = 1.30 m
xl = xlim;
yl = ylim;

hShade = fill([0, xl(2), xl(2), 0], ...
              [0, 0, ref_dist(2)*1000, ref_dist(2)*1000], ...
              [0.8 0.9 1.0], ...
              'FaceAlpha', 0.25, ...
              'EdgeColor', 'none', ...
              'HandleVisibility', 'off');

% Send shading behind lines and markers
uistack(hShade, 'bottom');

ylim([0, yl(2)]);

legend('Location', 'northwest', ...
    'Interpreter', 'latex', ...
    'FontSize', 12);

%% ------------------------------------------------------------------------
%  Summary table for both reference distances
%  ------------------------------------------------------------------------

fprintf('--- Stopping Distance & Velocity at Reference Distances ---\n');

fprintf('%5s  %5s  %8s', 'd', 'v0', 'x_stop');

for ri = 1:length(ref_dist)
    fprintf('  v@%.2fm', ref_dist(ri));
end

fprintf('\n');

fprintf('%5s  %5s  %8s', '---', '---', '------');

for ri = 1:length(ref_dist)
    fprintf('  ------');
end

fprintf('\n');

for i = 1:length(diameters)

    d = fit_results{i}.d;
    mass = fit_results{i}.mass;
    a_fit = fit_results{i}.a_fit;
    b_fit = fit_results{i}.b_fit;

    for v = velocities'

        if b_fit > 0
            x_stop = (mass / b_fit) * log(1 + b_fit * v / a_fit);
        else
            x_stop = (mass / a_fit) * v;
        end

        fprintf('%5.0f nm  %5.0f  %5.0f mm', d*1e9, v, x_stop*1000);

        for ri = 1:length(ref_dist)

            if x_stop < ref_dist(ri)
                fprintf('  stop  ');
            else
                if b_fit > 0
                    v_at = (1/b_fit) * ((a_fit + b_fit*v) * ...
                        exp(-b_fit*ref_dist(ri)/mass) - a_fit);
                else
                    v_at = v - (a_fit/mass) * ref_dist(ri);
                end

                fprintf('  %4.0f m/s', v_at);
            end
        end

        fprintf('\n');
    end
end

%% ------------------------------------------------------------------------
%  Critical velocities
%  ------------------------------------------------------------------------

fprintf('\n--- Critical Velocity to Reach Each Reference ---\n');

for ri = 1:length(ref_dist)

    fprintf('At %s:\n', ref_label{ri});

    for i = 1:length(diameters)

        d = fit_results{i}.d;
        mass = fit_results{i}.mass;
        a_fit = fit_results{i}.a_fit;
        b_fit = fit_results{i}.b_fit;

        if b_fit > 0
            v_crit = (a_fit/b_fit) * (exp(ref_dist(ri)*b_fit/mass) - 1);
        else
            v_crit = ref_dist(ri) * a_fit / mass;
        end

        fprintf('  d = %.0f nm: v_crit = %.0f m/s\n', d*1e9, v_crit);
    end
end

fprintf('---\n');