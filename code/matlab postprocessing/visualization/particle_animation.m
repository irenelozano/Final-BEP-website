%% ========================================================================
%  SPARTA Flow Visualization — 2D Cross-Section of DSMC Simulation
%  Shows: Cartesian grid, H₂ molecules flowing through cells,
%         surface collisions on the tin sphere (no intermolecular collisions)
%  ========================================================================
clear; clc; close all;

%% Domain (2D slice through the 3D box)
Lx = 10; Ly = 10;        % domain size [arb. units]
cx = 5; cy = 5; r = 1.5; % tin sphere (cross-section)
N_grid = 20;             % cells per side
dx_cell = Lx / N_grid;

% Flow parameters
v_flow = 1.0;             % mean drift (x-direction)
v_therm = 2.5;            % thermal speed
alpha = 0.3;              % CLL accommodation coefficient
dt = 0.04;                % timestep per frame

% Injection: particles per timestep at steady state ≈ injection_rate
% Crossing time = Lx/v_flow = 10 steps at v_flow=1 with dt=0.04 → 250 steps
% Target ~600 particles in domain → inject ~600/250 ≈ 2-3 per step
injection_rate = 3;

% Dynamic particle array (grows from empty)
x = []; y = [];
vx = []; vy = [];

%% Setup figure
fig = figure('Color', 'w', 'Position', [100, 100, 900, 700]);
ax = axes('Position', [0.08, 0.08, 0.85, 0.85]);
hold(ax, 'on'); axis equal; axis([0, Lx, 0, Ly]);
box on; set(ax, 'XTick', [], 'YTick', []);

% 1. Draw grid cells
for i = 0:N_grid
    plot(ax, [i*Lx/N_grid, i*Lx/N_grid], [0, Ly], 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
    plot(ax, [0, Lx], [i*Ly/N_grid, i*Ly/N_grid], 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
end

% 2. Draw tin sphere
theta = linspace(0, 2*pi, 80);
fill(ax, cx + r*cos(theta), cy + r*sin(theta), [0.5 0.7 0.9], ...
    'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.3);
text(ax, cx, cy, 'Tin', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14);

% 3. Boundary labels
text(ax, -0.3, Ly/2, 'Inflow', 'Rotation', 90, 'FontSize', 11, 'Color', [0.5 0.5 0.5]);
text(ax, Lx+0.3, Ly/2, 'Outflow', 'Rotation', 90, 'FontSize', 11, 'Color', [0.5 0.5 0.5]);
text(ax, Lx/2, -0.3, 'Symmetry (reflect)', 'FontSize', 10, 'Color', [0.5 0.5 0.5], ...
    'HorizontalAlignment', 'center');
text(ax, Lx/2, Ly+0.3, 'Symmetry (reflect)', 'FontSize', 10, 'Color', [0.5 0.5 0.5], ...
    'HorizontalAlignment', 'center');

% 4. Particle scatter (initially empty)
particle_h = scatter(ax, [], [], 8, [], 'filled');
colormap(ax, jet);
caxis(ax, [0, v_flow + v_therm]);
cb = colorbar(ax, 'Position', [0.92, 0.08, 0.025, 0.85]);
ylabel(cb, 'x-velocity', 'FontSize', 10);

% 5. Info text
info_h = text(ax, 0.02, 0.97, 'Steps: 0  |  Particles: 0', ...
    'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k');

title(ax, 'SPARTA DSMC: H$_2$ Flow Around Tin Sphere (2D Slice)', ...
    'Interpreter', 'latex', 'FontSize', 14);

%% Animation loop
for step = 1:3000
    % --- Ballistic advection ---
    x = x + vx * dt;
    y = y + vy * dt;

    % --- Outflow boundaries (oo): particles that reach either x-face leave ---
    exit_xhi = (x > Lx);
    exit_xlo = (x < 0) | (x < dx_cell & vx < 0);
    exit_mask = (exit_xhi | exit_xlo);
    x(exit_mask) = []; y(exit_mask) = [];
    vx(exit_mask) = []; vy(exit_mask) = [];

    % --- Symmetry (rr): specular reflection at y-faces ---
    reflect_ylo = (y < 0);
    reflect_yhi = (y > Ly);
    y(reflect_ylo) = -y(reflect_ylo);
    vy(reflect_ylo) = -vy(reflect_ylo);
    y(reflect_yhi) = 2*Ly - y(reflect_yhi);
    vy(reflect_yhi) = -vy(reflect_yhi);

    % --- Surface collisions with tin sphere ---
    dx_s = x - cx;
    dy_s = y - cy;
    dist = sqrt(dx_s.^2 + dy_s.^2);
    hit = (dist < r);
    if any(hit)
        hit_idx = find(hit)';
        for idx = hit_idx
            nx = (x(idx) - cx) / dist(idx);
            ny = (y(idx) - cy) / dist(idx);
            x(idx) = cx + nx * (r + 0.02);
            y(idx) = cy + ny * (r + 0.02);
            if rand < alpha
                angle = rand * 2*pi;
                speed = sqrt(vx(idx)^2 + vy(idx)^2);
                vx(idx) = speed * cos(angle);
                vy(idx) = speed * sin(angle);
            else
                vn = vx(idx)*nx + vy(idx)*ny;
                vx(idx) = vx(idx) - 2*vn*nx;
                vy(idx) = vy(idx) - 2*vn*ny;
            end
        end
    end

    % --- Inflow (fix in emit/face xlo): inject fresh particles ---
    n_inject = injection_rate;
    x_new = zeros(n_inject, 1);
    y_new = rand(n_inject, 1) * Ly;
    vx_new = v_flow + v_therm * randn(n_inject, 1);
    vy_new = v_therm * randn(n_inject, 1);

    x = [x; x_new]; y = [y; y_new];
    vx = [vx; vx_new]; vy = [vy; vy_new];

    % --- Update plot ---
    set(particle_h, 'XData', x, 'YData', y, 'CData', vx);
    set(info_h, 'String', sprintf('Steps: %d  |  Particles: %d', step, length(x)));
    drawnow;
    pause(0.03);
end
