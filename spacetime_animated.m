clc; clear; close all;
% schwarzchild radius
r_s = 5;
%% Create meshgrid for spacetime fabric
grid_size = 600; % smoothness
x_range = linspace(-30, 30, grid_size);
y_range = linspace(-30, 30, grid_size);
[X, Y] = meshgrid(x_range, y_range);
R = sqrt(X.^2 + Y.^2);
R(R < 0.1) = 0.1;
%% singularity
depth_scale = 250;
smooth_power = 1;
Z = -depth_scale ./ (R.^smooth_power);
%% heatmap plot
figure('Color', 'k', 'Position', [100 100 1200 800]);
surf(X, Y, Z, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.95);
hold on;
%% Wireframe overlay
mesh_density = 8;
surf(X(1:mesh_density:end, 1:mesh_density:end), ...
     Y(1:mesh_density:end, 1:mesh_density:end), ...
     Z(1:mesh_density:end, 1:mesh_density:end), ...
'FaceColor', 'none', 'EdgeColor', [0.7 0.7 0.7], ...
'LineWidth', 1.2, 'EdgeAlpha', 0.6);
%% event horizon
theta = linspace(0, 2*pi, 300);
eh_r = r_s;
eh_x = eh_r*cos(theta);
eh_y = eh_r*sin(theta);
eh_z = -depth_scale ./ (eh_r.^smooth_power) * ones(size(theta));
plot3(eh_x, eh_y, eh_z, 'c', 'LineWidth', 4);
%% photon sphere
ps_r = 1.5*r_s;
ps_x = ps_r*cos(theta);
ps_y = ps_r*sin(theta);
ps_z = -depth_scale ./ (ps_r.^smooth_power) * ones(size(theta));
plot3(ps_x, ps_y, ps_z, 'y--', 'LineWidth', 3);
%% Styling
colormap(turbo);
shading interp;
c = colorbar('Color', 'w');
c.Label.String = 'Gravitational Potential';
c.Label.Color = 'w';
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
xlabel('X (units)', 'Color', 'w');
ylabel('Y (units)', 'Color', 'w');
zlabel('Spacetime Curvature', 'Color', 'w');
title('Spacetime Curvature around a Black Hole', 'Color', 'w');
view(-30, 30);
axis tight;
grid on;
lighting phong;
camlight('right');
camlight('left');
material shiny;

%% Generate and animate light rays
num_rays = 8;
colors = {'r', 'g', 'b', 'm', 'w', 'c', [1 0.5 0], [0.5 1 0.5]}; % Various colors
ray_handles = cell(num_rays, 1);
marker_handles = cell(num_rays, 1);

% Pre-compute all trajectories
trajectories = cell(num_rays, 1);
for i = 1:num_rays
    % Random starting position (far from black hole)
    angle_start = rand * 2 * pi;
    r_start = 20 + rand * 8; % Between 20 and 28 units
    x0 = r_start * cos(angle_start);
    y0 = r_start * sin(angle_start);
    
    % Velocity aimed somewhat toward the black hole with variation
    target_angle = atan2(-y0, -x0) + (rand - 0.5) * pi/3;
    speed = 100 + rand * 100;
    vx = speed * cos(target_angle);
    vy = speed * sin(target_angle);
    
    % Convert to polar
    r0 = sqrt(x0^2 + y0^2);
    phi0 = atan2(y0, x0);
    dr0 = vx * cos(phi0) + vy * sin(phi0);
    dphi0 = (-vx * sin(phi0) + vy * cos(phi0)) / r0;
    X0 = [r0; phi0; dr0; dphi0];
    
    % Simulate trajectory
    t_span = [0 3];
    [~, X_traj] = ode89(@(t,X) geodesic_simple(t, X, r_s), t_span, X0);
    
    % Convert to Cartesian and calculate z
    sim_r = X_traj(:, 1);
    sim_phi = X_traj(:, 2);
    x_ray = sim_r .* cos(sim_phi);
    y_ray = sim_r .* sin(sim_phi);
    z_ray = -depth_scale ./ (sim_r.^smooth_power);
    
    trajectories{i} = [x_ray, y_ray, z_ray];
    
    % Initialize empty plot handles
    if iscell(colors{i})
        color = colors{i};
    else
        color = colors{i};
    end
    ray_handles{i} = plot3(nan, nan, nan, 'Color', color, 'LineWidth', 2.5);
    marker_handles{i} = plot3(nan, nan, nan, 'o', 'Color', color, ...
        'MarkerSize', 8, 'MarkerFaceColor', color);
end

%% Animation loop
max_points = max(cellfun(@(x) size(x, 1), trajectories));
animation_speed = 10; % Skip frames for faster animation

for frame = 1:animation_speed:max_points
    for i = 1:num_rays
        traj = trajectories{i};
        
        % Determine how many points to show
        end_idx = min(frame, size(traj, 1));
        
        if end_idx > 0
            % Update the trajectory line
            set(ray_handles{i}, 'XData', traj(1:end_idx, 1), ...
                'YData', traj(1:end_idx, 2), ...
                'ZData', traj(1:end_idx, 3));
            
            % Update the marker at the current position
            set(marker_handles{i}, 'XData', traj(end_idx, 1), ...
                'YData', traj(end_idx, 2), ...
                'ZData', traj(end_idx, 3));
        end
    end
    
    drawnow;
    pause(0.00001);
end

hold off;

%% geodesic func
function dXdt = geodesic_simple(~, X, r_s)
    r = X(1);
    phi = X(2);
    dr = X(3);
    dphi = X(4);
    
    if r <= r_s
        dXdt = zeros(4,1);
        return;
    end
    
    f = 1 - r_s/r;
    dXdt = zeros(4,1);
    dXdt(1) = dr;
    dXdt(2) = dphi;
    dXdt(3) = (r - r_s)*(dphi^2) - (r_s/(2*r^2))*dr^2;
    dXdt(4) = (-2*dr*dphi)/r;
end

%% --- Video setup ---