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

%% singulairty
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

