clc; clear; close all;

%% Black Hole Parameters
r_s = 5;   % Schwarzschild radius


%% Number of rays
N = 75;
angles = linspace(0, 360-(360/N), N) * pi/180;


%% Video Writer
v = VideoWriter('light_rays_line.mp4', 'MPEG-4');
v.FrameRate = 5;
open(v);
%% Figure Setup
fig = figure('Color','k', 'Position', [100 100 1000 650]);
hold on; axis equal;
view_width = 50; 
view_height = 50;
axis([-view_width view_width -view_height view_height]);
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('X (units)','Color','w'); 
ylabel('Y (units)','Color','w');
title(sprintf('%d Noisy Light Rays Around a Moving Point', N), 'Color','w');
grid on;

%% Draw black hole & photon sphere
theta = linspace(0, 2*pi, 200);
fill(r_s*cos(theta), r_s*sin(theta), 'k', 'EdgeColor', 'r', 'LineWidth', 1.5);
plot(1.5*r_s*cos(theta), 1.5*r_s*sin(theta), 'y--', 'LineWidth', 1.5);

%% MOVING SOURCE PATH
T = 10;                    % number of video frames
tpath = linspace(-10, -9, T);
% path_x = 25*cos(tpath) - 10;       % circular
% path_y = 25*sin(tpath);            % circular
path_y = repmat(25, size(tpath)); % constant y position
path_x = tpath;

%% Speed magnitude
speed = 80;

%% Time settings
t_span = linspace(0, 100, 10000);
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);

%% Noise parameters
noise_scale  = 1e-2;
noise_accel  = 1e0;

%% Colors
colors = hsv(N);


%% MAIN LOOP — ONE FRAME PER SOURCE POSITION
for frame = 1:T
    frame
    
    cla; hold on;
    axis([-view_width view_width -view_height view_height]);
    
    % redraw static objects
    fill(r_s*cos(theta), r_s*sin(theta), 'k', 'EdgeColor', 'r', 'LineWidth', 1.5);
    plot(1.5*r_s*cos(theta), 1.5*r_s*sin(theta), 'y--', 'LineWidth', 1.5);

    % current source position
    x0 = path_x(frame);
    y0 = path_y(frame);
    plot(x0, y0, 'go', 'MarkerSize', 10, 'MarkerFaceColor','g');

    % compute and draw all rays at this frame
    for k = 1:N
        
        vx = speed * cos(angles(k));
        vy = speed * sin(angles(k));
        
        % convert to (r,phi)
        r0 = sqrt(x0^2 + y0^2);
        phi0 = atan2(y0, x0);

        dr0 = vx*cos(phi0) + vy*sin(phi0);
        dphi0 = (-vx*sin(phi0) + vy*cos(phi0)) / r0;

        X0 = [r0; phi0; dr0; dphi0];

        geod = @(t,X) geodesic_noisy(t, X, r_s, noise_scale, noise_accel);
        
        [~, X] = ode45(geod, t_span, X0, opts);

        [xx, yy] = pol2cart(X(:,2), X(:,1));

        plot(xx, yy, 'Color', colors(k,:), 'LineWidth', 0.7);
    end

    % write frame
    frame_img = getframe(fig);
    writeVideo(v, frame_img);
end

close(v);
disp("Video saved as light_rays.mp4");


%% =========== NOISY GEODESIC FUNCTION ===========
function dXdt = geodesic_noisy(~, X, r_s, noise_scale, noise_accel)

    r = X(1); 
    phi = X(2); 
    dr = X(3); 
    dphi = X(4);

    if r <= r_s
        dXdt = zeros(4,1);
        return
    end

    % noise
    mult_f = normrnd(1, noise_scale);
    mult_acc3 = normrnd(1, noise_accel);
    mult_acc4 = normrnd(1, noise_accel);

    f = (1 - r_s/r) * mult_f;

    E = f * sqrt((dr^2)/(f^2) + (r^2 * dphi^2)/f);
    dt = E / f;

    dXdt = zeros(4,1);

    dXdt(1) = dr;
    dXdt(2) = dphi;

    base_dr = -( (r_s/(2*r^2)) * f * dt^2 ) ...
              + ( (r_s/(2*r^2*f)) * (dr^2) ) ...
              + ( (r - r_s) * (dphi^2) );

    dXdt(3) = base_dr * mult_acc3;

    base_dphi = -(2 * dr * dphi) / r;
    dXdt(4) = base_dphi * mult_acc4;
end
