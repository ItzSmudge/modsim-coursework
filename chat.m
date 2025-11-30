clc; clear; close all;

%% Black Hole Parameters
r_s = 5;   % Schwarzschild radius

%% Create 12 launch directions (0°, 30°, …, 330°)
N = 20;
angles = linspace(0, 360-(360/N), N) * pi/180;

%% Figure Setup
figure('Color','k', 'Position', [100 100 1000 650]);
hold on; axis equal;
view_width = 50; view_height = 50;
axis([-view_width view_width -view_height view_height]);
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('X (units)','Color','w'); 
ylabel('Y (units)','Color','w');
title(string(N) + " Noisy Light Rays Around a Point", 'Color','w');
grid on;

%% Draw black hole & photon sphere
theta = linspace(0, 2*pi, 200);
fill(r_s*cos(theta), r_s*sin(theta), 'k', 'EdgeColor', 'r', 'LineWidth', 1.5);
plot(1.5*r_s*cos(theta), 1.5*r_s*sin(theta), 'y--', 'LineWidth', 1.5);

%% Starting point for all rays
x0 = 29; 
y0 = -20;
plot(x0, y0, 'go', 'MarkerSize', 10, 'MarkerFaceColor','g');

%% Speed magnitude
speed = 80;

%% Time settings
t_span = linspace(0, 100, 5000);
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);

%% Noise parameters
noise_scale  = 1e-2;   % mild metric noise
noise_accel  = 10e-1;   % stronger acceleration noise

%% Colors for rays
colors = hsv(N);

%% Precompute all trajectories
X_all = cell(N,1);
for k = 1:N

    % Velocity components
    vx = speed * cos(angles(k));
    vy = speed * sin(angles(k));

    % Convert initial state to (r,phi,dr,dphi)
    r0 = sqrt(x0^2 + y0^2);
    phi0 = atan2(y0, x0);
    dr0 = vx*cos(phi0) + vy*sin(phi0);
    dphi0 = (-vx*sin(phi0) + vy*cos(phi0)) / r0;

    X0 = [r0; phi0; dr0; dphi0];

    % Unique noise per ray
    geod = @(t,X) geodesic_noisy(t, X, r_s, noise_scale, noise_accel);

    % Integrate geodesic
    [~, X] = ode45(geod, t_span, X0, opts);

    % Save
    X_all{k} = X;
end

%% Convert all trajectories to Cartesian
xy = cell(N,2);
for k = 1:N
    [xx, yy] = pol2cart(X_all{k}(:,2), X_all{k}(:,1));
    xy{k,1} = xx;
    xy{k,2} = yy;
end

%% Animated lines (one per ray)
ray = cell(N,1);
for k = 1:N
    ray{k} = animatedline('Color', colors(k,:), 'LineWidth', 1.5);
end

%% Animation loop — reveal points step-by-step
max_len = max(cellfun(@(a) length(a(:,1)), X_all));

for i = 1:max_len
    for k = 1:N
        
        if i <= length(xy{k,1})
            xk = xy{k,1}(i);
            yk = xy{k,2}(i);

            if abs(xk) <= view_width && abs(yk) <= view_height
                addpoints(ray{k}, xk, yk);
            end
        end

    end

    drawnow;   % Remove limitrate for smooth playback
    pause(0.01);
end





%% ================== NOISY GEODESIC FUNCTION ==================
function dXdt = geodesic_noisy(~, X, r_s, noise_scale, noise_accel)

    r = X(1); 
    phi = X(2); 
    dr_dlambda = X(3); 
    dphi_dlambda = X(4);

    if r <= r_s
        dXdt = zeros(4,1);
        return
    end

    % Noise multipliers
    mult_f    = normrnd(1, noise_scale);
    mult_acc3 = normrnd(1, noise_accel);
    mult_acc4 = normrnd(1, noise_accel);

    % Metric factor with noise
    f = (1 - r_s / r) * mult_f;

    % Energy-like quantity
    E = f * sqrt((dr_dlambda^2)/(f^2) + (r^2 * dphi_dlambda^2)/f);
    dt_dlambda = E / f;

    dXdt = zeros(4,1);

    % State evolution
    dXdt(1) = dr_dlambda;
    dXdt(2) = dphi_dlambda;

    % Radial acceleration (noisy)
    base_dr = -( (r_s/(2*r^2)) * f * (dt_dlambda)^2 ) ...
              + ( (r_s/(2*r^2*f)) * (dr_dlambda^2) ) ...
              + ( (r - r_s) * (dphi_dlambda^2) );

    dXdt(3) = base_dr * mult_acc3;

    % Angular acceleration (noisy)
    base_dphi = -(2 * dr_dlambda * dphi_dlambda) / r;
    dXdt(4) = base_dphi * mult_acc4;
end