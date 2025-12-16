clc; clear; close all;

%% Black Hole Parameters
r_s = 5;   % Schwarzschild radius

%% Figure Setup
figure('Color','k', 'Position', [100 100 1000 650]);
hold on; axis equal;
view_width = 50; view_height = 50;
axis([-view_width view_width -view_height view_height]);
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('X (units)','Color','w'); 
ylabel('Y (units)','Color','w');
title('Black Hole Light Ray Simulation','Color','w');
grid on;

%% Draw black hole & photon sphere
theta = linspace(0, 2*pi, 200);
fill(r_s*cos(theta), r_s*sin(theta), 'k', 'EdgeColor', 'c', 'LineWidth', 1.5);
plot(1.5*r_s*cos(theta), 1.5*r_s*sin(theta), 'y--', 'LineWidth', 1.5);

%% Create video writer
v = VideoWriter('blackhole_simulation2.mp4','MPEG-4');
v.FrameRate = 30;     % Adjust as desired
open(v);

%% Initial photon conditions
x = 29; 
y = 20;   
vx = -40;  
vy = -60;  

r0 = sqrt(x^2 + y^2);
phi0 = atan2(y, x);
dr0 = vx * cos(phi0) + vy * sin(phi0);
dphi0 = (-vx * sin(phi0) + vy * cos(phi0)) / r0;

X0 = [r0;phi0;dr0;dphi0];

%% Noisy geodesic wrappers (different noise per trajectory)
geod1 = @(t,X) geodesic(t, X, r_s, 1e-2, 1e-1);
geod2 = @(t,X) geodesic(t, X, r_s, 1e-2, 5e-1);
geod3 = @(t,X) geodesic(t, X, r_s, 1e-2, 10e-1);

%% Integrate
t_span = linspace(0, 100, 10000);
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);

[~, X1] = ode45(geod1, t_span, X0, opts);
[~, X2] = ode45(geod2, t_span, X0, opts);
[~, X3] = ode45(geod3, t_span, X0, opts);

%% Convert to Cartesian
[x1, y1] = pol2cart(X1(:,2), X1(:,1));
[x2, y2] = pol2cart(X2(:,2), X2(:,1));
[x3, y3] = pol2cart(X3(:,2), X3(:,1));

%% Reset figure for animation
cla;
fill(r_s*cos(theta), r_s*sin(theta), 'k', 'EdgeColor', 'r', 'LineWidth', 1.5); % black hole
plot(1.5*r_s*cos(theta), 1.5*r_s*sin(theta), 'y--', 'LineWidth', 1.5); % event horizon
plot(x, y, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

axis equal;
axis([-view_width view_width -view_height view_height]);
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('X (units)','Color','w');
ylabel('Y (units)','Color','w');
title('Light Trajectory','Color','w');
grid on;

%% Animated lines (different colors)
ray1 = animatedline('Color', [1 0.3 0.3], 'LineWidth', 1.5);
ray2 = animatedline('Color', [1 0.8 0.2], 'LineWidth', 1.5);
ray3 = animatedline('Color', [0.2 0.8 1], 'LineWidth', 1.5);

%% Animation loop
for i = 1:length(x1)
    if any(abs([x1(i) y1(i)]) > [view_width view_height])
        break
    end

    addpoints(ray1, x1(i), y1(i));
    addpoints(ray2, x2(i), y2(i));
    addpoints(ray3, x3(i), y3(i));

    drawnow limitrate;

    % ---- Capture frame for video ----
    frame = getframe(gcf);
    writeVideo(v, frame);

    pause(0.03);
end



legend("Black Hole","Photon Sphere","Start Point", ...
       "Ray 1","Ray 2","Ray 3", ...
       'TextColor','w','Color','k');
hold off;

% Close the video file
close(v);

%% around a point




%% ================== GEODESIC FUNCTION ==================
function dXdt = geodesic(~, X, r_s, noise_scale, noise_accel)

    r = X(1); 
    phi = X(2); 
    dr_dlambda = X(3); 
    dphi_dlambda = X(4);

    if r <= r_s
        dXdt = zeros(4,1);
        return
    end

    % Noise multipliers (centered at 1, small variance)
    mult_f    = normrnd(1, noise_scale);   % affects metric
    mult_acc3 = normrnd(1, noise_accel);   % affects dXdt(3)
    mult_acc4 = normrnd(1, noise_accel);   % affects dXdt(4)

    % Metric factor with mild multiplicative distortion
    f = (1 - r_s / r) * mult_f;

    % Energy-like term
    E = f * sqrt((dr_dlambda^2)/(f^2) + (r^2 * dphi_dlambda^2)/f);
    dt_dlambda = E / f;

    dXdt = zeros(4,1);

    % State equations
    dXdt(1) = dr_dlambda;
    dXdt(2) = dphi_dlambda;

    % -------------------------------
    % Radial acceleration (NOISE ADDED HERE)
    % -------------------------------
    base_dr = -(   (r_s/(2*r^2)) * f * (dt_dlambda)^2   ) ...
              +(   (r_s/(2*r^2*f)) * (dr_dlambda^2)      ) ...
              +(   (r - r_s) * (dphi_dlambda^2)          );

    dXdt(3) = base_dr * mult_acc3;   % ← MULTIPLICATIVE RANDOMNESS

    % -------------------------------
    % Angular acceleration (NOISE ADDED HERE)
    % -------------------------------
    base_dphi = -(2 * dr_dlambda * dphi_dlambda) / r;

    dXdt(4) = base_dphi * mult_acc4; % ← MULTIPLICATIVE RANDOMNESS
end
