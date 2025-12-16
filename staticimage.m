%% Modelling a black hole and light rays in matlab (coursework file)
clc; clear;
% assuming G=c=1
% geometrized units (everything has the same unit)
r_s = 5; % means 5 units of length instead of 5m or 5Km

%% Set up figure
figure('Color','k', 'Position', [100 100 1000 650]);
hold on; axis equal;
view_width = 50; view_height = 50;
axis([-view_width view_width -view_height view_height]);
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('X (units)','Color','w');
ylabel('Y (units)','Color','w');
title('Black Hole Light Ray Simulation','Color','w');
grid on

%% Plot black hole
theta = linspace(0, 2*pi, 100);
% Event Horizon
fill(r_s*cos(theta), r_s*sin(theta), 'k', 'EdgeColor', 'c', 'LineWidth', 1.5);
% Photon Sphere (r = 1.5 * rs)
% light is most unstable in the photon sphere
plot(1.5*r_s*cos(theta), 1.5*r_s*sin(theta), 'y--', 'LineWidth', 1.5);

%% ============= Initial conditions ================
x = -30;  % horizontal position
y = 13;   % vertical position
vx = 185; % horizontal velocity
vy = -1;  % vertical velocity

% convert into polar coordinates
r0 = sqrt(x^2 + y^2);
phi0 = atan2(y, x);
dr0 = vx * cos(phi0) + vy * sin(phi0);
dphi0 = (-vx * sin(phi0) + vy * cos(phi0)) / r0;

% initial state array
X0 = [r0; phi0; dr0; dphi0];

%% Noise parameters (adjust these to control noise strength)
noise_scale = 1e-2;  % affects metric factor (smaller = less noise)
noise_accel = 1e-1;  % affects accelerations (smaller = less noise)

%% Simulation
t_span = [0 1000]; % longer simulation means more of the ray of light is drawn
geodesic_equation = @(t,X) geodesic(t, X, r_s, noise_scale, noise_accel);
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4);
[sim_time, X] = ode45(geodesic_equation, t_span, X0, opts);

%% Plotting
% state is in form [r phi dr dphi]
sim_r = X(:,1);
sim_phi = X(:,2);
x_traj = sim_r .* cos(sim_phi);
y_traj = sim_r .* sin(sim_phi);

% Plot the trajectory
plot(x_traj, y_traj, 'r--', 'LineWidth', 1.5);
plot(x, y, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start point

legend('Black Hole', 'Photon Sphere','Light path', 'Light starting point',  ...
       'TextColor', 'w', 'Color', 'k')
title(sprintf('Black Hole Simulation\nPos: (%.2f, %.2f)', x, y), 'Color', 'w');

saveas(gcf, 'Static_4.png')

%% ============= Geodesic Function ================
function dXdt = geodesic(~, X, r_s, noise_scale, noise_accel)
    % input X is the state array
    % break state array into each component (easier to understand)
    r = X(1); 
    phi = X(2); 
    dr_dlambda = X(3); 
    dphi_dlambda = X(4);
    
    % Stop if inside event horizon
    if r <= r_s
        dXdt = zeros(4,1);
        return
    end
    
    % Noise multipliers (centered at 1, small variance)
    % This creates multiplicative noise: value ≈ 1 ± noise_scale
    mult_f = normrnd(1, noise_scale);       % affects metric
    mult_acc3 = normrnd(1, noise_accel);    % affects radial acceleration
    mult_acc4 = normrnd(1, noise_accel);    % affects angular acceleration
    
    % Schwarzschild metric component with multiplicative noise
    f = (1 - r_s / r) * mult_f;
    
    % Conserved energy for null geodesics (light rays)
    E = f * sqrt((dr_dlambda^2)/(f^2) + (r^2 * dphi_dlambda^2)/f);
    dt_dlambda = E / f;
    
    % Initialize the dXdt array which is a 4x1 matrix
    % dX should be [dr dphi drdot dphidot]
    dXdt = zeros(4,1);
    
    dXdt(1) = dr_dlambda;
    dXdt(2) = dphi_dlambda;
    
    % Radial acceleration with multiplicative noise
    base_dr = -(   (r_s/(2*r^2)) * f * (dt_dlambda)^2   ) ...
              +(   (r_s/(2*r^2*f)) * (dr_dlambda^2)      ) ...
              +(   (r - r_s) * (dphi_dlambda^2)          );
    
    dXdt(3) = base_dr * mult_acc3;  % Multiplicative randomness
    
    % Angular acceleration with multiplicative noise
    base_dphi = -(2 * dr_dlambda * dphi_dlambda) / r;
    
    dXdt(4) = base_dphi * mult_acc4;  % Multiplicative randomness
end