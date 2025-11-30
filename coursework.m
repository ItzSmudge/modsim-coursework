%% Modelling a black hole and light rays in matlab (this is the coursework file we need to use) 
clc; clear;

% assuming G=c=1 
%geometrized units (everything has the same unit)
r_s = 5; %means 5 units of length instead of 5m or 5Km

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
fill(r_s*cos(theta), r_s*sin(theta), 'k',  'EdgeColor', 'c',"LineWidth",1.5);
% Photon Sphere (r = 1.5 * rs) 
% light is most unstable in the photon sphere
plot(1.5*r_s*cos(theta), 1.5*r_s*sin(theta), 'y--', 'LineWidth', 1.5);

%% ============= The hard stuff ================
% do we set initial polar conditions 
% or do we set initial cartesian condition 

%polar means we just do everything and convert to x/y at the end 
% r,phi,dr,dphi DO need cartesian position and directions in their
% equations
% essentially means we need to make initial conditions in cartesian 
% then convert them to polar coordinates 

x = -30; %horizontal position
y= 13;   %vertical position 
vx = 185;  % horizontal velocity (honestly dont know what changing this does but it does)
vy = -1;  % vertical velocity 

% convert into polar coordinates
r0 = sqrt(x^2 + y^2);
phi0 = atan2(y, x);
dr0 = vx * cos(phi0) + vy * sin(phi0);
dphi0 = (-vx * sin(phi0) + vy * cos(phi0)) / r0;
%initial state array 
X0 = [r0;phi0;dr0;dphi0];

%% create the geodesic function 
function dXdt = geodesic(t,X,r_s)
    % input X is the state array
    % break state array into each component (easier to understand)
    r = X(1); phi = X(2); dr_dlambda = X(3); dphi_dlambda = X(4);
    noise_scale = 1e-2;
    noise_accel = 1e-6;

    % Stop if inside event horizon
    if r <= r_s
        dXdt = zeros(4,1);
        return;
    end
    
    % Schwarzschild metric component
    f = (1 - r_s / r) * ( 1  + (noise_scale * randn));  % ---- Gaussian noise added here ---- %
    
    % Conserved energy for null geodesics (light rays)
    % ---- we could also add noise to these 2 equations potentially ----%
    % ---- Maybe not the dt_dlambda eq though ----%
    E = f * sqrt((dr_dlambda^2)/(f^2) + (r^2 * dphi_dlambda^2)/f);
    dt_dlambda = E / f;

    %initialise the dXdt array which is a 4x1 matrix 
    % dX should be [dr dphi drdot dphidot]
    dXdt = zeros(4,1);
    dXdt(1) = dr_dlambda;
    dXdt(2) = dphi_dlambda;
    dXdt(3) = -(   (r_s/(2*r^2)) * f * (dt_dlambda)^2   )...
              +(  (r_s/(2*r^2*f))  * (dr_dlambda)^2     )...
              +(  (r-r_s) * (dphi_dlambda)^2           )...
               +  (noise_accel*randn); 
    dXdt(4) = ( noise_accel*randn) + (-2 * dr_dlambda * dphi_dlambda)/r; %regular noise added here
    %returned array dXdt essentially return the next state 
end 

%% simulation 
t_span = [0 1000]; %longer simulation means more of the ray of light is drawn
geodesic_equation = @(t,X) geodesic(t,X,r_s);
[sim_time, X] = ode45(geodesic_equation,t_span,X0);

%% plotting 
% state is in form [r phi dr dphi]
sim_r = X(:,1);
sim_phi = X(:,2);
x_traj = sim_r .* cos(sim_phi);
y_traj = sim_r .* sin(sim_phi);

% Plot the trajectory
plot(x_traj, y_traj, 'r--', 'LineWidth', 1);
plot(x, y, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start point
legend("Black Hole","Photon Sphere","Light starting point","light path",'TextColor', 'w', 'Color', 'k')
title(sprintf('Black Hole Simulation\nPos: (%.2e, %.2e) | Vel: (%.2e, %.2e)', x, y, vx, vy), 'Color', 'w');

% what have i learnt from this 
% realistically a whole load of nothing