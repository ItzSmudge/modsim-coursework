%% Monte Carlo Simulation - Black Hole Escape Analysis
clc; clear;

% Black hole parameters
r_s = 5; % Schwarzschild radius

%% Noise parameters
noise_scale = 1e-2;  % affects metric factor
noise_accel = 1e-1;  % affects accelerations

%% Monte Carlo parameters
n_simulations = 1000; % number of trials with different noise realizations
t_span = [0 1000];    % simulation time
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4);

%% Fixed starting position and velocity
x_start = -30;  % horizontal position
y_start = 13.3;   % vertical position
vx_start = 185; % horizontal velocity
vy_start = -1;  % vertical velocity

%% Initialize result arrays
escaped = zeros(n_simulations, 1);  % 1 = escaped, 0 = captured
final_radii = zeros(n_simulations, 1);

%% Convert fixed starting position to polar coordinates
r0 = sqrt(x_start^2 + y_start^2);
phi0 = atan2(y_start, x_start);
dr0 = vx_start * cos(phi0) + vy_start * sin(phi0);
dphi0 = (-vx_start * sin(phi0) + vy_start * cos(phi0)) / r0;
X0 = [r0; phi0; dr0; dphi0];

%% Monte Carlo Loop
fprintf('Running Monte Carlo simulation with %d trials...\n', n_simulations);
fprintf('Starting position: (%.2f, %.2f)\n', x_start, y_start);
fprintf('Starting velocity: (%.2f, %.2f)\n\n', vx_start, vy_start);
tic;

for i = 1:n_simulations
    % Run simulation with same initial conditions but different noise
    geodesic_equation = @(t,X) geodesic(t, X, r_s, noise_scale, noise_accel);
    [~, X] = ode45(geodesic_equation, t_span, X0, opts);
    
    % Check final radius
    final_r = X(end, 1);
    final_radii(i) = final_r;
    
    % Determine if escaped (1) or captured (0)
    if final_r <= r_s * 1.5
        escaped(i) = 0; % Captured by black hole
    else
        escaped(i) = 1; % Escaped
    end
    
    % Progress indicator
    if mod(i, 100) == 0
        fprintf('Completed %d/%d simulations (%.1f%%)\n', i, n_simulations, 100*i/n_simulations);
    end
end

elapsed_time = toc;

%% Calculate and display results
n_escaped = sum(escaped);
n_captured = n_simulations - n_escaped;
escape_rate = 100 * n_escaped / n_simulations;

fprintf('\n========== RESULTS ==========\n');
fprintf('Total simulations: %d\n', n_simulations);
fprintf('Escaped: %d (%.2f%%)\n', n_escaped, escape_rate);
fprintf('Captured: %d (%.2f%%)\n', n_captured, 100 - escape_rate);
fprintf('Simulation time: %.2f seconds\n', elapsed_time);
fprintf('============================\n\n');

%% Additional statistics
fprintf('Starting position: (%.2f, %.2f)\n', x_start, y_start);
fprintf('Starting velocity: (%.2f, %.2f)\n', vx_start, vy_start);
fprintf('Starting radius: %.2f\n', r0);
fprintf('Schwarzschild radius: %.2f\n', r_s);
fprintf('Photon sphere radius: %.2f\n', 1.5*r_s);

% Statistics on final radii
fprintf('\nFinal radii statistics:\n');
fprintf('  Mean: %.2f\n', mean(final_radii));
fprintf('  Median: %.2f\n', median(final_radii));
fprintf('  Min: %.2f\n', min(final_radii));
fprintf('  Max: %.2f\n', max(final_radii));

%% Save results to file
results = struct();
results.n_simulations = n_simulations;
results.n_escaped = n_escaped;
results.n_captured = n_captured;
results.escape_rate = escape_rate;
results.starting_position = [x_start, y_start, vx_start, vy_start];
results.escaped = escaped;
results.final_radii = final_radii;
results.parameters.r_s = r_s;
results.parameters.r0 = r0;
results.parameters.noise_scale = noise_scale;
results.parameters.noise_accel = noise_accel;

save('monte_carlo_results.mat', 'results');
fprintf('\nResults saved to monte_carlo_results.mat\n');

%% ============= Geodesic Function ================
function dXdt = geodesic(~, X, r_s, noise_scale, noise_accel)
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
    mult_f = normrnd(1, noise_scale);       % affects metric
    mult_acc3 = normrnd(1, noise_accel);    % affects radial acceleration
    mult_acc4 = normrnd(1, noise_accel);    % affects angular acceleration
    
    % Schwarzschild metric component with multiplicative noise
    f = (1 - r_s / r) * mult_f;
    
    % Conserved energy for null geodesics (light rays)
    E = f * sqrt((dr_dlambda^2)/(f^2) + (r^2 * dphi_dlambda^2)/f);
    dt_dlambda = E / f;
    
    % Initialize the dXdt array
    dXdt = zeros(4,1);
    
    dXdt(1) = dr_dlambda;
    dXdt(2) = dphi_dlambda;
    
    % Radial acceleration with multiplicative noise
    base_dr = -(   (r_s/(2*r^2)) * f * (dt_dlambda)^2   ) ...
              +(   (r_s/(2*r^2*f)) * (dr_dlambda^2)      ) ...
              +(   (r - r_s) * (dphi_dlambda^2)          );
    
    dXdt(3) = base_dr * mult_acc3;
    
    % Angular acceleration with multiplicative noise
    base_dphi = -(2 * dr_dlambda * dphi_dlambda) / r;
    
    dXdt(4) = base_dphi * mult_acc4;
end