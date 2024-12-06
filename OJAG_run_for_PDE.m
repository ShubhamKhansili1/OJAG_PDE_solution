% Parameters
N = 42; % Number of spatial divisions
num_leaders = 7; % Number of leaders
M = 1000; % Number of time divisions
T = 4; % Total time
a = 0.001; % Diffusion coefficient
kappa = 2; % Coefficient in boundary condition
K = 1.4; % Coefficient in boundary condition
sigma = 0.8; % Coefficient in boundary condition
gamma = @(x) 1.5*sin(x);
f = @(t, z) 0.3*z;

% Solve the PDEs
[e_L, t1, x1] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders);
[e_C, t2, x2] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders);

[e_C13, t213, x213] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, sigma, gamma, f, 14);


dx = 1 / N; % Spatial step


% Compute L2 norm squared for linear and constant control
L2_norm_sq_L = sum(e_L.^2, 1) * dx; % Linear control
L2_norm_sq_C = sum(e_C.^2, 1) * dx; % Constant control

L2_norm_sq_C13 = sum(e_C13.^2, 1) * dx; % Constant control


% Compute the logarithm of the solutions element-wise
log_e_L = log(abs(e_L) + 1e-12); % Adding a small value to avoid log(0)
log_e_C = log(abs(e_C) + 1e-12); % Adding a small value to avoid log(0)

log_e_C13 = log(abs(e_C13) + 1e-12); % Adding a small value to avoid log(0)

% Compute the average of log(e) over space for each time step
log_e_L_avg = mean(log_e_L, 1); % Average over spatial points
log_e_C_avg = mean(log_e_C, 1); % Average over spatial points

log_e_C_avg13 = mean(log_e_C13, 1); % Average over spatial points


% Plot the \|e(t,\cdot)\|^2 versus time.
  figure;
    plot(t1, L2_norm_sq_L, 'Color', [1, 0.5, 0],'LineWidth', 1, 'DisplayName', 'Linear'); hold on;
    plot(t2, L2_norm_sq_C, 'Color', 'blue', 'LineWidth', 1, 'DisplayName', 'Constant');
    hold on;
    plot(t213, L2_norm_sq_C13, 'Color', 'red', 'LineWidth', 1, 'DisplayName', 'Constant13');
    xlabel('Time (t)');
    ylabel('Square of L^2 norm');
    title('Comparison of Square of L^2 norm: Linear vs. Constant');
    legend('show');
    grid on;


% Plot log(e) versus time
figure;
plot(t1, log_e_L_avg, 'Color', [1, 0.5, 0], 'LineWidth', 1, 'DisplayName', 'Log Linear Control'); hold on;
plot(t2, log_e_C_avg, 'Color', 'blue', 'LineWidth', 2, 'DisplayName', 'Log Constant Control'); hold on;
plot(t213, log_e_C_avg13, 'Color', 'red', 'LineWidth', 2, 'DisplayName', 'Log Constant Control_13'); 
xlabel('Time (t)');
ylabel('Average Logarithm of e');
title('Log(e) vs Time: Linear vs Constant Control');
legend('show');
grid on;

% 2D comparison of logarithmic values at the final time step
figure;
plot(x1, log_e_L(:, end), 'LineWidth', 1, 'DisplayName', 'Log Linear Control'); hold on;
plot(x2, log_e_C(:, end), '--', 'LineWidth', 2, 'DisplayName', 'Log Constant Control');
xlabel('Spatial Domain (x)');
ylabel('Logarithmic Value of e(t,x)');
title('Comparison of Logarithmic Values at Final Time Step');
legend('show');
grid on;


% 3D Surface plot for Linear Control
figure;
surf(t1, x1, e_L, 'EdgeColor', 'none');
xlabel('Time (t)');
ylabel('Spatial Domain (x)');
zlabel('Solution e_L(x,t)');
title('3D Surface Plot: Linear Control');
colorbar;
view(3);

% 3D Surface plot for Constant Control
figure;
surf(t2, x2, e_C, 'EdgeColor', 'none');
xlabel('Time (t)');
ylabel('Spatial Domain (x)');
zlabel('Solution e_C(x,t)');
title('3D Surface Plot: Constant Control');
colorbar;
view(3);




