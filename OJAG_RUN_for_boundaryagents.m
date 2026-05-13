%% System Parameters
N = 40; % Number of spatial divisions
num_leaders = 7; % Number of leaders
M = 1000; % Number of time divisions
T = 4; % Total time
a = 0.004; % Diffusion coefficient
kappa =1; % Coefficient in boundary condition
K = 1.4; % Coefficient in boundary condition
sigma = 6; % Coefficient in boundary condition
gamma = @(x) 1.5*sin(x) ;
f = @(t, z) 0.3*z;
dx = 1 / N; % Spatial step

% Solving the PDEs
[Error_wentzell, t1, x1] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders); % Wentzell ( \sigma \neq 0) + communication

[Error_dirichlet, t11, x11] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, 0, gamma, f, num_leaders); % \sigma =0 + communication

% Indices for boundary agents
left_agent_index = 1; % x=0
right_agent_index = N+1; % x=1

% Square of errors over time
error_left_sigma_nonzero = Error_wentzell(left_agent_index, :).^2;
error_left_sigma_zero = Error_dirichlet(left_agent_index, :).^2;

error_right_sigma_nonzero = Error_wentzell(right_agent_index, :).^2;
error_right_sigma_zero = Error_dirichlet(right_agent_index, :).^2;

% Create figure
figure;

% Left subplot: Left boundary agent
subplot(1,2,1);
plot(t1, error_left_sigma_nonzero, 'LineWidth', 1.5);
hold on;
plot(t1, error_left_sigma_zero, '--', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Squared Error');
title('Left Boundary Agent');
legend('\sigma \neq 0', '\sigma = 0');
grid on;

% Right subplot: Right boundary agent
subplot(1,2,2);
plot(t1, error_right_sigma_nonzero, 'LineWidth', 1.5);
hold on;
plot(t1, error_right_sigma_zero, '--', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Squared Error');
title('Right Boundary Agent');
legend('\sigma \neq 0', '\sigma = 0');
grid on;