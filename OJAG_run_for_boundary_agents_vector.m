%% System Parameters
N = 40; % Number of spatial divisions
num_leaders = 7; % Number of leaders
M = 1000; % Number of time divisions
T = 4; % Total time
a = 0.004; % Diffusion coefficient
kappa = 1; % Boundary coefficient
K = 1.4; % Coupling coefficient
sigma = 6; % Wentzell coefficient
dx = 1 / N;

% Vector-valued initial condition gamma: returns 3 × (N+1)
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];

% Vector-valued nonlinearity (elementwise)
f = @(t, z) 0.3 * z;

%% Solve PDEs
[Error_wentzell, t1, x1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders);
[Error_dirichlet, t11, x11] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, 0, gamma, f, num_leaders);

%% Extract Boundary Errors
idx1 = 1:N+1;                           % Indices for e1
idx2 = N+2:2*(N+1);                     % Indices for e2
idx3 = 2*(N+1)+1:3*(N+1);               % Indices for e3

% Time vector
time = t1;

%% Compute ||e|| and log(||e||) over time for boundaries

% Left boundary (x = 0)
l_w = Error_wentzell([idx1(1), idx2(1), idx3(1)], :);
l_d = Error_dirichlet([idx1(1), idx2(1), idx3(1)], :);
norm_left_w = sqrt(sum(l_w.^2, 1));
norm_left_d = sqrt(sum(l_d.^2, 1));
log_norm_left_w = log(norm_left_w);
log_norm_left_d = log(norm_left_d);

% Right boundary (x = 1)
r_w = Error_wentzell([idx1(end), idx2(end), idx3(end)], :);
r_d = Error_dirichlet([idx1(end), idx2(end), idx3(end)], :);
norm_right_w = sqrt(sum(r_w.^2, 1));
norm_right_d = sqrt(sum(r_d.^2, 1));
log_norm_right_w = log(norm_right_w);
log_norm_right_d = log(norm_right_d);
%% Plot: ||e(t,x)|| vs Time
figure;
subplot(1,2,1);
plot(time, norm_left_w, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5); hold on;
plot(time, norm_left_d, '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\|e(t,0)\|$', 'Interpreter', 'latex');
title('Left Boundary Agent');
legend({'$\sigma \neq 0$', '$\sigma = 0$'}, 'Interpreter', 'latex');
grid on;

subplot(1,2,2);
plot(time, norm_right_w, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5); hold on;
plot(time, norm_right_d, '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\|e(t,1)\|$', 'Interpreter', 'latex');
title('Right Boundary Agent');
legend({'$\sigma \neq 0$', '$\sigma = 0$'}, 'Interpreter', 'latex');
grid on;

%% Plot: log(||e(t,x)||) vs Time
figure;
subplot(1,2,1);
plot(time, log_norm_left_w, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5); hold on;
plot(time, log_norm_left_d, '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\log(\|e(t,0)\|)$', 'Interpreter', 'latex');
title('Left Boundary Agent ');
legend({'$\sigma \neq 0$', '$\sigma = 0$'}, 'Interpreter', 'latex');
grid on;

subplot(1,2,2);
plot(time, log_norm_right_w, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5); hold on;
plot(time, log_norm_right_d, '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\log(\|e(t,1)\|)$', 'Interpreter', 'latex');
title('Right Boundary Agent');
legend({'$\sigma \neq 0$', '$\sigma = 0$'}, 'Interpreter', 'latex');
grid on;


