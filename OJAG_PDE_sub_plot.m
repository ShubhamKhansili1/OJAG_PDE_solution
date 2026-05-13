% Parameters
N = 40; % Number of spatial divisions
num_leaders = 7; % Number of leaders
M = 1000; % Number of time divisions
T = 4; % Total time
a = 0.004; % Diffusion coefficient
kappa =1; % Coefficient in boundary condition
K = 1.4; % Coefficient in boundary condition
sigma = 6; % Coefficient in boundary condition
gamma = @(x) 1.5*sin(x);
f = @(t, z) 0.3*z;

% Solve the PDEs
[e_L, t1, x1] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders); % Wentzell + communication
[e_C, t2, x2] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders); % wentzell + no communication

[e_L1, t11, x11] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, 0, gamma, f, num_leaders); % Dirichlet + communication
[e_C1, t22, x22] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, 0, gamma, f, num_leaders); % Dirichlet + no communication

%% For a = 0.001, solving PDEs

[e_LA, t1A, x1A] = OJAG_solveParabolicPDE(N, M, T, 0.001, kappa, K, sigma, gamma, f, num_leaders); % Wentzell + communication
[e_CA, t2A, x2A] = OJAG_solveParabolicPDE_constant(N, M, T, 0.001, kappa, K, sigma, gamma, f, num_leaders); % wentzell + no communication

[e_L1A, t11A, x11A] = OJAG_solveParabolicPDE(N, M, T, 0.001, kappa, K, 0, gamma, f, num_leaders); % Dirichlet + communication
[e_C1A, t22A, x22A] = OJAG_solveParabolicPDE_constant(N, M, T, 0.001, kappa, K, 0, gamma, f, num_leaders); % Dirichlet + no communication

dx = 1 / N; % Spatial step

% Compute L2 norm squared for linear and constant control for Wentzell
L2_norm_L = sqrt(sum(e_L.^2, 1) * dx); % Linear control
L2_norm_C = sqrt(sum(e_C.^2, 1) * dx); % Constant control

% Compute L2 norm squared for linear and constant control for Dirichlet 
L2_norm_L1 = sqrt(sum(e_L1.^2, 1) * dx); % Linear control
L2_norm_C1 = sqrt(sum(e_C1.^2, 1) * dx); % Constant control



% Compute log of squared L2 norm for Wentzell
log_L2_norm_L = log(L2_norm_L); % Log of L2 norm squared for linear control
log_L2_norm_C = log(L2_norm_C); % Log of L2 norm squared for constant control

% Compute log of squared L2 norm for Dirichlet 
log_L2_norm_L1 = log(L2_norm_L1); % Log of L2 norm squared for linear control
log_L2_norm_C1 = log(L2_norm_C1); % Log of L2 norm squared for constant control



%% For a = 0.004

% Compute L2 norm squared for linear and constant control for Wentzell
L2_norm_LA = sqrt(sum(e_LA.^2, 1) * dx); % Linear control
L2_norm_CA = sqrt(sum(e_CA.^2, 1) * dx); % Constant control

% Compute L2 norm squared for linear and constant control for Dirichlet 
L2_norm_L1A = sqrt(sum(e_L1A.^2, 1) * dx); % Linear control
L2_norm_C1A = sqrt(sum(e_C1A.^2, 1) * dx); % Constant control

% Compute log of squared L2 norm for Wentzell
log_L2_norm_LA = log(L2_norm_LA); % Log of L2 norm squared for linear control
log_L2_norm_CA = log(L2_norm_CA); % Log of L2 norm squared for constant control

% Compute log of squared L2 norm for Dirichlet 
log_L2_norm_L1A = log(L2_norm_L1A); % Log of L2 norm squared for linear control
log_L2_norm_C1A = log(L2_norm_C1A); % Log of L2 norm squared for constant control


%% Figure 1: Plot showing log of squared L2 norm versus time with focus on 
% comparing Dirichlet and Wentzell BCs.
figure;

% First subplot (top or left plot)
subplot(2,1,1);     % 2 rows, 1 column, first plot
plot(t1, log_L2_norm_L, 'Color', [1, 0.5, 0], 'LineWidth', 2); 
hold on;
plot(t2, log_L2_norm_C, 'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2); 
hold on;
plot(t11, log_L2_norm_L1, 'Color', 'red', 'LineWidth', 2); 
hold on;
plot(t22, log_L2_norm_C1, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs. $t$: $a = 0.004$','Interpreter','latex');

legend({'$\sigma = 6$ with communication', ...
        '$\sigma = 6$ with no communication', ...
        '$\sigma = 0$ with communication', ...
        '$\sigma = 0$ with no communication'}, ...
        'Location', 'best', 'FontSize', 25, 'NumColumns', 2, 'Interpreter', 'latex');
grid on;

% Second subplot (bottom or right plot)
subplot(2,1,2);     % 2 rows, 1 column, second plot

plot(t1A, log_L2_norm_LA, 'Color', [1, 0.5, 0], 'LineWidth', 2); 
hold on;
plot(t2A, log_L2_norm_CA, 'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2); 
hold on;
plot(t11A, log_L2_norm_L1A, 'Color', 'red', 'LineWidth', 2); 
hold on;
plot(t22A, log_L2_norm_C1A, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs. $t$: $a = 0.001$','Interpreter','latex');

legend({'$\sigma = 6$ with communication', ...
        '$\sigma = 6$ with no communication', ...
        '$\sigma = 0$ with communication', ...
        '$\sigma = 0$ with no communication'}, ...
        'Location', 'best', 'FontSize', 25, 'NumColumns', 2, 'Interpreter', 'latex');

grid on;








