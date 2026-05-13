% Code for Chapter 2, to show two leaders are not sufficient for Stabilit
% when interaction strength is low. 

%% System Parameter's
N = 70; num_leaders = 0; M = 100000; T = 5; kappa = 1; K = 0; 
sigma= 0.5; dx = 1 / N; h = 1/N; delta =0 ; f = @(t, z) 0.3* z; 
% --- Target curve -------------------------------------------
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];   
gamma_in = @(x)[ 0*ones(size(x));  0*ones(size(x));  0*ones(size(x))];
% For dirichlet conditions, this is the initial position
gamma_in_dirichlet = @(x) gamma(x) .* ((x == 0) | (x == 1));


% ---- Parameters for a = 1 ----
a1 = 1;
[t16_1, e16_1] = OJAG_solveParabolicPDE_vec(N, M, T, a1, kappa, K, sigma, gamma, gamma_in, f, num_leaders, delta);      % Wentzell
[t16_2, e16_2] = OJAG_solveParabolicPDE_vec(N, M, T, a1, kappa, K, 0, gamma, gamma_in_dirichlet, f, num_leaders, delta); % Dirichlet
[t16_3, e16_3] = OJAG_solveParabolicPDE_vec(N, M, T, a1, kappa, K, 0, gamma, gamma_in, f, num_leaders, delta);          % Dynamic Dirichlet

L2_norm16_1 = compute_vec_L2(e16_1, dx);
L2_norm16_2 = compute_vec_L2(e16_2, dx);
L2_norm16_3 = compute_vec_L2(e16_3, dx);

% ---- Parameters for a = 0.1 ----
a2 = 0.01;
[t17_1, e17_1] = OJAG_solveParabolicPDE_vec(N, M, T, a2, kappa, K, sigma, gamma, gamma_in, f, num_leaders, delta);      % Wentzell
[t17_2, e17_2] = OJAG_solveParabolicPDE_vec(N, M, T, a2, kappa, K, 0, gamma, gamma_in_dirichlet, f, num_leaders, delta); % Dirichlet
[t17_3, e17_3] = OJAG_solveParabolicPDE_vec(N, M, T, a2, kappa, K, 0, gamma, gamma_in, f, num_leaders, delta);          % Dynamic Dirichlet

L2_norm17_1 = compute_vec_L2(e17_1, dx);
L2_norm17_2 = compute_vec_L2(e17_2, dx);
L2_norm17_3 = compute_vec_L2(e17_3, dx);

figure('Name', 'Comparison for a = 1 and a = 0.1', 'NumberTitle', 'off');

light_maroon = [0.7, 0.3, 0.3];

% Subplot 1: a = 1
subplot(1,2,1);
hold on;
h1 = plot(t16_1, L2_norm16_1, '-',  'Color', light_maroon, 'LineWidth', 2);  % Wentzell
h2 = plot(t16_2, L2_norm16_2, '--', 'Color', light_maroon, 'LineWidth', 2);  % Dirichlet
h3 = plot(t16_3, L2_norm16_3, ':',  'Color', light_maroon, 'LineWidth', 2);  % Dynamic Dirichlet
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\|e(t,\cdot)\|$$', 'Interpreter', 'latex');
title('$a = 1$', 'Interpreter', 'latex');
grid on;

% Subplot 2: a = 0.1
subplot(1,2,2);
hold on;
plot(t17_1, L2_norm17_1, '-',  'Color', light_maroon, 'LineWidth', 2);  % Wentzell
plot(t17_2, L2_norm17_2, '--', 'Color', light_maroon, 'LineWidth', 2);  % Dirichlet
plot(t17_3, L2_norm17_3, ':',  'Color', light_maroon, 'LineWidth', 2);  % Dynamic Dirichlet
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\|e(t,\cdot)\|$$', 'Interpreter', 'latex');
title('$a = 0.01$', 'Interpreter', 'latex');
grid on;

% Add ONE shared legend below both subplots
% Use invisible axes to position it nicely
hl = legend([h1, h2, h3], ...
    {'Wentzell Boundary Control', ...
     'Dirichlet Boundary Control', ...
     'Dynamic Dirichlet Boundary Control'}, ...
     'Orientation', 'horizontal', ...
     'Interpreter', 'latex', ...
     'Location', 'southoutside');
