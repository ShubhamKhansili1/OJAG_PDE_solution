
N = 140; num_leaders = 7; M = 10000; T = 4; kappa = 1; K = 1.4; 
sigma= 6; dx = 1 / N; h = 1/N; a = 0.004; delta =1 ; f = @(t, z) 0.3* z; 
% --- Target curve -------------------------------------------
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];   
gamma_in = @(x)[ 0*ones(size(x));  0*ones(size(x));  0*ones(size(x))];

%Solving the error PDE-system
[t14_1, e14_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
[t14_2, e14_2] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma-4, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
[t14_3, e14_3] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma+5, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
[t14_4, e14_4] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma+10, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
[t14_5, e14_5] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + Comm. leader
[t14_6, e14_6] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, 80, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + Comm. leader

%Computing the Log of L^2-norm for case 1
L2_norm14_1 = compute_vec_L2(e14_1,dx);  log_L2_norm14_1  = log(L2_norm14_1);
L2_norm14_2 = compute_vec_L2(e14_2,dx);  log_L2_norm14_2  = log(L2_norm14_2);
L2_norm14_3 = compute_vec_L2(e14_3,dx);  log_L2_norm14_3  = log(L2_norm14_3);
L2_norm14_4 = compute_vec_L2(e14_4,dx);  log_L2_norm14_4  = log(L2_norm14_4);
L2_norm14_5 = compute_vec_L2(e14_5,dx);  log_L2_norm14_5  = log(L2_norm14_5);
L2_norm14_6 = compute_vec_L2(e14_6,dx);  log_L2_norm14_6  = log(L2_norm14_6);

% Define custom colors
colors = [
    1, 0.5, 0;               % Orange (σ = 6)
    0.5, 0, 0;               % Maroon (σ - 1)
    1, 0, 0;                 % Red (σ + 1)
   0.8, 0.6, 0;           % Dark Yellow (σ + 20)
   0.3010, 0.7450, 0.9330  % Cyan → σ = 0
];

% Labels for legend
labels = {'$\sigma = 6$', '$\sigma=2$', '$\sigma = 11$', '$\sigma = 16$', '$\sigma = 0$'};

% Plot both subplots
figure('Name', 'Comparison of Log L2 Norms', 'NumberTitle', 'off');

% ------------------
% Subplot 1: All 5 curves
% ------------------
% ------------------
% Subplot 1
% ------------------
subplot(1,2,1);
hold on;
plot(t14_1, log_L2_norm14_1, 'Color', colors(1,:), 'LineWidth', 3);
plot(t14_2, log_L2_norm14_2, 'Color', colors(2,:), 'LineWidth', 3);
plot(t14_3, log_L2_norm14_3, 'Color', colors(3,:), 'LineWidth', 3);
plot(t14_4, log_L2_norm14_4, 'Color', colors(4,:), 'LineWidth', 3);
plot(t14_5, log_L2_norm14_5, 'Color', colors(5,:), 'LineWidth', 3);

xlabel('Time $t$', 'Interpreter', 'latex', 'FontSize', 33);
ylabel('log $\|e(t,\cdot)\|$', 'Interpreter', 'latex', 'FontSize', 33);

legend(labels, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 25);

grid on;
box on;
set(gca, 'FontSize', 25);  % axes ticks font size

% ------------------
% Subplot 2
% ------------------
subplot(1,2,2);
plot(t14_6, log_L2_norm14_6, 'Color', colors(2,:), 'LineWidth', 3);

xlabel('Time $t$', 'Interpreter', 'latex', 'FontSize', 33);
ylabel('log $\|e(t,\cdot)\|$', 'Interpreter', 'latex', 'FontSize', 33);

legend('$\sigma= 80$', 'Interpreter', 'latex', 'FontSize', 25, 'Location', 'northeast');

grid on;
box on;
set(gca, 'FontSize', 25);  % axes ticks font size
