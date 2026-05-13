%% System Parameters
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
dx = 1 / N; % Spatial step

% Solving the PDEs
[e_L, t1, x1] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders); % Wentzell ( \sigma \neq 0) + communication
[e_C, t2, x2] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders); % wentzell ( \sigma \neq 0) + no communication

% Solving for non-communicating leaders for 14, and 12 leaders respectively
[e_C13, t213, x213] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, sigma, gamma, f, 14); % for 14 leaders
[e_C12, t212, x212] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, sigma, gamma, f, 12);  % for 12 leaders


[e_L1, t11, x11] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, 0, gamma, f, num_leaders); % \sigma =0 + communication
[e_C1, t22, x22] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, 0, gamma, f, num_leaders); % \sigma =0 + no communication
%% Comparision of \|e(t, \cdot)\|^2 for wentzell with and without communication

% Compute L2 norm squared for linear and constant control for Wentzell
L2_norm_L = sqrt(sum(e_L.^2, 1) * dx); % Linear control
L2_norm_C = sqrt(sum(e_C.^2, 1) * dx); % Constant control

% Compute L2 norm squared for linear and constant control for Dirichlet 
L2_norm_L1 = sqrt(sum(e_L1.^2, 1) * dx); % Linear control
L2_norm_C1 = sqrt(sum(e_C1.^2, 1) * dx); % Constant control

L2_norm_C13 = sqrt(sum(e_C13.^2, 1) * dx); % Constant control for 13 leaders
L2_norm_C12 = sqrt(sum(e_C12.^2, 1) * dx); % Constant control for 12 leaders

% Compute log of squared L2 norm for Wentzell
log_L2_norm_L = log(L2_norm_L); % Log of L2 norm squared for linear control
log_L2_norm_C = log(L2_norm_C); % Log of L2 norm squared for constant control

% Compute log of squared L2 norm for Dirichlet 
log_L2_norm_L1 = log(L2_norm_L1); % Log of L2 norm squared for linear control
log_L2_norm_C1 = log(L2_norm_C1); % Log of L2 norm squared for constant control

log_L2_norm_C13 = log(L2_norm_C13);
log_L2_norm_C12 = log(L2_norm_C12);

%% Figure 1: log(\|e(t,\cdot\|)) versus t for Dirichlet and Wentzell with or without communication 
% OR Comparitive analysis of Boudary Control strategies + Cimmunication
% strategy  (This can be used for sub-plot file)
figure;
plot(t1, log_L2_norm_L, 'Color', [1, 0.5, 0], 'LineWidth', 2); 
hold on;
plot(t2, log_L2_norm_C, 'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2); 
hold on;
plot(t11, log_L2_norm_L1, 'Color', 'red', 'LineWidth', 2); 
hold on;
plot(t22, log_L2_norm_C1, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|^2\right)$$', 'Interpreter', 'latex');
title('Comparison of log of Square of L^2 Norm: Dirichlet vs Wentzell');
legend({'Wentzell + Communication', 'Wentzell + No Communication',...
    'Dirichlet + Communication', 'Dirichlet + No Communication'}, ...
    'Location', 'best', 'FontSize', 28, 'NumColumns', 2);
grid on;


%% Figure 2:   \log(\|e(t, \cdot)\|^2)\) vs \(t\) Wentzell with (7 leaders)...
% or without leader (7,12 and 13 leaders) communication.
% Compartive analysis of communication strategy for \sigma \neq 0
  figure;
    plot(t1, log_L2_norm_L, 'Color', [1, 0.5, 0],'LineWidth', 3, 'DisplayName', 'Communicating 7 leaders'); hold on;
    plot(t2, log_L2_norm_C, 'Color', 'blue', 'LineWidth', 3, 'DisplayName', 'Non-communicating 7 leaders');
    hold on;
    plot(t213, log_L2_norm_C13, 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Non-communicating 13 leaders');
    hold on;
    plot(t213, log_L2_norm_C12, 'Color', 'red','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Non-communicating 12 leaders');
   xlabel('$$ t $$', 'Interpreter', 'latex');
    ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
    title('Comparison of Square of L^2 norm: communicating and non communicating with different leaders');
    legend('show');
    grid on;

%% Figure 3: log(e) vs t for different values of sigma = 0.01, 0.5, 1
% For s1 sigma = 0.01

%Effect of sigma on the PDE system
[e_Ls1, t1s1, x1s1] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, 0.01, gamma, f, num_leaders); % Wentzell + communication
[e_Cs1, t2s1, x2s1] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, 0.01, gamma, f, num_leaders); % wentzell + no communication
% For s2 sigma = 0.5
[e_Ls2, t1s2, x1s2] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, 0.5, gamma, f, num_leaders); % Wentzell + communication
[e_Cs2, t2s2, x2s2] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, 0.5, gamma, f, num_leaders); % wentzell + no communication
% For s3 sigma = 1
[e_Ls3, t1s3, x1s3] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, 1, gamma, f, num_leaders); % Wentzell + communication
[e_Cs3, t2s3, x2s3] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, 1, gamma, f, num_leaders); % wentzell + no communication

% Compute L2 norm squared for linear and constant control for Wentzell
L2_norm_Ls1 = sqrt(sum(e_Ls1.^2, 1) * dx); % Linear control
L2_norm_Cs1 = sqrt(sum(e_Cs1.^2, 1) * dx); % Constant control

L2_norm_Ls2 = sqrt(sum(e_Ls2.^2, 1) * dx); % Linear control
L2_norm_Cs2 = sqrt(sum(e_Cs2.^2, 1) * dx); % Constant control

L2_norm_Ls3 = sqrt(sum(e_Ls3.^2, 1) * dx) ; % Linear control
L2_norm_Cs3 = sqrt(sum(e_Cs3.^2, 1) * dx); % Constant control

% Compute log of squared L2 norm for Wentzell
log_L2_norm_Ls1 = log(L2_norm_Ls1); % Log of L2 norm squared for linear control
log_L2_norm_Cs1 = log(L2_norm_Cs1); % Log of L2 norm squared for constant control


log_L2_norm_Ls2 = log(L2_norm_Ls2); % Log of L2 norm squared for linear control
log_L2_norm_Cs2 = log(L2_norm_Cs2); % Log of L2 norm squared for constant control


log_L2_norm_Ls3 = log(L2_norm_Ls3); % Log of L2 norm squared for linear control
log_L2_norm_Cs3 = log(L2_norm_Cs3); % Log of L2 norm squared for constant control

figure;
% Communicating (sigma = 0.01)
plot(t1s1, log_L2_norm_Ls1, 'Color', [1, 0.5, 0], 'LineWidth', 1);  
hold on;
% Communicating (sigma = 0.5)
plot(t1s2, log_L2_norm_Ls2, 'Color', 'red', 'LineWidth', 1); 
hold on;
% Communicating (sigma = 1)
plot(t1s3, log_L2_norm_Ls3, 'Color', 'green', 'LineWidth', 1);
hold on;

% Non-Communicating (sigma = 0.01)
plot(t2s1, log_L2_norm_Cs1, 'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 4); 
hold on;
% Non-Communicating (sigma = 0.5)
plot(t2s2, log_L2_norm_Cs2, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 4); 
hold on;
% Non-Communicating (sigma = 1)
plot(t2s3, log_L2_norm_Cs3, 'Color', 'green', 'LineStyle', '--', 'LineWidth', 4); 

% Axis labels and title
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
title('Effect of $$\sigma$$ on the system', 'Interpreter', 'latex');

% Legend customization
legend({'Communicating, \sigma = 0.01', 'Communicating, \sigma = 0.5', 'Communicating, \sigma = 1',...
    'Non-communicating, \sigma = 0.01', 'Non-communicating, \sigma = 0.5', 'Non-communicating, \sigma = 1'}, ...
    'Location', 'best', 'FontSize', 28, 'NumColumns', 2);
grid on;


%% Fourth Figure

% Indices for boundary agents
left_agent_index = 1; % x=0
right_agent_index = N+1; % x=1

% Square of errors over time
error_left_sigma_nonzero = e_L(left_agent_index, :).^2;
error_left_sigma_zero = e_L1(left_agent_index, :).^2;

error_right_sigma_nonzero = e_L(right_agent_index, :).^2;
error_right_sigma_zero = e_L1(right_agent_index, :).^2;

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

   






