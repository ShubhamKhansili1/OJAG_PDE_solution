% Parameters
N = 42; % Number of spatial divisions
num_leaders = 7; % Number of leaders
M = 1000; % Number of time divisions
T = 4; % Total time
a = 0.001; % Diffusion coefficient
kappa = 2; % Coefficient in boundary condition
K = 1.4; % Coefficient in boundary condition
sigma = 6; % Coefficient in boundary condition
gamma = @(x) 1.5*sin(x);
f = @(t, z) 0.3*z;

% Solve the PDEs
[e_L, t1, x1] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders); % Wentzell + communication
[e_C, t2, x2] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders); % wentzell + no communication

[e_C13, t213, x213] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, sigma, gamma, f, 13);


[e_L1, t11, x11] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, 0, gamma, f, num_leaders); % Dirichlet + communication
[e_C1, t22, x22] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, 0, gamma, f, num_leaders); % Dirichlet + no communication


dx = 1 / N; % Spatial step

%% Comparision of \|e(t, \cdot)\|^2 for wentzell with and without communication

% Compute L2 norm squared for linear and constant control for Wentzell
L2_norm_sq_L = sum(e_L.^2, 1) * dx; % Linear control
L2_norm_sq_C = sum(e_C.^2, 1) * dx; % Constant control

% Compute L2 norm squared for linear and constant control for Dirichlet 
L2_norm_sq_L1 = sum(e_L1.^2, 1) * dx; % Linear control
L2_norm_sq_C1 = sum(e_C1.^2, 1) * dx; % Constant control

L2_norm_sq_C13 = sum(e_C13.^2, 1) * dx; % Constant control for 13 leaders

% Compute log of squared L2 norm for Wentzell
log_L2_norm_sq_L = log(L2_norm_sq_L); % Log of L2 norm squared for linear control
log_L2_norm_sq_C = log(L2_norm_sq_C); % Log of L2 norm squared for constant control

% Compute log of squared L2 norm for Dirichlet 
log_L2_norm_sq_L1 = log(L2_norm_sq_L1); % Log of L2 norm squared for linear control
log_L2_norm_sq_C1 = log(L2_norm_sq_C1); % Log of L2 norm squared for constant control

log_L2_norm_sq_C13 = log(L2_norm_sq_C13);

%% Figure 1: Plot showing log of squared L2 norm versus time with focus on 
% comparing Dirichlet and Wentzell BCs.
figure;
plot(t1, log_L2_norm_sq_L, 'Color', [1, 0.5, 0], 'LineWidth', 2); 
hold on;
plot(t2, log_L2_norm_sq_C, 'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2); 
hold on;
plot(t11, log_L2_norm_sq_L1, 'Color', 'red', 'LineWidth', 2); 
hold on;
plot(t22, log_L2_norm_sq_C1, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|^2\right)$$', 'Interpreter', 'latex');
title('Comparison of log of Square of L^2 Norm: Dirichlet vs Wentzell');
legend({'Wentzell + Communication', 'Wentzell + No Communication',...
    'Dirichlet + Communication', 'Dirichlet + No Communication'}, ...
    'Location', 'best', 'FontSize', 28, 'NumColumns', 2);
grid on;


%% Figure 2:  Plot showing log(\|e(t,\cdot)\|^2) versus time focusing only
% on the wentzell boundary with and without communication but with
% different number of leaders
  figure;
    plot(t1, log_L2_norm_sq_L, 'Color', [1, 0.5, 0],'LineWidth', 3, 'DisplayName', 'Communicating 7 leaders'); hold on;
    plot(t2, log_L2_norm_sq_C, 'Color', 'blue', 'LineWidth', 3, 'DisplayName', 'Non-communicating 7 leaders');
    hold on;
    plot(t213, log_L2_norm_sq_C13, 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Non-communicating 13 leaders');
   

    xlabel('$$ t $$', 'Interpreter', 'latex');
    ylabel('$$\log\left(\|e(t, \cdot)\|^2\right)$$', 'Interpreter', 'latex');
    title('Comparison of Square of L^2 norm: communicating and non communicating with different leaders');
    legend('show');
    grid on;

 %% Effect of \sigma
% For s1 sigma = 0.01
[e_Ls1, t1s1, x1s1] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, 0.01, gamma, f, num_leaders); % Wentzell + communication
[e_Cs1, t2s1, x2s1] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, 0.01, gamma, f, num_leaders); % wentzell + no communication
% For s2 sigma = 4
[e_Ls2, t1s2, x1s2] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, 0.5, gamma, f, num_leaders); % Wentzell + communication
[e_Cs2, t2s2, x2s2] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, 0.5, gamma, f, num_leaders); % wentzell + no communication
% For s3 sigma = 8
[e_Ls3, t1s3, x1s3] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, 8, gamma, f, num_leaders); % Wentzell + communication
[e_Cs3, t2s3, x2s3] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, 8, gamma, f, num_leaders); % wentzell + no communication

% Compute L2 norm squared for linear and constant control for Wentzell
L2_norm_sq_Ls1 = sum(e_Ls1.^2, 1) * dx; % Linear control
L2_norm_sq_Cs1 = sum(e_Cs1.^2, 1) * dx; % Constant control

L2_norm_sq_Ls2 = sum(e_Ls2.^2, 1) * dx; % Linear control
L2_norm_sq_Cs2 = sum(e_Cs2.^2, 1) * dx; % Constant control

L2_norm_sq_Ls3 = sum(e_Ls3.^2, 1) * dx; % Linear control
L2_norm_sq_Cs3 = sum(e_Cs3.^2, 1) * dx; % Constant control

% Compute log of squared L2 norm for Wentzell
log_L2_norm_sq_Ls1 = log(L2_norm_sq_Ls1); % Log of L2 norm squared for linear control
log_L2_norm_sq_Cs1 = log(L2_norm_sq_Cs1); % Log of L2 norm squared for constant control


log_L2_norm_sq_Ls2 = log(L2_norm_sq_Ls2); % Log of L2 norm squared for linear control
log_L2_norm_sq_Cs2 = log(L2_norm_sq_Cs2); % Log of L2 norm squared for constant control


log_L2_norm_sq_Ls3 = log(L2_norm_sq_Ls3); % Log of L2 norm squared for linear control
log_L2_norm_sq_Cs3 = log(L2_norm_sq_Cs3); % Log of L2 norm squared for constant control

%% FIGURE 3: System behavior across various sigma values
figure;
% Communicating (sigma = 0.01)
plot(t1s1, log_L2_norm_sq_Ls1, 'Color', [1, 0.5, 0], 'LineWidth', 1);  
hold on;
% Communicating (sigma = 4)
plot(t1s2, log_L2_norm_sq_Ls2, 'Color', 'red', 'LineWidth', 1); 
hold on;
% Communicating (sigma = 8)
plot(t1s3, log_L2_norm_sq_Ls3, 'Color', 'green', 'LineWidth', 1);
hold on;

% Non-Communicating (sigma = 0.01)
plot(t2s1, log_L2_norm_sq_Cs1, 'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 4); 
hold on;
% Non-Communicating (sigma = 4)
plot(t2s2, log_L2_norm_sq_Cs2, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 4); 
hold on;
% Non-Communicating (sigma = 8)
plot(t2s3, log_L2_norm_sq_Cs3, 'Color', 'green', 'LineStyle', '--', 'LineWidth', 4); 

% Axis labels and title
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|^2\right)$$', 'Interpreter', 'latex');
title('Effect of $$\sigma$$ on the system', 'Interpreter', 'latex');

% Legend customization
legend({'Communicating, \sigma = 0.01', 'Communicating, \sigma = 0.5', 'Communicating, \sigma = 8',...
    'Non-communicating, \sigma = 0.01', 'Non-communicating, \sigma = 0.5', 'Non-communicating, \sigma = 8'}, ...
    'Location', 'best', 'FontSize', 28, 'NumColumns', 2);
grid on;

   
%% 


% Compute the logarithm of the solutions element-wise for Wentzell
log_e_L = log(abs(e_L) + 1e-12); % Adding a small value to avoid log(0)
log_e_C = log(abs(e_C) + 1e-12); % Adding a small value to avoid log(0)

log_e_C13 = log(abs(e_C13) + 1e-12); % Adding a small value to avoid log(0)

% Compute the logarithm of the solutions element-wise for Dirichlet
log_e_L1 = log(abs(e_L1) + 1e-12); % Adding a small value to avoid log(0)
log_e_C1 = log(abs(e_C1) + 1e-12); % Adding a small value to avoid log(0)

% Compute the average of log(e) over space for each time step
log_e_L_avg = mean(log_e_L, 1); % Average over spatial points
log_e_C_avg = mean(log_e_C, 1); % Average over spatial points

% Compute the average of log(e) over space for each time step
log_e_L1_avg = mean(log_e_L1, 1); % Average over spatial points
log_e_C1_avg = mean(log_e_C1, 1); % Average over spatial points

log_e_C_avg13 = mean(log_e_C13, 1); % Average over spatial points



%% Figure 3 Plot log(e) versus time
figure;
plot(t1, log_e_L_avg, 'Color', [1, 0.5, 0], 'LineWidth', 1, 'DisplayName', 'Communicating 7 leaders'); hold on;
plot(t2, log_e_C_avg, 'Color', 'blue', 'LineWidth', 2, 'DisplayName', 'Not communicating 7 leaders'); hold on;
plot(t213, log_e_C_avg13, 'Color', 'red', 'LineWidth', 2, 'DisplayName', 'Not communicating 13 leaders'); 
xlabel('Time (t)');
ylabel('Log(e)');
title('Log(e) vs Time: communicating vs not communicating leaders');
legend('show');
grid on;


%% Figure 4 2D comparison of logarithmic values at the final time step
figure;
plot(x1, log_e_L(:, end), 'LineWidth', 1, 'DisplayName', 'Log Linear Control'); hold on;
plot(x2, log_e_C(:, end), '--', 'LineWidth', 2, 'DisplayName', 'Log Constant Control');
xlabel('Spatial Domain (x)');
ylabel('Logarithmic Value of e(t,x)');
title('Comparison of Logarithmic Values at Final Time Step');
legend('show');
grid on;


%% Figure 5 3D Surface plot for Linear Control
figure;
surf(t1, x1, e_L, 'EdgeColor', 'none');

% Set the colormap to transition from red to orange and slightly to yellow
nColors = 256; % Number of colors
redChannel = ones(nColors, 1); % Full red
greenChannel = linspace(0, 0.7, nColors)'; % Gradually increase green towards yellow
blueChannel = zeros(nColors, 1); % No blue
colormap([redChannel, greenChannel, blueChannel]);

xlabel('$t$', 'Interpreter', 'latex');
ylabel('$x$', 'Interpreter', 'latex');
zlabel('$e(t,x)$', 'Interpreter', 'latex');
title('3D Surface Plot: communicating leaders');
colorbar;
view(3);

%% Figure 6 3D Surface plot for Constant Control
figure;
surf(t2, x2, e_C, 'EdgeColor', 'none');
% Set the colormap to transition from red to orange and slightly to yellow
nColors = 256; % Number of colors
redChannel = ones(nColors, 1); % Full red
greenChannel = linspace(0, 0.7, nColors)'; % Gradually increase green towards yellow
blueChannel = zeros(nColors, 1); % No blue
colormap([redChannel, greenChannel, blueChannel]);




xlabel('$t$', 'Interpreter', 'latex');
ylabel('$x$', 'Interpreter', 'latex');
zlabel('$e(t,x)$', 'Interpreter', 'latex');
title('3D Surface Plot: Non-communicating leaders');
colorbar;
view(3);


%% Figure 7:  Plot showing log(\|e(t,\cdot)\|^2) versus time focusing only
% on the wentzell boundary with and without communication but with
% different number of leaders
  figure;
    plot(t1, L2_norm_sq_L, 'Color', [1, 0.5, 0],'LineWidth', 1, 'DisplayName', 'Communicating 7 leaders'); hold on;
    plot(t2, L2_norm_sq_C, 'Color', 'blue', 'LineWidth', 1, 'DisplayName', 'Non-communicating 7 leaders');
    hold on;
    plot(t213, L2_norm_sq_C13, 'Color', 'red', 'LineWidth', 1, 'DisplayName', 'Non-communicating 13 leaders');
    xlabel('$$ t $$', 'Interpreter', 'latex');
    ylabel('$$\log\left(\|e(t, \cdot)\|^2\right)$$', 'Interpreter', 'latex');
    title('Comparison of Square of L^2 norm: communicating and non communicating with different leaders');
    legend('show');
    grid on;





