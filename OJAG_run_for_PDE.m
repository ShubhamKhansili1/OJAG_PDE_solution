N = 42; % Number of spatial divisions
num_leaders =7; %Number of leaders
M = 1000; % Number of time divisions
T = 4; % Total time
a = 0.0001; % Diffusion coefficient
kappa = 2; % Coefficient in boundary condition
K = 1.4; % Coefficient in boundary condition
sigma = 0.8; % Coefficient in boundary condition
gamma = @(x) 1.5*sin(x);
f=@(t,z) 0.3*sin(z);


% Solve the PDE
[e_L, t1, x1] = OJAG_solveParabolicPDE(N, M, T, a, kappa, K, sigma, gamma, f,num_leaders);
[e_C, t2, x2] = OJAG_solveParabolicPDE_constant(N, M, T, a, kappa, K, sigma, gamma, f,num_leaders);
dx = 1 / N; % Spatial step

L2_norm_sq_L = sum(e_L.^2, 1) * dx; % Sum across spatial points, scale by dx

L2_norm_sq_C = sum(e_C.^2, 1) * dx; % Sum across spatial points, scale by dx

 % Plot the comparison
    figure;
    plot(t1, L2_norm_sq_L, 'LineWidth', 1, 'DisplayName', 'Linear'); hold on;
    plot(t2, L2_norm_sq_C, '--', 'LineWidth', 2, 'DisplayName', 'Constant');
    xlabel('Time (t)');
    ylabel('Square of L^2 norm');
    title('Comparison of Square of L^2 norm: Linear vs. Constant');
    legend('show');
    grid on;




