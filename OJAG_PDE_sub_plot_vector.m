%% Parameters
N = 40; % Number of spatial divisions
num_leaders = 14; %8
M = 10000;
T = 4;
kappa = 1; kappa1= 1; % 0.5
K = 1.4; K1= 1.4; %8
sigma= 6; sigma1= 6; %8
dx = 1 / N;
h = 1/N;
a = 0.001; a1 = 0.004;
delta =1 ; 

% --- Target curve -------------------------------------------
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];   
%gamma = @(x) [1.5*sin(x); 1.5*sin(x);1.5*sin(x)];  

%gamma_in = @(x) gamma(x) .* (x == 0); 
%gamma_in = @(x) gamma(x) .* ((x == 0) | (x == 1));
%gamma_in = @(x) repmat(gamma(0), 1, numel(x));  % All agents at gamma(0)
gamma_in = @(x)[ 0*ones(size(x));  0*ones(size(x));  0*ones(size(x))];
%gamma = @(x) [h*cos(pi*N*x/3); h*sin(pi*N*x/3); 0*ones(size(x))]; 
% Vector-valued nonlinearity (element-wise)
f = @(t, z) 0.3* z; f1 = @(t, z) 0.3* z; % 3

%% Solve the PDEs for K = 1.4

[  t1, e_L  ] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in,f,num_leaders,delta); % sigma≠0, with comm
[  t2, e_C ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, sigma, gamma,gamma_in, f, num_leaders,delta); % sigma≠0, no comm
[ t11, e_L1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);     % sigma=0, with comm
[ t22, e_C1] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);     % sigma=0, no comm

% Solve the PDEs for K = 3
[  t1_K3, e_L_K3  ] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 3, sigma, gamma, gamma_in,f,14,delta); % sigma≠0, with comm
[  t2_K3, e_C_K3 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 3, sigma, gamma,gamma_in, f, 14,delta); % sigma≠0, no comm
[ t11_K3, e_L1_K3] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 3, 0, gamma,gamma_in, f, 14,delta);     % sigma=0, with comm
[ t22_K3, e_C1_K3] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 3, 0, gamma,gamma_in, f, 14,delta);     % sigma=0, no

% Solve the PDEs for K = 5

[  t1_K5, e_L_K5  ] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 5, sigma, gamma, gamma_in,f,16,delta); % sigma≠0, with comm
[  t2_K5, e_C_K5 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, sigma, gamma,gamma_in, f, 16,delta); % sigma≠0, no comm
[ t11_K5, e_L1_K5] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 5, 0, gamma,gamma_in, f,16,delta);     % sigma=0, with comm
[ t22_K5, e_C1_K5] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, 0, gamma,gamma_in, f, 16,delta);     % sigma=0, no

% Solve the PDEs for K = 8
[  t1_K8, e_L_K8  ] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 8, sigma, gamma, gamma_in,f,18,delta); % sigma≠0, with comm
[  t2_K8, e_C_K8 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, sigma, gamma,gamma_in, f, 18,delta); % sigma≠0, no comm
[ t11_K8, e_L1_K8] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 8, 0, gamma,gamma_in, f, 18,delta);     % sigma=0, with comm
[ t22_K8, e_C1_K8] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, 0, gamma,gamma_in, f, 18,delta);     % sigma=0, no

%%  For Different number of Leaders

[  t2_K14_15, e_C_K14_15 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 1.4, sigma, gamma,gamma_in, f, 15,delta);
[  t2_K14_17, e_C_K14_17 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 1.4, sigma, gamma,gamma_in, f, 17,delta);
[  t2_K14_19, e_C_K14_19] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 1.4, sigma, gamma,gamma_in, f, 19,delta);


[  t2_K3_14, e_C_K3_14 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 3, sigma, gamma,gamma_in, f, 14,delta);
[  t2_K3_19, e_C_K3_19 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 3, sigma, gamma,gamma_in, f, 19,delta);
[  t2_K3_25, e_C_K3_25 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 3, sigma, gamma,gamma_in, f, 25,delta);

[  t2_K5_16, e_C_K5_16 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, sigma, gamma,gamma_in, f, 16,delta);
[  t2_K5_19, e_C_K5_19 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, sigma, gamma,gamma_in, f, 19,delta);
[  t2_K5_25, e_C_K5_25 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, sigma, gamma,gamma_in, f, 25,delta);

[  t2_K8_18, e_C_K8_18 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, sigma, gamma,gamma_in, f, 18,delta);
[  t2_K8_20, e_C_K8_20 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, sigma, gamma,gamma_in, f, 20,delta);
[  t2_K8_25, e_C_K8_25 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, sigma, gamma,gamma_in, f, 25,delta);


%%
[ t1A, e_LA ] = OJAG_solveParabolicPDE_vec(N, M, T, a1, kappa1, K1, sigma1, gamma, gamma_in,f1, num_leaders,delta);
[  t2A,  e_CA] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a1, kappa1, K1, sigma1, gamma,gamma_in, f1, num_leaders,delta);
[ t11A, e_L1A] = OJAG_solveParabolicPDE_vec(N, M, T, a1, kappa1, K1, 0, gamma,gamma_in, f1, num_leaders,delta);
[t22A, e_C1A ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a1, kappa1, K1, 0, gamma, gamma_in,f1, num_leaders,delta);



% For Figure 2, we need to solve PDE with different no of leaders 14, 12,
% 10
[t22A_14, e_C1A_14 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a1, kappa1, K1, 0, gamma, gamma_in,f1, 15,delta);
[t22A_13, e_C1A_13 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a1, kappa1, K1, 0, gamma, gamma_in,f1,16,delta);
[t22A_12, e_C1A_12 ] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a1, kappa1, K1,0, gamma, gamma_in,f1, 19,delta);

% For Figure 5, we are analyzing the benefit of sigma
[ t_sigma_2, e_sigma_2 ] = OJAG_solveParabolicPDE_vec(N, M, T, a1, kappa1, K1, sigma1+1, gamma, gamma_in,f1, num_leaders,delta);
[ t_sigma_4, e_sigma_4 ] = OJAG_solveParabolicPDE_vec(N, M, T, a1, kappa1, K1, sigma1 +0.5, gamma, gamma_in,f1, num_leaders,delta);
[ t_sigma_6, e_sigma_6 ] = OJAG_solveParabolicPDE_vec(N, M, T, a1, kappa1, K1, sigma1-0.5, gamma, gamma_in,f1, num_leaders,delta);

%% Correct vector-valued L2 norm computation for agent-major format
%compute_vec_L2 = @(E) arrayfun(@(t) ...
    %sqrt(sum(arrayfun(@(i) norm(E((i-1)*3+1:i*3, t))^2, 1:N+1)) * dx), ...
    %1:size(E,2));
% K = 1.4
L2_norm_L   = compute_vec_L2(e_L,dx);  log_L2_norm_L   = log(L2_norm_L);
L2_norm_C   = compute_vec_L2(e_C,dx);  log_L2_norm_C   = log(L2_norm_C);
L2_norm_L1  = compute_vec_L2(e_L1,dx);  log_L2_norm_L1  = log(L2_norm_L1);
L2_norm_C1  = compute_vec_L2(e_C1,dx);  log_L2_norm_C1  = log(L2_norm_C1);

% K = 3
L2_norm_L_K3   = compute_vec_L2(e_L_K3,dx);  log_L2_norm_L_K3   = log(L2_norm_L_K3);
L2_norm_C_K3   = compute_vec_L2(e_C_K3,dx);  log_L2_norm_C_K3   = log(L2_norm_C_K3);
L2_norm_L1_K3  = compute_vec_L2(e_L1_K3,dx);  log_L2_norm_L1_K3  = log(L2_norm_L1_K3);
L2_norm_C1_K3  = compute_vec_L2(e_C1_K3,dx);  log_L2_norm_C1_K3  = log(L2_norm_C1_K3);

% K = 5
L2_norm_L_K5   = compute_vec_L2(e_L_K5,dx);  log_L2_norm_L_K5   = log(L2_norm_L_K5);
L2_norm_C_K5   = compute_vec_L2(e_C_K5,dx);  log_L2_norm_C_K5   = log(L2_norm_C_K5);
L2_norm_L1_K5  = compute_vec_L2(e_L1_K5,dx);  log_L2_norm_L1_K5  = log(L2_norm_L1_K5);
L2_norm_C1_K5  = compute_vec_L2(e_C1_K5,dx);  log_L2_norm_C1_K5  = log(L2_norm_C1_K5);


% K = 8
L2_norm_L_K8   = compute_vec_L2(e_L_K8,dx);  log_L2_norm_L_K8   = log(L2_norm_L_K8);
L2_norm_C_K8   = compute_vec_L2(e_C_K8,dx);  log_L2_norm_C_K8   = log(L2_norm_C_K8);
L2_norm_L1_K8  = compute_vec_L2(e_L1_K8,dx);  log_L2_norm_L1_K8  = log(L2_norm_L1_K8);
L2_norm_C1_K8  = compute_vec_L2(e_C1_K8,dx);  log_L2_norm_C1_K8  = log(L2_norm_C1_K8);

%% For different number of Leaders

L2_norm_C_K14_15   = compute_vec_L2(e_C_K14_15,dx);  log_L2_norm_C_K14_15   = log(L2_norm_C_K14_15);
L2_norm_C_K14_17   = compute_vec_L2(e_C_K14_17,dx);  log_L2_norm_C_K14_17   = log(L2_norm_C_K14_17);
L2_norm_C_K14_19   = compute_vec_L2(e_C_K14_19,dx);  log_L2_norm_C_K14_19   = log(L2_norm_C_K14_19);

L2_norm_C_K3_14   = compute_vec_L2(e_C_K3_14,dx);  log_L2_norm_C_K3_14   = log(L2_norm_C_K3_14);
L2_norm_C_K3_19   = compute_vec_L2(e_C_K3_19,dx);  log_L2_norm_C_K3_19   = log(L2_norm_C_K3_19);
L2_norm_C_K3_25   = compute_vec_L2(e_C_K3_25,dx);  log_L2_norm_C_K3_25  = log(L2_norm_C_K3_25);

L2_norm_C_K5_16   = compute_vec_L2(e_C_K5_16,dx);  log_L2_norm_C_K5_16   = log(L2_norm_C_K5_16);
L2_norm_C_K5_19   = compute_vec_L2(e_C_K5_19,dx);  log_L2_norm_C_K5_19   = log(L2_norm_C_K5_19);
L2_norm_C_K5_25   = compute_vec_L2(e_C_K5_25,dx);  log_L2_norm_C_K5_25   = log(L2_norm_C_K5_25);

L2_norm_C_K8_18   = compute_vec_L2(e_C_K8_18,dx);  log_L2_norm_C_K8_18   = log(L2_norm_C_K8_18);
L2_norm_C_K8_20   = compute_vec_L2(e_C_K8_20,dx);  log_L2_norm_C_K8_20   = log(L2_norm_C_K8_20);
L2_norm_C_K8_25   = compute_vec_L2(e_C_K8_25,dx);  log_L2_norm_C_K8_25   = log(L2_norm_C_K8_25);



%%
% a = 0.001
L2_norm_LA   = compute_vec_L2(e_LA,dx); log_L2_norm_LA   = log(L2_norm_LA);
L2_norm_CA   = compute_vec_L2(e_CA,dx); log_L2_norm_CA   = log(L2_norm_CA);
L2_norm_L1A  = compute_vec_L2(e_L1A,dx); log_L2_norm_L1A  = log(L2_norm_L1A);
L2_norm_C1A  = compute_vec_L2(e_C1A,dx); log_L2_norm_C1A  = log(L2_norm_C1A);

L2_norm_C1A_14  = compute_vec_L2(e_C1A_14,dx); log_L2_norm_C1A_14  = log(L2_norm_C1A_14); % for 14 Leaders
L2_norm_C1A_13  = compute_vec_L2(e_C1A_13,dx); log_L2_norm_C1A_13  = log(L2_norm_C1A_13); % for 13 Leaders
L2_norm_C1A_12  = compute_vec_L2(e_C1A_12,dx); log_L2_norm_C1A_12  = log(L2_norm_C1A_12); % for 12 Leaders

% For Sigma = 2,4,6,8
L2_norm_sigma_2  = compute_vec_L2(e_sigma_2,dx); log_L2_norm_sigma_2  = log(L2_norm_sigma_2);
L2_norm_sigma_4  = compute_vec_L2(e_sigma_4,dx); log_L2_norm_sigma_4  = log(L2_norm_sigma_4);
L2_norm_sigma_6  = compute_vec_L2(e_sigma_6,dx); log_L2_norm_sigma_6  = log(L2_norm_sigma_6);

%% Figure 1: Plotting log(||e||) over time
figure;
% First subplot: K = 1.4
subplot(2,2,1);
plot(t1, log_L2_norm_L,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t2, log_L2_norm_C,  'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
plot(t11, log_L2_norm_L1, 'Color', 'red', 'LineWidth', 2);
plot(t22, log_L2_norm_C1, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\log\left(\|e(t,\cdot)\|\right)$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs.\ $t$: $K = 1.4$, and $14$ leaders', 'Interpreter', 'latex');
legend({'$\sigma \neq 0$ with comm.', ...
        '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.', ...
        '$\sigma = 0$ no comm.'}, ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
grid on;

% Second subplot: K = 3
subplot(2,2,2);
plot(t1_K3, log_L2_norm_L_K3,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t2_K3, log_L2_norm_C_K3,  'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
plot(t11_K3, log_L2_norm_L1_K3, 'Color', 'red', 'LineWidth', 2);
plot(t22_K3, log_L2_norm_C1_K3, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\log\left(\|e(t,\cdot)\|\right)$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs.\ $t$: $K = 3$, and $14$ leaders', 'Interpreter', 'latex');
legend({'$\sigma \neq 0$ with comm.', ...
        '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.', ...
        '$\sigma = 0$ no comm.'}, ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
grid on;

subplot(2,2,3);
plot(t1_K5, log_L2_norm_L_K5,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t2_K5, log_L2_norm_C_K5,  'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
plot(t11_K5, log_L2_norm_L1_K5, 'Color', 'red', 'LineWidth', 2);
plot(t22_K5, log_L2_norm_C1_K5, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\log\left(\|e(t,\cdot)\|\right)$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs.\ $t$: $K = 5$, and $16$ leaders', 'Interpreter', 'latex');
legend({'$\sigma \neq 0$ with comm.', ...
        '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.', ...
        '$\sigma = 0$ no comm.'}, ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
grid on;

subplot(2,2,4);
plot(t1_K8, log_L2_norm_L_K8,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t2_K8, log_L2_norm_C_K8,  'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
plot(t11_K8, log_L2_norm_L1_K8, 'Color', 'red', 'LineWidth', 2);
plot(t22_K8, log_L2_norm_C1_K8, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\log\left(\|e(t,\cdot)\|\right)$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs.\ $t$: $K = 8$, and $18$ leaders', 'Interpreter', 'latex');
legend({'$\sigma \neq 0$ with comm.', ...
        '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.', ...
        '$\sigma = 0$ no comm.'}, ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
grid on;

%% Figure 1: Plotting log(||e||) over time
figure;

% First subplot: K = 1.4
subplot(1,2,1);
plot(t1, log_L2_norm_L,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t2, log_L2_norm_C,  'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
plot(t11, log_L2_norm_L1, 'Color', 'red', 'LineWidth', 2);
plot(t22, log_L2_norm_C1, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\log\left(\|e(t,\cdot)\|\right)$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs.\ $t$: $a = 0.004$', 'Interpreter', 'latex');
legend({'$\sigma \neq 0$ with comm.', ...
        '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.', ...
        '$\sigma = 0$ no comm.'}, ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
grid on;

% Second subplot: K = 3
subplot(1,2,2);
plot(t1A, log_L2_norm_LA,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t2A, log_L2_norm_CA,  'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
plot(t11A, log_L2_norm_L1A, 'Color', 'red', 'LineWidth', 2);
plot(t22A, log_L2_norm_C1A, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\log\left(\|e(t,\cdot)\|\right)$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs.\ $t$: $a = 0.01$', 'Interpreter', 'latex');
legend({'$\sigma \neq 0$ with comm.', ...
        '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.', ...
        '$\sigma = 0$ no comm.'}, ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
grid on;

%% Figure 2 Plotting ||e(t,.)|| over time (linear scale)
figure;

% First subplot: a = 0.004
subplot(2,1,1);
plot(t1, L2_norm_L,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t2, L2_norm_C,  'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
plot(t11, L2_norm_L1, 'Color', 'red', 'LineWidth', 2);
plot(t22, L2_norm_C1, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\|e(t,\cdot)\|$', 'Interpreter', 'latex');
title('$\|e(t,\cdot)\|$ vs.\ $t$: $a = 0.004$', 'Interpreter', 'latex');
legend({'$\sigma \neq 0$ with comm.', ...
        '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.', ...
        '$\sigma = 0$ no comm.'}, ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
grid on;

% Second subplot: a = 0.01
subplot(2,1,2);
plot(t1A, L2_norm_LA,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t2A, L2_norm_CA,  'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
plot(t11A, L2_norm_L1A, 'Color', 'red', 'LineWidth', 2);
plot(t22A, L2_norm_C1A, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);

xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\|e(t,\cdot)\|$', 'Interpreter', 'latex');
title('$\|e(t,\cdot)\|$ vs.\ $t$: $a = 0.01$', 'Interpreter', 'latex');
legend({'$\sigma \neq 0$ with comm.', ...
        '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.', ...
        '$\sigma = 0$ no comm.'}, ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
grid on;


%% Figure 3: Correct Phase Portrait: Plot actual agent states z = e + gamma

num_time_steps = size(e_L, 2);
Z = zeros(N+1, 3, num_time_steps); % Agent state z(t,x)

% Compute gamma(x) at spatial nodes
X = linspace(0, 1, N+1);
gamma_profile = gamma(X);  % 3 x (N+1)
gamma_in_profile = gamma_in(X);  % 3 x (N+1)

% Reconstruct agent states z(x,t) = e(x,t) + gamma(x)
for i = 1:num_time_steps
    e_i = e_L(:, i);  % Full error vector at time i
    e_vec = [e_i(1:3:end)'; e_i(2:3:end)'; e_i(3:3:end)'];  % 3 x (N+1)
    Z(:,:,i) = e_vec' + gamma_profile';  % (N+1) x 3
end

%% Plotting Phase Portrait
figure; 

% Plot target curve gamma(x) (static)
plot3(gamma_profile(1,:), gamma_profile(2,:), gamma_profile(3,:), ...
      'Color', [0.5, 0, 0], 'LineWidth', 2.5);  % Maroon
hold on
% Plot initial curve gamma_in(x) - dark yellow
%plot3(gamma_in_profile(1,:), gamma_in_profile(2,:), gamma_in_profile(3,:), ...
   %   'Color', [0.8, 0.6, 0], 'LineWidth', 2.5);
hold on;
% Plot interior agents' trajectories in black
for i = 2:N
    plot3(squeeze(Z(i,1,:)), squeeze(Z(i,2,:)), squeeze(Z(i,3,:)), ...
          'k', 'LineWidth', 1);
end
hold on
% Plot boundary agents' trajectories in orange
plot3(squeeze(Z(1,1,:)), squeeze(Z(1,2,:)), squeeze(Z(1,3,:)), ...
      'Color', [1, 0.5, 0], 'LineWidth', 2);  % Left boundary
plot3(squeeze(Z(end,1,:)), squeeze(Z(end,2,:)), squeeze(Z(end,3,:)), ...
      'Color', [1, 0.5, 0], 'LineWidth', 2);  % Right boundary
hold off;
% Labels and aesthetics
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
zlabel('$z$', 'Interpreter', 'latex');
title('Phase Portrait: Agent Trajectories $z(t,x)$', 'Interpreter', 'latex');
grid on;



%% Figure 4
%    \log(\|e(t, \cdot)\|^2)\) vs \(t\) Wentzell with (7 leaders)...
% or without leader (7,12 and 13 leaders) communication.
% Compartive analysis of communication strategy for \sigma \neq 0
  figure;
    plot(t1A, log_L2_norm_LA , 'Color', [1, 0.5, 0],'LineWidth', 3, 'DisplayName', 'Swarm Coupled Boundary + Comm. 7 leaders'); hold on;
    plot(t22A_14, log_L2_norm_C1A_14 , 'Color', 'blue', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 14 leaders');
    hold on;
    plot(t22A_13, log_L2_norm_C1A_13 , 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 16 leaders');
    hold on;
    plot(t22A_12, log_L2_norm_C1A_12 , 'Color', 'red','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 17 leaders');
   xlabel('$$ t $$', 'Interpreter', 'latex');
    ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
    title('Comparison of Square of L^2 norm: communicating and non communicating with different leaders');
    legend('show');
    grid on;

%% Figure 4. A
%    \log(\|e(t, \cdot)\|^2)\) vs \(t\) Wentzell with (7 leaders)...
% or without leader (7,12 and 13 leaders) communication.
% Compartive analysis of communication strategy for \sigma \neq 0
  figure;
  subplot(2,2,1);
    plot(t1,  log_L2_norm_L , 'Color', [1, 0.5, 0],'LineWidth', 3, 'DisplayName', 'Swarm Coupled Boundary + Comm. 7 leaders'); hold on;
    plot( t2_K14_15, log_L2_norm_C_K14_15 , 'Color', 'blue', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 14 leaders');
    hold on;
    plot( t2_K14_17, log_L2_norm_C_K14_17 , 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 16 leaders');
    hold on;
    plot( t2_K14_19, log_L2_norm_C_K14_19 , 'Color', 'red','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 17 leaders');
   xlabel('$$ t $$', 'Interpreter', 'latex');
    ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
    title('Comparison of Square of L^2 norm: communicating and non communicating with different leaders');
    legend('show');
    grid on;
subplot(2,2,2);
    plot( t1_K3,   log_L2_norm_L_K3  , 'Color', [1, 0.5, 0],'LineWidth', 3, 'DisplayName', 'Swarm Coupled Boundary + Comm. 7 leaders'); hold on;
    plot( t2_K3_14, log_L2_norm_C_K3_14 , 'Color', 'blue', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 14 leaders');
    hold on;
    plot( t2_K3_19, log_L2_norm_C_K3_19 , 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 16 leaders');
    hold on;
    plot( t2_K3_25, log_L2_norm_C_K3_25 , 'Color', 'red','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 17 leaders');
   xlabel('$$ t $$', 'Interpreter', 'latex');
    ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
    title('Comparison of Square of L^2 norm: communicating and non communicating with different leaders');
    legend('show');
    grid on;
subplot(2,2,3);
    plot( t1_K5,   log_L2_norm_L_K5  , 'Color', [1, 0.5, 0],'LineWidth', 3, 'DisplayName', 'Swarm Coupled Boundary + Comm. 7 leaders'); hold on;
    plot( t2_K5_16, log_L2_norm_C_K5_16 , 'Color', 'blue', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 14 leaders');
    hold on;
    plot( t2_K5_19, log_L2_norm_C_K5_19 , 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 16 leaders');
    hold on;
    plot( t2_K5_25, log_L2_norm_C_K5_25 , 'Color', 'red','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 17 leaders');
   xlabel('$$ t $$', 'Interpreter', 'latex');
    ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
    title('Comparison of Square of L^2 norm: communicating and non communicating with different leaders');
    legend('show');
    grid on;
subplot(2,2,4);
    plot( t1_K8,   log_L2_norm_L_K8 , 'Color', [1, 0.5, 0],'LineWidth', 3, 'DisplayName', 'Swarm Coupled Boundary + Comm. 7 leaders'); hold on;
    plot( t2_K8_18, log_L2_norm_C_K8_18 , 'Color', 'blue', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 14 leaders');
    hold on;
    plot( t2_K8_20, log_L2_norm_C_K8_20 , 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 16 leaders');
    hold on;
    plot( t2_K8_25, log_L2_norm_C_K8_25 , 'Color', 'red','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 17 leaders');
    xlabel('$$ t $$', 'Interpreter', 'latex');
    ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
    title('Comparison of Square of L^2 norm: communicating and non communicating with different leaders');
    legend('show');
    grid on;




%% Figure 5
 figure;
    plot( t1A, log_L2_norm_LA , 'Color', [1, 0.5, 0],'LineWidth', 3, 'DisplayName', 'sigma original'); hold on;
    plot( t_sigma_2, log_L2_norm_sigma_2 , 'Color', 'blue', 'LineWidth', 3, 'DisplayName', 'sigma +1');
    hold on;
    plot( t_sigma_4, log_L2_norm_sigma_4 , 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'sigma +0.5');
    hold on;
    plot( t_sigma_6, log_L2_norm_sigma_6 , 'Color', 'red','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'sigma -0.5');
   xlabel('$$ t $$', 'Interpreter', 'latex');
    ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
    title('Comparison of Square of L^2 norm: communicating and non communicating with different leaders');
    legend('show');
    grid on;

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







