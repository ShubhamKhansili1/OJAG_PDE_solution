%% System Parameter's
N = 70; num_leaders = 18; M = 10000; T = 4; kappa = 1; K = 8; 
sigma= 6; dx = 1 / N; h = 1/N; a = 0.001; delta =1 ; f = @(t, z) 0.3* z; 
% --- Target curve -------------------------------------------
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];   
%gamma = @(x) [1.5*sin(x); 1.5*sin(x);1.5*sin(x)];  

%gamma_in = @(x) gamma(x) .* (x == 0); 
gamma_in = @(x) gamma(x) .* ((x == 0) | (x == 1));
%gamma_in = @(x) repmat(gamma(0), 1, numel(x));  % All agents at gamma(0)
%gamma_in = @(x)[ 0*ones(size(x));  0*ones(size(x));  0*ones(size(x))];
%gamma = @(x) [h*cos(pi*N*x/3); h*sin(pi*N*x/3); 0*ones(size(x))]; 

%Solving the error PDE-system
% Case 1: When interaction strength is 0.0025  
[t1_1, e1_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
[t1_2, e1_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, sigma, gamma,gamma_in, f, num_leaders,delta);% Boundary intact + Nocom. Leader
[t1_3, e1_3] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + Comm. leader
[t1_4, e1_4] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + NoComm. Leader

% Case 2: When interaction strength is 0.001
[t2_1, e2_1] = OJAG_solveParabolicPDE_vec(N, M, T, 0.001, kappa, K, sigma, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
[t2_2, e2_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, 0.001, kappa, K, sigma, gamma,gamma_in, f, num_leaders,delta);% Boundary intact + Nocom. Leader
[t2_3, e2_3] = OJAG_solveParabolicPDE_vec(N, M, T, 0.001, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + Comm. leader
[t2_4, e2_4] = OJAG_solveParabolicPDE_constant_vec(N, M, T, 0.001, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + NoComm. Leader

%Computing the Log of L^2-norm for case 1
L2_norm1_1 = compute_vec_L2(e1_1,dx);  log_L2_norm1_1  = log(L2_norm1_1);
L2_norm1_2 = compute_vec_L2(e1_2,dx);  log_L2_norm1_2  = log(L2_norm1_2);
L2_norm1_3 = compute_vec_L2(e1_3,dx);  log_L2_norm1_3  = log(L2_norm1_3);
L2_norm1_4 = compute_vec_L2(e1_4,dx);  log_L2_norm1_4  = log(L2_norm1_4);

%Computing the Log of L^2-norm for case 2
L2_norm2_1 = compute_vec_L2(e2_1,dx);  log_L2_norm2_1  = log(L2_norm2_1);
L2_norm2_2 = compute_vec_L2(e2_2,dx);  log_L2_norm2_2  = log(L2_norm2_2);
L2_norm2_3 = compute_vec_L2(e2_3,dx);  log_L2_norm2_3  = log(L2_norm2_3);
L2_norm2_4 = compute_vec_L2(e2_4,dx);  log_L2_norm2_4  = log(L2_norm2_4);



% Figure 1: Plotting log(||e||) over time
figure;
plot(t1_1, log_L2_norm1_1,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t1_2, log_L2_norm1_2,  'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
plot(t1_3, log_L2_norm1_3, 'Color', 'red', 'LineWidth', 2);
plot(t1_4, log_L2_norm1_4, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs. $t$: $a = 0.0025$','Interpreter','latex');
legend({ ...
    'Strategy~3: $\sigma = 6$, with comm.', ...
    'Strategy~2: $\sigma = 6$, no comm.', ...
    'Strategy~1: $\sigma = 0$, with comm.', ...
    'Strategy~0: $\sigma = 0$, no comm.'}, ...
    'Location', 'best', 'FontSize', 25, 'NumColumns', 1, 'Interpreter', 'latex');


grid on;

% Second subplot (bottom or right plot)
figure;  % 2 rows, 1 column, second plot

plot(t2_1, log_L2_norm2_1, 'Color', [1, 0.5, 0], 'LineWidth', 2); 
hold on;
plot(t2_2, log_L2_norm2_2, 'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2); 
hold on;
plot(t2_3, log_L2_norm2_3, 'Color', 'red', 'LineWidth', 2); 
hold on;
plot(t2_4, log_L2_norm2_4, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs. $t$: $a = 0.001$','Interpreter','latex');
legend({ ...
    'Strategy~3: $\sigma = 6$, with comm.', ...
    'Strategy~2: $\sigma = 6$, no comm.', ...
    'Strategy~1: $\sigma = 0$, with comm.', ...
    'Strategy~0: $\sigma = 0$, no comm.'}, ...
    'Location', 'best', 'FontSize', 25, 'NumColumns', 2, 'Interpreter', 'latex');
grid on;

% Figure 1: Plotting log(||e||) over time
% Figure 1: Plotting log(||e||) over time (one shared legend)
figure;
tiledlayout(2,1);

% First subplot
nexttile;
p1 = plot(t1_1, log_L2_norm1_1,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
p2 = plot(t1_2, log_L2_norm1_2,  'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
p3 = plot(t1_3, log_L2_norm1_3, 'Color', 'red', 'LineWidth', 2);
p4 = plot(t1_4, log_L2_norm1_4, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\|e(t, \cdot)\|$$', 'Interpreter', 'latex');
title('$\log(\|e(t,\cdot)\|)$ vs. $t$: $a = 0.0025$', 'Interpreter', 'latex');
grid on;

% Second subplot
nexttile;
plot(t2_1, log_L2_norm2_1, 'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t2_2, log_L2_norm2_2, 'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 2);
plot(t2_3, log_L2_norm2_3, 'Color', 'red', 'LineWidth', 2);
plot(t2_4, log_L2_norm2_4, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\|e(t, \cdot)\|$$', 'Interpreter', 'latex');
title('$\log(\|e(t,\cdot)\|)$ vs. $t$: $a = 0.001$', 'Interpreter', 'latex');
grid on;

% Shared legend
lgd = legend([p1 p2 p3 p4], ...
    {'$\sigma = 6$ with communication', ...
     '$\sigma = 6$ with no communication', ...
     '$\sigma = 0$ with communication', ...
     '$\sigma = 0$ with no communication'}, ...
    'Interpreter', 'latex', 'NumColumns', 2, 'FontSize', 25);
lgd.Layout.Tile = 'south';



 