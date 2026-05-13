% This code shows that more leaders are needed to match the performance of the system when leaders are communicating and with wentzell boundary control.  

% System Parameter's
N = 70; num_leaders = 14; M = 10000; T = 4; kappa = 1; K = 1.4; 
sigma= 6; dx = 1 / N; h = 1/N; a = 0.001; delta =1 ; f = @(t, z) 0.3* z; 
% Do not change the value of a in this code. 
% --- Target curve -------------------------------------------
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];   
%gamma = @(x) [1.5*sin(x); 1.5*sin(x);1.5*sin(x)];  

%gamma_in = @(x) gamma(x) .* (x == 0); 
%gamma_in = @(x) gamma(x) .* ((x == 0) | (x == 1));
%gamma_in = @(x) repmat(gamma(0), 1, numel(x));  % All agents at gamma(0)
gamma_in = @(x)[ 0*ones(size(x));  0*ones(size(x));  0*ones(size(x))];
%gamma = @(x) [h*cos(pi*N*x/3); h*sin(pi*N*x/3); 0*ones(size(x))]; 

%Solving the error PDE-system for K = 1.4
% Wentzell + 14 communicating leaders
[t7_1, e7_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
% No Wentzell + 14 Non comm leaders
[t7_2, e7_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + NoComm. Leader
% No Wentzell +15 Non Comm Leaders
[t7_3, e7_3] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, 15,delta);% Independent Bo. + NoComm. Leader
% No Wentzell +17 Non Comm Leaders
[t7_4, e7_4] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, 17,delta);% Independent Bo. + NoComm. Leader
% No Wentzell +19 Non Comm Leaders
[t7_5, e7_5] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, 19,delta);% Independent Bo. + NoComm. Leader
% No Wentzell + 22 Non-comm Leaders
[t7_6, e7_6] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, 22,delta);% Independent Bo. + NoComm. Leader
% Computing the Log of L^2-norm for case 1
L2_norm7_1 = compute_vec_L2(e7_1,dx);  log_L2_norm7_1  = log(L2_norm7_1);
L2_norm7_2 = compute_vec_L2(e7_2,dx);  log_L2_norm7_2  = log(L2_norm7_2);
L2_norm7_3 = compute_vec_L2(e7_3,dx);  log_L2_norm7_3  = log(L2_norm7_3);
L2_norm7_4 = compute_vec_L2(e7_4,dx);  log_L2_norm7_4  = log(L2_norm7_4);
L2_norm7_5 = compute_vec_L2(e7_5,dx);  log_L2_norm7_5  = log(L2_norm7_5);
L2_norm7_6 = compute_vec_L2(e7_6,dx);  log_L2_norm7_6  = log(L2_norm7_6);
figure;
plot(t7_1, log_L2_norm7_1 , 'Color', [1, 0.5, 0],'LineWidth', 3, 'DisplayName', 'Swarm Coupled Boundary + Comm. 7 leaders'); hold on;
plot(t7_2, log_L2_norm7_2 , 'Color', 'blue', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 7 leaders');
hold on;
plot(t7_3, log_L2_norm7_3 , 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 15 leaders');
hold on;
plot(t7_4, log_L2_norm7_4 , 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 17 leaders');
hold on;
plot(t7_5, log_L2_norm7_5 , 'Color', 'green','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 19 leaders');
hold on;
plot(t7_6, log_L2_norm7_6, 'Color', 'red','LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Swarm Uncoupled Boundary + Non-comm. 22 leaders');
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
title('Comparison of Square of L^2 norm: communicating and non communicating with different leaders');
legend('show');
grid on;

% Now we do same for different values of K = 3, 5, 8
%Solving the error PDE-system for K = 3
% Wentzell + 14 communicating leaders
[t8_1, e8_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 3, sigma, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
% No Wentzell + 14 Non comm leaders
[t8_2, e8_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 3, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + NoComm. Leader
% No Wentzell +15 Non Comm Leaders
[t8_3, e8_3] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 3, 0, gamma,gamma_in, f, 15,delta);% Independent Bo. + NoComm. Leader
% No Wentzell +17 Non Comm Leaders
[t8_4, e8_4] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 3, 0, gamma,gamma_in, f, 17,delta);% Independent Bo. + NoComm. Leader
% No Wentzell +19 Non Comm Leaders
[t8_5, e8_5] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 3, 0, gamma,gamma_in, f, 19,delta);% Independent Bo. + NoComm. Leader
% No Wentzell + 22 Non-comm Leaders
[t8_6, e8_6] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 3, 0, gamma,gamma_in, f, 22,delta);% Independent Bo. + NoComm. Leader
% Computing the Log of L^2-norm for case 1
L2_norm8_1 = compute_vec_L2(e8_1,dx);  log_L2_norm8_1  = log(L2_norm8_1);
L2_norm8_2 = compute_vec_L2(e8_2,dx);  log_L2_norm8_2  = log(L2_norm8_2);
L2_norm8_3 = compute_vec_L2(e8_3,dx);  log_L2_norm8_3  = log(L2_norm8_3);
L2_norm8_4 = compute_vec_L2(e8_4,dx);  log_L2_norm8_4  = log(L2_norm8_4);
L2_norm8_5 = compute_vec_L2(e8_5,dx);  log_L2_norm8_5  = log(L2_norm8_5);
L2_norm8_6 = compute_vec_L2(e8_6,dx);  log_L2_norm8_6  = log(L2_norm8_6);

%Solving the error PDE-system for K = 5
% Wentzell + 16 communicating leaders
[t9_1, e9_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 5, sigma, gamma, gamma_in,f,16,delta);% Boundary intact + Comm. Leader
% No Wentzell + 16 Non comm leaders
[t9_2, e9_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, 0, gamma,gamma_in, f, 16,delta);% Independent Bo. + NoComm. Leader
% No Wentzell +17 Non Comm Leaders
[t9_3, e9_3] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, 0, gamma,gamma_in, f, 17,delta);% Independent Bo. + NoComm. Leader
% No Wentzell +19 Non Comm Leaders
[t9_4, e9_4] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, 0, gamma,gamma_in, f, 19,delta);% Independent Bo. + NoComm. Leader
% No Wentzell + 21 Non Comm Leaders
[t9_5, e9_5] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, 0, gamma,gamma_in, f, 21,delta);% Independent Bo. + NoComm. Leader
% No Wentzell + 23 Non-comm Leaders
[t9_6, e9_6] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, 0, gamma,gamma_in, f, 23,delta);% Independent Bo. + NoComm. Leader
% Computing the Log of L^2-norm for case 1
L2_norm9_1 = compute_vec_L2(e9_1,dx);  log_L2_norm9_1  = log(L2_norm9_1);
L2_norm9_2 = compute_vec_L2(e9_2,dx);  log_L2_norm9_2  = log(L2_norm9_2);
L2_norm9_3 = compute_vec_L2(e9_3,dx);  log_L2_norm9_3  = log(L2_norm9_3);
L2_norm9_4 = compute_vec_L2(e9_4,dx);  log_L2_norm9_4  = log(L2_norm9_4);
L2_norm9_5 = compute_vec_L2(e9_5,dx);  log_L2_norm9_5  = log(L2_norm9_5);
L2_norm9_6 = compute_vec_L2(e9_6,dx);  log_L2_norm9_6  = log(L2_norm9_6);

%Solving the error PDE-system for K = 8
% Wentzell + 18 communicating leaders
[t10_1, e10_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 8, sigma, gamma, gamma_in,f,18,delta);% Boundary intact + Comm. Leader
% No Wentzell + 18 Non comm leaders
[t10_2, e10_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, 0, gamma,gamma_in, f, 18,delta);% Independent Bo. + NoComm. Leader
% No Wentzell +19 Non Comm Leaders
[t10_3, e10_3] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, 0, gamma,gamma_in, f, 19,delta);% Independent Bo. + NoComm. Leader
% No Wentzell + 21 Non Comm Leaders
[t10_4, e10_4] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, 0, gamma,gamma_in, f, 21,delta);% Independent Bo. + NoComm. Leader
% No Wentzell + 23 Non Comm Leaders
[t10_5, e10_5] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, 0, gamma,gamma_in, f, 23,delta);% Independent Bo. + NoComm. Leader
% No Wentzell + 25 Non-comm Leaders
[t10_6, e10_6] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, 0, gamma,gamma_in, f, 25,delta);% Independent Bo. + NoComm. Leader
% Computing the Log of L^2-norm for case 1
L2_norm10_1 = compute_vec_L2(e10_1,dx);  log_L2_norm10_1  = log(L2_norm10_1);
L2_norm10_2 = compute_vec_L2(e10_2,dx);  log_L2_norm10_2  = log(L2_norm10_2);
L2_norm10_3 = compute_vec_L2(e10_3,dx);  log_L2_norm10_3  = log(L2_norm10_3);
L2_norm10_4 = compute_vec_L2(e10_4,dx);  log_L2_norm10_4  = log(L2_norm10_4);
L2_norm10_5 = compute_vec_L2(e10_5,dx);  log_L2_norm10_5  = log(L2_norm10_5);
L2_norm10_6 = compute_vec_L2(e10_6,dx);  log_L2_norm10_6  = log(L2_norm10_6);

% Figure 
figure;
colors = [
    0, 0.4470, 0.7410;    % Blue
    0.4660, 0.6740, 0.1880;  % Green
    0.4940, 0.1840, 0.5560;  % Purple
    0.9290, 0.6940, 0.1250;  % Yellow-ish
    0.3010, 0.7450, 0.9330   % Cyan-ish
];
% Case 1 K = 1.4 for a = 0.001 and 14 leaders
subplot(2,2,1);     % 2 rows, 1 column, second plot
plot(t7_1, log_L2_norm7_1,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t7_2, log_L2_norm7_2, 'Color', colors(1,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t7_3, log_L2_norm7_3, 'Color', colors(2,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t7_4, log_L2_norm7_4, 'Color', colors(3,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t7_5, log_L2_norm7_5, 'Color', colors(4,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t7_6, log_L2_norm7_6, 'Color', colors(5,:), 'LineStyle', '--', 'LineWidth', 2);
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs. $t$: $K = 1.4$','Interpreter','latex');
legend({'Strategy 3 with 14 leaders', ...
         'Strategy 0 with 14 leaders', ...
         'Strategy 0 with 15 leaders', ...
         'Strategy 0 with 17 leaders', ...
        'Strategy 0 with 19 leaders', ...
        'Strategy 0 with 22 leaders'}, ...
        'Location', 'best', 'FontSize', 12, 'NumColumns', 1, 'Interpreter', 'latex');
grid on;

%case 2 K = 3 , a = 0.001 and 14 leaders
subplot(2,2,2);    

plot(t8_1, log_L2_norm8_1, 'Color', [1, 0.5, 0], 'LineWidth', 2); 
hold on;
plot(t8_2, log_L2_norm8_2, 'Color', colors(1,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t8_3, log_L2_norm8_3, 'Color', colors(2,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t8_4, log_L2_norm8_4, 'Color', colors(3,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t8_5, log_L2_norm8_5, 'Color', colors(4,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t8_6, log_L2_norm8_6, 'Color', colors(5,:), 'LineStyle', '--', 'LineWidth', 2);
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs. $t$: $K = 3$ ','Interpreter','latex');
legend({'Strategy 3 with 14 leaders', ...
         'Strategy 0 with 14 leaders', ...
         'Strategy 0 with 15 leaders', ...
         'Strategy 0 with 17 leaders', ...
        'Strategy 0 with 19 leaders', ...
        'Strategy 0 with 22 leaders'}, ...
        'Location', 'best', 'FontSize', 12, 'NumColumns', 1, 'Interpreter', 'latex');
grid on;

%Case 3; K = 5 for a = 0.001, with 16 leaders
subplot(2,2,3);     % 2 rows, 1 column, second plot
plot(t9_1, log_L2_norm9_1,  'Color', [1, 0.5, 0], 'LineWidth', 2); hold on;
plot(t9_2, log_L2_norm9_2, 'Color', colors(1,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t9_3, log_L2_norm9_3, 'Color', colors(2,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t9_4, log_L2_norm9_4, 'Color', colors(3,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t9_5, log_L2_norm9_5, 'Color', colors(4,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t9_6, log_L2_norm9_6, 'Color', colors(5,:), 'LineStyle', '--', 'LineWidth', 2);
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs. $t$: $K = 5$ ','Interpreter','latex');
legend({'Strategy 3 with 16 leaders', ...
         'Strategy 0 with 16 leaders', ...
         'Strategy 0 with 17 leaders', ...
         'Strategy 0 with 19 leaders', ...
        'Strategy 0 with 21 leaders', ...
        'Strategy 0 with 23 leaders'}, ...
        'Location', 'best', 'FontSize', 12, 'NumColumns', 1, 'Interpreter', 'latex');
grid on;

% Case 4: K =8 with a = 0.001 with 18 leaders
subplot(2,2,4);     % 2 rows, 1 column, second plot

plot(t10_1, log_L2_norm10_1, 'Color', [1, 0.5, 0], 'LineWidth', 2); 
hold on;
plot(t10_2, log_L2_norm10_2, 'Color', colors(1,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t10_3, log_L2_norm10_3, 'Color', colors(2,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t10_4, log_L2_norm10_4, 'Color', colors(3,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t10_5, log_L2_norm10_5, 'Color', colors(4,:), 'LineStyle', '--', 'LineWidth', 2);
plot(t10_6, log_L2_norm10_6, 'Color', colors(5,:), 'LineStyle', '--', 'LineWidth', 2);
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\log\left(\|e(t, \cdot)\|\right)$$', 'Interpreter', 'latex');
title('$\log\bigl(\|e(t,\cdot)\|\bigr)$ vs. $t$: $K = 8$ ','Interpreter','latex');
legend({'Strategy 3 with 18 leaders', ...
         'Strategy 0 with 18 leaders', ...
         'Strategy 0 with 19 leaders', ...
         'Strategy 0 with 21 leaders', ...
        'Strategy 0 with 23 leaders', ...
        'Strategy 0 with 25 leaders'}, ...
        'Location', 'best', 'FontSize', 12, 'NumColumns', 1, 'Interpreter', 'latex');
grid on;




