% This code shows the effect of the global strength K for K = 1.4, 3, 5 and 8 

% System Parameter's
N = 70; num_leaders = 14; M = 10000; T = 4; kappa = 1; K = 1.4; 
sigma= 6; dx = 1 / N; h = 1/N; a = 0.001; delta =1 ; f = @(t, z) 0.3* z; 
% --- Target curve -------------------------------------------
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];   
%gamma = @(x) [1.5*sin(x); 1.5*sin(x);1.5*sin(x)];  

%gamma_in = @(x) gamma(x) .* (x == 0); 
%gamma_in = @(x) gamma(x) .* ((x == 0) | (x == 1));
%gamma_in = @(x) repmat(gamma(0), 1, numel(x));  % All agents at gamma(0)
gamma_in = @(x)[ 0*ones(size(x));  0*ones(size(x));  0*ones(size(x))];
%gamma = @(x) [h*cos(pi*N*x/3); h*sin(pi*N*x/3); 0*ones(size(x))]; 

%Solving the error PDE-system
% Case 1: When Global strength is K = 1.4 (14 leaders)
[t3_1, e3_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
[t3_2, e3_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, sigma, gamma,gamma_in, f, num_leaders,delta);% Boundary intact + Nocom. Leader
[t3_3, e3_3] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + Comm. leader
[t3_4, e3_4] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + NoComm. Leader

% Case 2: When interaction strength is K = 3 (14 leaders)
[t4_1, e4_1] = OJAG_solveParabolicPDE_vec(N, M, T, 0.001, kappa, 3, sigma, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
[t4_2, e4_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, 0.001, kappa, 3, sigma, gamma,gamma_in, f, num_leaders,delta);% Boundary intact + Nocom. Leader
[t4_3, e4_3] = OJAG_solveParabolicPDE_vec(N, M, T, 0.001, kappa, 3, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + Comm. leader
[t4_4, e4_4] = OJAG_solveParabolicPDE_constant_vec(N, M, T, 0.001, kappa, 3, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + NoComm. Leader

% Case 3: When Global strength is K = 5 (16 leaders)
[t5_1, e5_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 5, sigma, gamma, gamma_in,f,16,delta);% Boundary intact + Comm. Leader
[t5_2, e5_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, sigma, gamma,gamma_in, f, 16,delta);% Boundary intact + Nocom. Leader
[t5_3, e5_3] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 5, 0, gamma,gamma_in, f, 16,delta);% Independent Bo. + Comm. leader
[t5_4, e5_4] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 5, 0, gamma,gamma_in, f, 16,delta);% Independent Bo. + NoComm. Leader

% Case 4: When Global strength is K = 8 (18 leaders)
[t6_1, e6_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 8, sigma, gamma, gamma_in,f,18,delta);% Boundary intact + Comm. Leader
[t6_2, e6_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, sigma, gamma,gamma_in, f, 18,delta);% Boundary intact + Nocom. Leader
[t6_3, e6_3] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, 8, 0, gamma,gamma_in, f, 18,delta);% Independent Bo. + Comm. leader
[t6_4, e6_4] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, 8, 0, gamma,gamma_in, f, 18,delta);% Independent Bo. + NoComm. Leader

%Computing the Log of L^2-norm for case 1
L2_norm3_1 = compute_vec_L2(e3_1,dx);  log_L2_norm3_1  = log(L2_norm3_1);
L2_norm3_2 = compute_vec_L2(e3_2,dx);  log_L2_norm3_2  = log(L2_norm3_2);
L2_norm3_3 = compute_vec_L2(e3_3,dx);  log_L2_norm3_3  = log(L2_norm3_3);
L2_norm3_4 = compute_vec_L2(e3_4,dx);  log_L2_norm3_4  = log(L2_norm3_4);

%Computing the Log of L^2-norm for case 2
L2_norm4_1 = compute_vec_L2(e4_1,dx);  log_L2_norm4_1  = log(L2_norm4_1);
L2_norm4_2 = compute_vec_L2(e4_2,dx);  log_L2_norm4_2  = log(L2_norm4_2);
L2_norm4_3 = compute_vec_L2(e4_3,dx);  log_L2_norm4_3  = log(L2_norm4_3);
L2_norm4_4 = compute_vec_L2(e4_4,dx);  log_L2_norm4_4  = log(L2_norm4_4);

%Computing the Log of L^2-norm for case 3
L2_norm5_1 = compute_vec_L2(e5_1,dx);  log_L2_norm5_1  = log(L2_norm5_1);
L2_norm5_2 = compute_vec_L2(e5_2,dx);  log_L2_norm5_2  = log(L2_norm5_2);
L2_norm5_3 = compute_vec_L2(e5_3,dx);  log_L2_norm5_3  = log(L2_norm5_3);
L2_norm5_4 = compute_vec_L2(e5_4,dx);  log_L2_norm5_4  = log(L2_norm5_4);

%Computing the Log of L^2-norm for case 4
L2_norm6_1 = compute_vec_L2(e6_1,dx);  log_L2_norm6_1  = log(L2_norm6_1);
L2_norm6_2 = compute_vec_L2(e6_2,dx);  log_L2_norm6_2  = log(L2_norm6_2);
L2_norm6_3 = compute_vec_L2(e6_3,dx);  log_L2_norm6_3  = log(L2_norm6_3);
L2_norm6_4 = compute_vec_L2(e6_4,dx);  log_L2_norm6_4  = log(L2_norm6_4);

figure;
tiledlayout(2,2);

% Case 1
nexttile;
p1 = plot(t3_1,log_L2_norm3_1,'Color',[1 0.5 0],'LineWidth',2); hold on;
p2 = plot(t3_2,log_L2_norm3_2,'Color',[1 0.5 0],'LineStyle','--','LineWidth',2);
p3 = plot(t3_3,log_L2_norm3_3,'r','LineWidth',2);
p4 = plot(t3_4,log_L2_norm3_4,'r--','LineWidth',2);
xlabel('$$ t $$','Interpreter','latex');
ylabel('$$\log\|e(t,\cdot)\|$$','Interpreter','latex');
title('$K=1.4$, 14 leaders','Interpreter','latex'); grid on;

% Case 2
nexttile;
plot(t4_1,log_L2_norm4_1,'Color',[1 0.5 0],'LineWidth',2); hold on;
plot(t4_2,log_L2_norm4_2,'Color',[1 0.5 0],'LineStyle','--','LineWidth',2);
plot(t4_3,log_L2_norm4_3,'r','LineWidth',2);
plot(t4_4,log_L2_norm4_4,'r--','LineWidth',2);
xlabel('$$ t $$','Interpreter','latex');
ylabel('$$\log\|e(t,\cdot)\|$$','Interpreter','latex');
title('$K=3$, 14 leaders','Interpreter','latex'); grid on;

% Case 3
nexttile;
plot(t5_1,log_L2_norm5_1,'Color',[1 0.5 0],'LineWidth',2); hold on;
plot(t5_2,log_L2_norm5_2,'Color',[1 0.5 0],'LineStyle','--','LineWidth',2);
plot(t5_3,log_L2_norm5_3,'r','LineWidth',2);
plot(t5_4,log_L2_norm5_4,'r--','LineWidth',2);
xlabel('$$ t $$','Interpreter','latex');
ylabel('$$\log\|e(t,\cdot)\|$$','Interpreter','latex');
title('$K=5$, 16 leaders','Interpreter','latex'); grid on;

% Case 4
nexttile;
plot(t6_1,log_L2_norm6_1,'Color',[1 0.5 0],'LineWidth',2); hold on;
plot(t6_2,log_L2_norm6_2,'Color',[1 0.5 0],'LineStyle','--','LineWidth',2);
plot(t6_3,log_L2_norm6_3,'r','LineWidth',2);
plot(t6_4,log_L2_norm6_4,'r--','LineWidth',2);
xlabel('$$ t $$','Interpreter','latex');
ylabel('$$\log\|e(t,\cdot)\|$$','Interpreter','latex');
title('$K=8$, 18 leaders','Interpreter','latex'); grid on;

% One common legend
lgd = legend([p1 p2 p3 p4], ...
    {'Strategy~3: $\sigma = 6$, with comm.', ...
     'Strategy~2: $\sigma = 6$, no comm.', ...
     'Strategy~1: $\sigma = 0$, with comm.', ...
     'Strategy~0: $\sigma = 0$, no comm.'}, ...
    'Interpreter','latex','NumColumns',2,'FontSize',14);
lgd.Layout.Tile = 'south';
