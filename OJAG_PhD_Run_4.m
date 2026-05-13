% This code checks the performance of the boundary agent
% when they are coupled and not with the swarm 

% System Parameter's
N = 70; num_leaders = 7; M = 10000; T = 4; kappa = 1; K = 1.4; 
sigma= 6; dx = 1 / N; h = 1/N; a = 0.004; delta =1 ; f = @(t, z) 0.3* z; 
% --- Target curve -------------------------------------------
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];   
gamma_in = @(x)[ 0*ones(size(x));  0*ones(size(x));  0*ones(size(x))];


% Solving error dynamics PDE under leader-leader communication
[t11_1, e11_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in,f,num_leaders,delta); % Wentzell boundary
[t11_2, e11_2] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, 0, gamma, gamma_in,f,num_leaders,delta); % Not-Wentzell
% Time vector
t = t11_1;
% Total number of agents
n_agents = N + 1;

% Row indices for left and right boundary agents
left_agent_rows  = [1, 2, 3];                                % Agent 1: rows 1–3
right_agent_rows = 3*(n_agents - 1) + [1, 2, 3];             % Agent N+1: rows end-2:end

% Compute log L2 norm at each time step
log_L2_left_wentzell  = log( sqrt( sum( e11_1(left_agent_rows,:).^2, 1 ) ) );
log_L2_right_wentzell = log( sqrt( sum( e11_1(right_agent_rows,:).^2, 1 ) ) );

log_L2_left_nonw      = log( sqrt( sum( e11_2(left_agent_rows,:).^2, 1 ) ) );
log_L2_right_nonw     = log( sqrt( sum( e11_2(right_agent_rows,:).^2, 1 ) ) );

% Plotting
figure;
% Left boundary
subplot(1,2,1);
plot(t, log_L2_left_wentzell, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5); hold on;
plot(t, log_L2_left_nonw, '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\log(\|e(t,0)\|)$', 'Interpreter', 'latex');
title('Left Boundary Agent ');
legend({'$\sigma \neq 0$', '$\sigma = 0$'}, 'Interpreter', 'latex');
grid on;
% Right boundary
subplot(1,2,2);
plot(t, log_L2_right_wentzell, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5); hold on;
plot(t, log_L2_right_nonw, '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\log(\|e(t,1)\|)$', 'Interpreter', 'latex');
title('Right Boundary Agent');
legend({'$\sigma \neq 0$', '$\sigma = 0$'}, 'Interpreter', 'latex');
grid on;