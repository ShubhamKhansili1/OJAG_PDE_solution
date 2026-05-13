% This code is for the Phase Portrait

% System Parameter's
N = 70; num_leaders = 7; M = 10000; T = 4; kappa = 1; K = 1.4; 
sigma= 6; dx = 1 / N; h = 1/N; a = 0.004; delta =1 ; f = @(t, z) 0.3* z; 
% --- Target curve -------------------------------------------
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];   
%gamma = @(x) [1.5*sin(x); 1.5*sin(x);1.5*sin(x)];  

%gamma_in = @(x) gamma(x) .* (x == 0); 
%gamma_in = @(x) gamma(x) .* ((x == 0) | (x == 1));
%gamma_in = @(x) repmat(gamma(0), 1, numel(x));  % All agents at gamma(0)
gamma_in = @(x)[ 0*ones(size(x));  0*ones(size(x));  0*ones(size(x))];
%gamma = @(x) [h*cos(pi*N*x/3); h*sin(pi*N*x/3); 0*ones(size(x))]; 


% Solving PDe error dynamics
[t12_1, e12_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
%  Plot actual agent states z = e + gamma
num_time_steps = size(e12_1, 2);
Z = zeros(N+1, 3, num_time_steps); % Agent state z(t,x)
% Compute gamma(x) at spatial nodes
X = linspace(0, 1, N+1);
gamma_profile = gamma(X);  % 3 x (N+1)
gamma_in_profile = gamma_in(X);  % 3 x (N+1)

% Reconstruct agent states z(x,t) = e(x,t) + gamma(x)
for i = 1:num_time_steps
    e_i = e12_1(:, i);  % Full error vector at time i
    e_vec = [e_i(1:3:end)'; e_i(2:3:end)'; e_i(3:3:end)'];  % 3 x (N+1)
    Z(:,:,i) = e_vec' + gamma_profile';  % (N+1) x 3
end

% Plotting Phase Portrait
figure; 

% Define the indices of the columns to be plotted in orange
yellowColumns = [1, 13,24,36,48,59,71];

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
