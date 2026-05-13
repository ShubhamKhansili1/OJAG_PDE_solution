% Code for Chapter 2 in thesis to show, the system is unstable without
% control. 
% System Parameter's
N = 70; num_leaders = 0; M = 10000; T = 4; kappa = 0; K = 0; 
sigma= 0; dx = 1 / N; h = 1/N; a = 0; delta =0 ; f = @(t, z) 0.5* sin(z); 
% --- Target curve -------------------------------------------
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];   
gamma_in = @(x)[ 0*ones(size(x));  0*ones(size(x));  0*ones(size(x))];
%Solving the error PDE-system

[t15_1, e15_1] = OJAG_chapter2_PDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader

L2_norm15_1 = compute_vec_L2(e15_1,dx);  

% Time and space info
num_time_steps = size(e15_1, 2);
X = linspace(0, 1, N + 1);  % spatial grid

% Evaluate target and initial curves
gamma_profile = gamma(X);         % 3 x (N+1)
gamma_in_profile = gamma_in(X);   % 3 x (N+1)

% Initialize reconstructed state matrix: z(x,t) = e(x,t) + gamma(x)
Z = zeros(N + 1, 3, num_time_steps);  % (agents) x (components) x (time)

for i = 1:num_time_steps
    e_i = e15_1(:, i);  % Full error vector at time step i
    e_vec = [e_i(1:3:end)'; e_i(2:3:end)'; e_i(3:3:end)'];  % 3 x (N+1)
    Z(:, :, i) = e_vec' + gamma_profile';  % (N+1) x 3
end

% -----------------------------
% Subplot: L2 norm + Phase Portrait
% -----------------------------
figure('Name', 'Error Norm and Phase Portrait', 'NumberTitle', 'off');

% ---- Subplot 1: L2 Norm vs t ----
subplot(1,2,1);
plot(t15_1, L2_norm15_1, 'Color', [0.7, 0.3, 0.3], 'LineWidth', 2);
xlabel('$$ t $$', 'Interpreter', 'latex');
ylabel('$$\|e(t)\|$$', 'Interpreter', 'latex');
title('$\|e(t, \cdot)\|$ vs. $t$: $a = 0.004$', 'Interpreter', 'latex');
grid on;

% ---- Subplot 2: Phase Portrait ----
subplot(1,2,2);
hold on;
% Target curve
plot3(gamma_profile(1,:), gamma_profile(2,:), gamma_profile(3,:), ...
      'Color', [0.8, 0.6, 0], 'LineWidth', 2.5);  % Dark yellow

% Agent trajectories
for i = 1:N+1
    plot3(squeeze(Z(i,1,:)), squeeze(Z(i,2,:)), squeeze(Z(i,3,:)), ...
        'Color', [0.7, 0.3, 0.3], 'LineWidth', 1);  % Maroon-like
end

xlabel('$z_i^1$', 'Interpreter', 'latex');
ylabel('$z_i^2$', 'Interpreter', 'latex');
zlabel('$z_i^3$', 'Interpreter', 'latex');
title('Phase Portrait: Agent Trajectories $z_i(t)$', 'Interpreter', 'latex');
grid on;
axis tight;
view(3);





