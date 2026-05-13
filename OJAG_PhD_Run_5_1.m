% -------------------------------
% Phase Portrait of Agent Trajectories z(t,x)
% -------------------------------

% System Parameters
N = 70; num_leaders = 14; M = 10000; T = 4; kappa = 1; K = 1.4; 
sigma = 6; dx = 1 / N; h = 1 / N; a = 0.001; delta = 1;
f = @(t, z) 0.3 * z;

% Target curve
gamma = @(x) [1.5 * sin(x); 1.5 * cos(x); 6 * ones(size(x))];

% Initial condition (zero offset)
%gamma_in = @(x) [zeros(size(x)); zeros(size(x)); zeros(size(x))];

%gamma_in = @(x) gamma(x) .* (x == 0); % If left boundary agent is on Target Curve
gamma_in = @(x) gamma(x) .* ((x == 0) | (x == 1)); % If both boundary
%agent on target curve
%gamma_in = @(x) repmat(gamma(0), 1, numel(x));  % All agents at gamma(0)

% Solve the PDE for error dynamics
[t12_1, e12_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in, f, num_leaders, delta);

% Time and space info
num_time_steps = size(e12_1, 2);
X = linspace(0, 1, N + 1);  % spatial grid

% Evaluate target and initial curves
gamma_profile = gamma(X);         % 3 x (N+1)
gamma_in_profile = gamma_in(X);   % 3 x (N+1)

% Initialize reconstructed state matrix: z(x,t) = e(x,t) + gamma(x)
Z = zeros(N + 1, 3, num_time_steps);  % (agents) x (components) x (time)

for i = 1:num_time_steps
    e_i = e12_1(:, i);  % Full error vector at time step i
    e_vec = [e_i(1:3:end)'; e_i(2:3:end)'; e_i(3:3:end)'];  % 3 x (N+1)
    Z(:, :, i) = e_vec' + gamma_profile';  % Transpose to (N+1) x 3
end

% ===============================
% 3D Trajectory Plot with Arrows
% ===============================

figure;

% Plot target curve (static)
plot3(gamma_profile(1,:), gamma_profile(2,:), gamma_profile(3,:), ...
      'Color', [0.5, 0, 0], 'LineWidth', 2.5);  % Maroon
hold on;

% Plot boundary agents in orange
plot3(squeeze(Z(1,1,:)), squeeze(Z(1,2,:)), squeeze(Z(1,3,:)), ...
      'Color', [1, 0.5, 0], 'LineWidth', 2);  % Left boundary
plot3(squeeze(Z(end,1,:)), squeeze(Z(end,2,:)), squeeze(Z(end,3,:)), ...
      'Color', [1, 0.5, 0], 'LineWidth', 2);  % Right boundary

% -------------------------------
% Add arrows along trajectories
% -------------------------------
arrow_len = 0.06;
num_arrows = 6;

% Left trajectory
left_traj3 = squeeze(Z(1,:,:));   % 3 x T
add_arrows_arc_length_3D(left_traj3, num_arrows, arrow_len, [0 0 0]);

% Right trajectory
right_traj3 = squeeze(Z(end,:,:)); % 3 x T
add_arrows_arc_length_3D(right_traj3, num_arrows, arrow_len, [0 0 0]);

% -------------------------------
% Axis labels and view
% -------------------------------
xlabel('$z_i^1$', 'Interpreter', 'latex');
ylabel('$z_i^2$', 'Interpreter', 'latex');
hz = zlabel('$z_i^3$', 'Interpreter', 'latex');

% Rotate zlabel
set(hz, 'Rotation', 0);   % 0° keeps it horizontal, 90° vertical

title('3D Phase Portrait with Arc-Length-Based Arrows', 'Interpreter', 'latex');
grid on;
axis equal;
view(3);

% ------------------------------------
% Plotting Phase Portrait (3D)
% ------------------------------------

figure;

% Selected agents to highlight in orange
% orangeColumns = [1, 12, 23, 33, 44, 55, 66, 76, 87, 98, 109, 119, 130, 141]; % for N = 140

%orangeColumns = [1, 4, 8, 10, 13, 16, 19, 23, 26, 29, 32, 35, 38, 41]; % For N = 40

 orangeColumns = [1, 6, 12, 17, 23, 28, 33, 39, 44, 49, 55, 60, 66, 71]; % For N = 70

% Plot target curve (static) – MAROON
plot3(gamma_profile(1,:), gamma_profile(2,:), gamma_profile(3,:), ...
      'Color', [0.5, 0, 0], 'LineWidth', 4);  % maroon
hold on;

% Plot interior agents
for i = 1:N+1
    if ismember(i, orangeColumns)
        % Highlighted ORANGE
        plot3(squeeze(Z(i,1,:)), squeeze(Z(i,2,:)), squeeze(Z(i,3,:)), ...
              'Color', [1, 0.4, 0], 'LineWidth', 3);  % orange solid
    else
        % Black for the rest
        plot3(squeeze(Z(i,1,:)), squeeze(Z(i,2,:)), squeeze(Z(i,3,:)), ...
              'k', 'LineWidth', 2);
    end
end

% Axis labels
xlabel('$z_i^1$', 'Interpreter', 'latex');
ylabel('$z_i^2$', 'Interpreter', 'latex');
hz = zlabel('$z_i^3$', 'Interpreter', 'latex');
set(hz, 'Rotation', 0);

grid on;
axis tight;
view(3);  % 3D perspective

% Export to high-quality vector PDF
print(gcf, 'PhasePortrait3D.pdf', '-dpdf', '-r600');


% ------------------------------------
% Plotting Phase Portrait (2D projection)
% ------------------------------------

figure;

% Extract boundary trajectories
left_traj = squeeze(Z(1, :, :));    % 3 x T
right_traj = squeeze(Z(end, :, :)); % 3 x T

% Plot target curve projection – MAROON
plot(gamma_profile(1,:), gamma_profile(2,:), ...
     'Color', [0.5, 0, 0], 'LineWidth', 3.5);  % maroon
hold on;

% Plot left & right boundaries – ORANGE
plot(left_traj(1,:), left_traj(2,:), ...
     'Color', [1, 0.4, 0], 'LineWidth', 1);  % orange
plot(right_traj(1,:), right_traj(2,:), ...
     'Color', [1, 0.4, 0], 'LineWidth', 1);  % orange

% Add arrows for clarity
arrow_len = 0.06;
num_arrows = 5;
add_arrows_arc_length(left_traj(1:2,:), num_arrows, arrow_len, [0 0 0]);
add_arrows_arc_length(right_traj(1:2,:), num_arrows, arrow_len, [0 0 0]);

% Axis settings
xlabel('$z_i^1$', 'Interpreter', 'latex');
hy = ylabel('$z_i^2$', 'Interpreter', 'latex');
set(hy, 'Rotation', 0);

axis equal;
grid on;

% Export to high-quality vector PDF
print(gcf, 'PhasePortrait2D.pdf', '-dpdf', '-r600');












