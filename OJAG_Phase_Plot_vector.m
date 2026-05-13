function OJAG_Phase_Plot_vector(e, gamma, gamma_in, N)
% OJAG_Phase_Plot_vector
% Plots agent trajectories from vector-valued error solution e.
% Inputs:
%   e         - error matrix, size: 3*(N+1) x T
%   gamma     - target curve function handle (3 x positions)
%   gamma_in  - initial curve function handle (3 x positions)
%   N         - number of spatial divisions (N+1 agents)

% Spatial step and gamma values
h = 1 / N;
x_pos = (0:N) * h;
gammaVal = gamma(x_pos);  % 3 x (N+1)

% Reconstruct agent trajectories: z = e + gamma
x1 = gammaVal(1,:)' + e(1:3:end, :);  % (N+1) x T
x2 = gammaVal(2,:)' + e(2:3:end, :);  % (N+1) x T
x3 = gammaVal(3,:)' + e(3:3:end, :);  % (N+1) x T

% For reference curves
gammaPlot = gamma(0:0.01:1); 
gammainPlot = gamma_in(0:0.01:1);

% Determine how many time steps
[numTimeSteps,~] = size(x1);

% Define selected columns to highlight (e.g., leaders)
highlightCols = [1, 8, 14, 21, 28, 34, numTimeSteps];

% Plot target and initial curves

plot3(gammaPlot(1,:), gammaPlot(2,:), gammaPlot(3,:), ...
      'Color', [1, 0.5, 0], 'LineWidth', 5);  % Target (orange)
hold on
plot3(gammainPlot(1,:), gammainPlot(2,:), gammainPlot(3,:), ...
      'Color', [0, 0.5, 0], 'LineWidth', 2);  % Initial (green)
hold on

% Highlight selected agent trajectories
for i = highlightCols
    plot3(x1(i,:), x2(i,:), x3(i,:), 'g', 'LineWidth', 2);
end

% Plot remaining agent trajectories in black
for i = 1:numTimeSteps
    if ~ismember(i,highlightCols)
        plot3(x1(i,:), x2(i,:), x3(i,:), 'k');
    end
end
hold off;
% Final touches
xlabel('$z_i^1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$z_i^2$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$z_i^3$', 'Interpreter', 'latex', 'FontSize', 20);
title('Agent Phase Portrait (Vector Format)', 'Interpreter', 'latex');
grid on;
legend({'Target Curve $\gamma$', 'Initial Curve $\gamma_{in}$', ...
        'Highlighted Agents', 'Other Agents'}, ...
        'Interpreter', 'latex', 'Location', 'best');
hold off;
end

