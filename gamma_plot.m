% Define the parameter range
x = linspace(0, 1, 100);

% Define the parametric curve gamma
gamma = @(x) [5 * sin(x); 1.5 * sin(2*x); 6 *cos(2*x)];

% Evaluate gamma over the range of x
gamma_vals = gamma(x);

% Plot the curve in 3D
figure;
plot3(gamma_vals(1, :), gamma_vals(2, :), gamma_vals(3, :), 'LineWidth', 2);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Plot of \gamma(x)');
grid on;
axis equal;
