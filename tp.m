% Define the grid
N = 100;
x = linspace(0, 1, N+1);              % 101 grid points
f_true = x;            % Example true function

% Select 33 uniformly spaced known points from x
x_known = linspace(0, 1, 33);         % 33 uniform known locations
f_known =  x_known;      % Evaluate the true function at known points

% Piecewise constant approximation
f_pc = zeros(size(x));
for i = 1:length(x)
    idx = find(x_known <= x(i), 1, 'last');
    if isempty(idx)
        idx = 1;
    end
    f_pc(i) = f_known(idx);
end

% Piecewise linear interpolation
f_pl = interp1(x_known, f_known, x, 'linear', 'extrap');

% Plot
figure;
plot(x, f_true, 'k', 'LineWidth', 2); hold on;
plot(x, f_pc, '--b', 'LineWidth', 1.5);
plot(x, f_pl, '-.r', 'LineWidth', 1.5);
plot(x_known, f_known, 'ko', 'MarkerFaceColor', 'g');
legend('True Function', 'Piecewise Constant', 'Piecewise Linear', 'Known Points');
xlabel('x'); ylabel('f(x)');
title('Piecewise Constant vs Linear Approximation (33 Uniform Points)');
grid on;

% Compute and display error norms
err_pc = norm(f_true - f_pc, 2);     % L2 norm
err_pl = norm(f_true - f_pl, 2);     % L2 norm
fprintf('L2 Error (Piecewise Constant): %.4f\n', err_pc);
fprintf('L2 Error (Piecewise Linear):   %.4f\n', err_pl);
