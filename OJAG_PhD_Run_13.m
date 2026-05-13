% System Parameter's
N = 70; num_leaders = 7; M = 10000; T = 4; kappa = 1; K = 1.4; 
sigma= 6; dx = 1 / N; h = 1/N; a = 0.001; delta =1 ; f = @(t, z) 0.3* z; 
% --- Target curve -------------------------------------------
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];   
gamma_in = @(x)[ 0*ones(size(x));  0*ones(size(x));  0*ones(size(x))];

% Solving error PDE system
[t13_1, e13_1] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in,f,num_leaders,delta);% Boundary intact + Comm. Leader
[t13_2, e13_2] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, 0, gamma,gamma_in, f, num_leaders,delta);% Independent Bo. + NoComm. Leader

% Prepare space and time grids
x = linspace(0, 1, N+1);
t = t13_1;
[X, T] = meshgrid(x, t);

% ===================================
% ======= Wentzell Case =====
% ===================================

% Extract components (time x space)
e1_w = e13_1(1:3:end, :)';
e2_w = e13_1(2:3:end, :)';
e3_w = e13_1(3:3:end, :)';

figure('Name', 'Wentzell Boundary Condition', 'NumberTitle', 'off');

% e1 subplot
subplot(1,3,1);
surf(T, X, e1_w, 'EdgeColor', 'none');
nColors = 256; % Number of colors
redChannel = ones(nColors, 1); % Full red
greenChannel = linspace(0, 0.7, nColors)'; % Gradually increase green towards yellow
blueChannel = zeros(nColors, 1); % No blue
colormap([redChannel, greenChannel, blueChannel]);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$x$', 'Interpreter', 'latex');
zlabel('$e^1(t,x)$', 'Interpreter', 'latex');
title('3D Surface Plot: communicating leaders');
%colorbar;
view(3);
% e2 subplot
subplot(1,3,2);
surf(T, X, e2_w, 'EdgeColor', 'none');
nColors = 256; % Number of colors
redChannel = ones(nColors, 1); % Full red
greenChannel = linspace(0, 0.7, nColors)'; % Gradually increase green towards yellow
blueChannel = zeros(nColors, 1); % No blue
colormap([redChannel, greenChannel, blueChannel]);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$x$', 'Interpreter', 'latex');
zlabel('$e^2(t,x)$', 'Interpreter', 'latex');
title('3D Surface Plot: communicating leaders');
%colorbar;
view(3);

% e3 subplot
subplot(1,3,3);
surf(T, X, e3_w, 'EdgeColor', 'none');
nColors = 256; % Number of colors
redChannel = ones(nColors, 1); % Full red
greenChannel = linspace(0, 0.7, nColors)'; % Gradually increase green towards yellow
blueChannel = zeros(nColors, 1); % No blue
colormap([redChannel, greenChannel, blueChannel]);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$x$', 'Interpreter', 'latex');
zlabel('$e^3(t,x)$', 'Interpreter', 'latex');
title('3D Surface Plot: communicating leaders');
%colorbar;
view(3);


% ===================================
% ==== Non-Wentzell Case  ====
% ===================================

% Extract components (time x space)
e1_nw = e13_2(1:3:end, :)';
e2_nw = e13_2(2:3:end, :)';
e3_nw = e13_2(3:3:end, :)';

figure('Name', 'Non-Wentzell Boundary Condition', 'NumberTitle', 'off');

% e1 subplot
subplot(1,3,1);
surf(T, X, e1_nw, 'EdgeColor', 'none');
nColors = 256; % Number of colors
redChannel = ones(nColors, 1); % Full red
greenChannel = linspace(0, 0.7, nColors)'; % Gradually increase green towards yellow
blueChannel = zeros(nColors, 1); % No blue
colormap([redChannel, greenChannel, blueChannel]);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$x$', 'Interpreter', 'latex');
zlabel('$e^1(t,x)$', 'Interpreter', 'latex');
title('3D Surface Plot: Non-communicating leaders');
%colorbar;
view(3);
% e2 subplot
subplot(1,3,2);
surf(T, X, e2_nw, 'EdgeColor', 'none');
nColors = 256; % Number of colors
redChannel = ones(nColors, 1); % Full red
greenChannel = linspace(0, 0.7, nColors)'; % Gradually increase green towards yellow
blueChannel = zeros(nColors, 1); % No blue
colormap([redChannel, greenChannel, blueChannel]);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$x$', 'Interpreter', 'latex');
zlabel('$e^2(t,x)$', 'Interpreter', 'latex');
title('3D Surface Plot: Non-communicating leaders');
%colorbar;
view(3);

% e3 subplot
subplot(1,3,3);
surf(T, X, e3_nw, 'EdgeColor', 'none');
nColors = 256; % Number of colors
redChannel = ones(nColors, 1); % Full red
greenChannel = linspace(0, 0.7, nColors)'; % Gradually increase green towards yellow
blueChannel = zeros(nColors, 1); % No blue
colormap([redChannel, greenChannel, blueChannel]);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$x$', 'Interpreter', 'latex');
zlabel('$e^3(t,x)$', 'Interpreter', 'latex');
title('3D Surface Plot: Non-communicating leaders');
%colorbar;
view(3);




