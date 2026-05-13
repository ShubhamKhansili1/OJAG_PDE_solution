function [t, e] = OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, gamma_in, f, num_leaders,delta)
    % Solves a vector-valued parabolic PDE over 1D space
    % e(t,x) in R^3 at each (t,x)

    % Spatial and temporal discretization
  
    dx = 1 / N;
    dt = T / M;
    x = linspace(0, 1, N+1); % Space grid
    t = linspace(0, T, M+1); % Time grid

    % Stability check
    if dt > dx^2 / (2 * a)
        error('Stability condition not satisfied: Reduce dt or increase N.');
    end

    % Build scalar A matrix (size (N+1)x(N+1)) for 1D Laplacian
    A = zeros(N+1, N+1);
    alpha = a / dx^2;
    for i = 2:N
        A(i,i-1) = alpha;
        A(i,i)   = -2 * alpha;
        A(i,i+1) = alpha;
    end
    % Boundary conditions
    A(1,1)     = -sigma/dx - (kappa); % (K+kappa)
    A(1,2)     =  sigma/dx;
    A(N+1,N)   =  sigma/dx;
    A(N+1,N+1) = -sigma/dx - (kappa); %(K+kappa)

    % Expand A and B to act on R^3 vector fields
    A3 = kron(A,eye(3)); % (3(N+1)) x (3(N+1))
    B  = Ram_OJAG_Linear_connection_mixed(N, num_leaders); % size (N+1)x(N+1)
    B3 = kron(B,eye(3)); % (3(N+1)) x (3(N+1))

    % Compute gamma and gamma_in at initial time
    gamma_matrix    = gamma(x);      % 3 x (N+1)
    gamma_in_matrix = gamma_in(x);   % 3 x (N+1)

    gamma_vec    = reshape(gamma_matrix, [], 1);    % (3*(N+1)) x 1
    gamma_in_vec = reshape(gamma_in_matrix, [], 1); % (3*(N+1)) x 1

    % Initialize error as difference between gamma_in and gamma
    e = zeros(3*(N+1), M+1); % (space x vector) × time
    e(:,1) = gamma_in_vec - gamma_vec;

    % Define nonlinear term F (assume f outputs 3 x (N+1))
    F = @(t, gamma_vec, e_vec) reshape( ...
        f(t, reshape(gamma_vec + e_vec, 3, [])) - ...
        f(t, reshape(gamma_vec, 3, [])), [], 1);

    % Time-stepping loop
    for n = 1:M
        A_times_e = A3 * e(:, n);
        B_times_e = B3 * e(:, n);
       
        e(:, n+1) = e(:, n) + dt * (A_times_e + F(t(n), gamma_vec, e(:,n)) - (K+delta) * B_times_e);
    end
end

