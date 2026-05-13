function [t,e] = OJAG_solveParabolicPDE_constant_vec(N, M, T, a, kappa, K, sigma, gamma,gamma_in, f, num_leaders,delta)
    % Solves a vector-valued parabolic PDE with constant control
    % State: e(t,x) ∈ R^3, space is 1D, time evolves over [0, T]

    % Spatial and temporal discretization
    L = 1;
    dx = L / N;
    dt = T / M;
    x = linspace(0, L, N+1);
    t = linspace(0, T, M+1);

    % Stability check
    if dt > dx^2 / (2 * a)
        error('Stability condition not satisfied: Reduce dt or increase N.');
    end

    % Build scalar Laplacian A (for 1D)
    A = zeros(N+1, N+1);
    alpha = a / dx^2;
    for i = 2:N
        A(i,i-1) = alpha;
        A(i,i)   = -2*alpha;
        A(i,i+1) = alpha;
    end
    % Boundary conditions
    A(1,1) = -sigma/dx - ( kappa); % (K + kappa);
    A(1,2) = sigma/dx;
    A(N+1,N) = sigma/dx;
    A(N+1,N+1) = -sigma/dx - ( kappa); % (K + kappa);

    % Expand A and C to act on R^3-valued vectors
    A3 = kron(A,eye(3));  % (3N+3) × (3N+3)
   % C = OJAG_Constant_approximation_matrix(N, num_leaders);  % (N+1)x(N+1)
     C = Ram_OJAG_Constant_approximation_matrix(N, num_leaders);  % (N+1)x(N+1)
    C3 = kron(C,eye(3));  % (3N+3) × (3N+3)
 
    % Compute gamma and gamma_in at initial time
    gamma_matrix    = gamma(x);      % 3 x (N+1)
    gamma_in_matrix = gamma_in(x);   % 3 x (N+1)

    gamma_vec    = reshape(gamma_matrix, [], 1);    % (3*(N+1)) x 1
    gamma_in_vec = reshape(gamma_in_matrix, [], 1); % (3*(N+1)) x 1

    % Initialize error as difference between gamma_in and gamma
    e = zeros(3*(N+1), M+1); % (space x vector) × time
    e(:,1) = gamma_in_vec - gamma_vec;

% Nonlinear function
    F = @(t, gamma_vec, e_vec) reshape( ...
        f(t, reshape(gamma_vec + e_vec, 3, [])) - ...
        f(t, reshape(gamma_vec, 3, [])), [], 1);

    % Time-stepping loop
    for n = 1:M
        A_term = A3 * e(:,n);
        C_term = C3 * e(:,n);
     
        e(:,n+1) = e(:,n) + dt * (A_term + F(t(n), gamma_vec, e(:,n)) - (K+delta)* C_term );
    end
end
