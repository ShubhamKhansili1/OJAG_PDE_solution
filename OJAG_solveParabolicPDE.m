function [e, t, x]= OJAG_solveParabolicPDE(N, M, T, a, kappa, K, sigma, gamma, f,num_leaders)
    % solveParabolicPDE solves a parabolic PDE with nonlinear F.
    %
    % INPUTS:
    % N      : Number of spatial divisions
    % M      : Number of time divisions
    % T      : Total time
    % a      : Diffusion coefficient
    % kappa  : Coefficient in boundary condition
    % K      : Coefficient in boundary condition
    % sigma  : Coefficient in boundary condition
    % gamma  : Function handle for initial condition gamma(x)
    % f      : Function handle for nonlinearity f(t, x)
    
    % Spatial and temporal grid
    L = 1; % Length of spatial domain
    dx = L / N; % Spatial step
    dt = T / M; % Time step
    x = linspace(0, L, N+1); % Spatial grid
    t = linspace(0, T, M+1); % Time grid
    
    % Stability condition check
    if dt > dx^2 / (2 * a)
        error('Stability condition not satisfied: Reduce dt or increase N.');
    end
    
    % Initialize solution matrix
    e = zeros(N+1, M+1); % Solution matrix (rows: spatial, cols: time)
    
    % Initial condition: e(x, 0) = -gamma(x)
    e(:, 1) = -gamma(x); % Apply gamma(x) for initial condition
    
    % Construct full matrix A
    A = zeros(N+1, N+1);
    
    % Fill interior rows for the second derivative
    alpha = a / dx^2; % Coefficient for second derivative
    for i = 2:N
        A(i, i-1) = alpha; % Lower diagonal
        A(i, i) = -2 * alpha; % Main diagonal
        A(i, i+1) = alpha; % Upper diagonal
    end
    
    % Boundary conditions in matrix form
    % Left boundary
    A(1, 1) = -sigma / dx - (K+kappa);
    A(1, 2) = sigma / dx;
    
    % Right boundary
    A(N+1, N) = sigma / dx;
    A(N+1, N+1) = -sigma / dx - (K+kappa );
    
    % The control i.e the linear connection is stored in the B matrix_linear
      B = OJAG_Linear_connection_mixed(N, num_leaders);
      gamma_vec = gamma(x); % Size (N+1, 1)
      F=@(t,gamma,e) f(t,gamma+e)-f(t,gamma); 
      % Time-stepping loop
      for n = 1:M
    
        % Ensure matrix-vector multiplication for A * e(:, n)
         A_times_e = A * e(:, n); % Size is (N+1, 1)
         B_times_e = B * e(:, n);

        % Update rule: e^{n+1} = e^n + dt * (A * e^n + F_n)
        e(:, n+1) = e(:, n) + dt * (A_times_e +F(t,gamma_vec(:),e(:,n)) -K*B_times_e);
      end 
end

