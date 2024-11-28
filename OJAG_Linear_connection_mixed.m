function M = OJAG_Linear_connection_mixed(N, num_leaders)
% This function creates a (N+1)x(N+1) matrix M for linear interpolation.
% Inputs:
% N: Total number of points (N+1 grid points).
% num_leaders: Number of leaders in the domain.

% Step size
h = 1 / N;

% Initialize the matrix M with zeros
M = zeros(N + 1, N + 1);

% Compute leader indices
leader_indices = round(linspace(0, N, num_leaders));

% Iterate over each leader interval
for i = 1:num_leaders - 1
    % Start and end indices for the current interval
    start_idx = leader_indices(i) + 1;    % Convert 0-based to 1-based indexing
    end_idx = leader_indices(i + 1) + 1; % Convert 0-based to 1-based indexing
    
    % Loop through grid points in the current interval
    for j = start_idx:end_idx
        % Current position in the grid
        x = (j - 1) * h;
        
        % Calculate weights for the current interval
        M(j, start_idx) = (leader_indices(i + 1) * h - x) / (leader_indices(i + 1) * h - leader_indices(i) * h);
        M(j, end_idx) = (x - leader_indices(i) * h) / (leader_indices(i + 1) * h - leader_indices(i) * h);
    end
end
end
