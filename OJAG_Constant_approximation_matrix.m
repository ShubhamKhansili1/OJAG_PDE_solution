function M = OJAG_Constant_approximation_matrix(N, num_leaders)
% This function creates a (N+1)x(N+1) matrix M for constant approximation
% with split influence in each interval.
% Inputs:
% N: Total number of points (N+1 grid points).
% num_leaders: Number of leaders in the domain.



% Initialize the matrix M with zeros
M = zeros(N + 1, N + 1);

% Compute leader indices
leader_indices = round(linspace(0, N, num_leaders));

% Iterate over each leader interval
for i = 1:num_leaders - 1
    % Start and end indices for the current interval
    start_idx = leader_indices(i) + 1;    % Convert 0-based to 1-based indexing
    end_idx = leader_indices(i + 1) + 1; % Convert 0-based to 1-based indexing
    
    % Split the interval into two equal halves
    midpoint_idx = floor((start_idx + end_idx) / 2);
    
    % Assign the first half to the starting leader
    M(start_idx:midpoint_idx, start_idx) = 1;
    
    % Assign the second half to the ending leader
    M(midpoint_idx+1:end_idx, end_idx) = 1;
end

end
