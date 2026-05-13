function C = Ram_OJAG_Constant_approximation_matrix(N, num_leaders)
% Constructs a constant approximation matrix C ∈ R^{(N+1)x(N+1)}
% Each row of C selects the nearest leader's value
% Used in e_t = ... - K * C * e

    % Total number of spatial nodes
    num_nodes = N + 1;

    % Initialize C as a zero matrix
    C = zeros(num_nodes);

    % Uniformly place leaders among the grid points
    leader_indices = round(linspace(1, num_nodes, num_leaders));

    % Assign each node to its nearest leader
    for i = 1:num_nodes
        % Find nearest leader index
        [~, idx] = min(abs(i - leader_indices));
        leader = leader_indices(idx);

        % Set row i to select value at leader index
        C(i, leader) = 1;
    end
end
