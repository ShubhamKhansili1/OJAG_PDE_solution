function B = Ram_OJAG_Linear_connection_mixed(N, num_leaders)
% Constructs a piecewise linear approximation matrix B ∈ R^{(N+1)×(N+1)}
% Each node's value is linearly interpolated from its two nearest leaders

    % Grid and leader indices
    num_nodes = N + 1;
    x = linspace(0, 1, num_nodes);
    leader_indices = round(linspace(1, num_nodes, num_leaders));
    leader_positions = x(leader_indices);

    % Initialize matrix
    B = zeros(num_nodes);

    % Loop over each node in the domain
    for i = 1:num_nodes
        xi = x(i);

        % If i is a leader, set identity
        if ismember(i, leader_indices)
            B(i, i) = 1;
        else
            % Find closest two leaders surrounding this node
            idx_left = find(leader_positions <= xi, 1, 'last');
            idx_right = find(leader_positions >= xi, 1, 'first');

            % If only one leader is found (at boundary), use constant
            if isempty(idx_left)
                B(i, leader_indices(idx_right)) = 1;
            elseif isempty(idx_right)
                B(i, leader_indices(idx_left)) = 1;
            elseif idx_left == idx_right
                B(i, leader_indices(idx_left)) = 1;
            else
                % Perform linear interpolation
                xL = leader_positions(idx_left);
                xR = leader_positions(idx_right);
                wL = (xR - xi) / (xR - xL);
                wR = (xi - xL) / (xR - xL);
                B(i, leader_indices(idx_left)) = wL;
                B(i, leader_indices(idx_right)) = wR;
            end
        end
    end
end
