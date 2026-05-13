% ----------------------------------------
% Function: Add arrows spaced by arc length
% ----------------------------------------
function add_arrows_arc_length(traj2d, num_arrows, arrow_len, color)
    % Compute arc length along trajectory
    ds = sqrt(sum(diff(traj2d,1,2).^2, 1));
    s = [0, cumsum(ds)];

    % Sample evenly along arc length
    arrow_arc_positions = linspace(0, s(end), num_arrows);
    arrow_indices = zeros(1, num_arrows);
    for i = 1:num_arrows
        [~, idx] = min(abs(s - arrow_arc_positions(i)));
        arrow_indices(i) = idx;
    end

    % Plot arrows
    for k = arrow_indices
        if k < size(traj2d,2)
            p1 = traj2d(:, k);
            p2 = traj2d(:, k+1);
            v = p2 - p1;
            if norm(v) > 1e-6
                v = arrow_len * v / norm(v);
                quiver(p1(1), p1(2), v(1), v(2), ...
                       'Color', color, ...
                       'LineWidth', 1.5, ...
                       'MaxHeadSize', 1.2, ...
                       'AutoScale', 'off');
            end
        end
    end
end