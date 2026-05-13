function add_arrows_arc_length_3D(traj, num_arrows, arrow_len, color)
    % traj: 3 x N
    % num_arrows: number of arrows to plot
    % arrow_len: scaling of arrows
    % color: RGB

    % Compute arc length
    diffs = diff(traj,1,2);
    seg_lengths = sqrt(sum(diffs.^2,1));
    arc_length = [0, cumsum(seg_lengths)];

    % Interpolate equally spaced positions
    total_len = arc_length(end);
    s_query = linspace(0, total_len, num_arrows+2); % avoid endpoints
    s_query = s_query(2:end-1);

    % Arrow positions
    arrow_pos = zeros(3, num_arrows);
    arrow_dir = zeros(3, num_arrows);

    for k = 1:num_arrows
        % Find segment index
        idx = find(arc_length >= s_query(k), 1);
        if idx > 1
            t = (s_query(k) - arc_length(idx-1)) / (arc_length(idx) - arc_length(idx-1));
            arrow_pos(:,k) = (1-t)*traj(:,idx-1) + t*traj(:,idx);
            arrow_dir(:,k) = diffs(:,idx-1);
        end
    end

    % Normalize directions
    arrow_dir = arrow_len * arrow_dir ./ vecnorm(arrow_dir);

    % Plot arrows using quiver3
    quiver3(arrow_pos(1,:), arrow_pos(2,:), arrow_pos(3,:), ...
            arrow_dir(1,:), arrow_dir(2,:), arrow_dir(3,:), ...
            0, 'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 1.2);
end
