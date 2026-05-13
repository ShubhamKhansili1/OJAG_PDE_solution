function draw_arrow_2D(traj2d, indices, len, color)
    for k = indices
        p1 = traj2d(:, k);
        p2 = traj2d(:, k + 1);
        v = p2 - p1;
        if norm(v) > 1e-6
            v = len * v / norm(v);  % Normalize and scale
            quiver(p1(1), p1(2), v(1), v(2), ...
                   'Color', color, ...
                   'LineWidth', 1.2, ...
                   'MaxHeadSize', 1, ...
                   'AutoScale', 'off');
        end
    end
end