% ----------------------------------------
% Local Function: Draw arrows using larger time step gaps
% ----------------------------------------
function draw_arrow_set(traj, indices, len, color)
    for k = indices
        p1 = traj(:, k);
        p2 = traj(:, k + 1);
        v = p2 - p1;
        if norm(v) > 1e-6
            v = len * v / norm(v);  % Normalize and scale
            quiver3(p1(1), p1(2), p1(3), ...
                    v(1), v(2), v(3), ...
                    'Color', color, ...
                    'LineWidth', 1.2, ...
                    'MaxHeadSize', 0.8, ...
                    'AutoScale', 'off');
        end
    end
end