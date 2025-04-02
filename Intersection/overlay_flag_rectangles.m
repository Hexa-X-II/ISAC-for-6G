function overlay_flag_rectangles(x_points, y_points, flag_grid)
    grid_size = size(flag_grid, 1);
    dx = x_points(2) - x_points(1);
    dy = y_points(2) - y_points(1);
    for i = 1:grid_size
        for j = 1:grid_size
            if flag_grid(i,j) == 1
                rectangle('Position', [x_points(j)-dx/2, y_points(i)-dy/2, dx, dy], 'FaceColor', 'w', 'EdgeColor', 'none');
            elseif flag_grid(i,j) == 2
                rectangle('Position', [x_points(j)-dx/2, y_points(i)-dy/2, dx, dy], 'FaceColor', 'k', 'EdgeColor', 'none');
            end
        end
    end
end