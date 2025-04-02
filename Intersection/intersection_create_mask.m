function mask = intersection_create_mask(x_points, y_points, tx_pos, intersection_size, alpha)
    num_x = length(x_points);
    num_y = length(y_points);
    mask = ones(num_y, num_x);
    
    % Compute building footprints in absolute coordinates.
    % Bottom-left building: [0, intersection_size(1)*alpha] x [0, intersection_size(2)*alpha]
    bl_poly_x = [0, intersection_size(1)*alpha, intersection_size(1)*alpha, 0];
    bl_poly_y = [0, 0, intersection_size(2)*alpha, intersection_size(2)*alpha];
    
    % Bottom-right building: [intersection_size(1)*(1-alpha), intersection_size(1)] x [0, intersection_size(2)*alpha]
    br_poly_x = [intersection_size(1)*(1-alpha), intersection_size(1), intersection_size(1), intersection_size(1)*(1-alpha)];
    br_poly_y = [0, 0, intersection_size(2)*alpha, intersection_size(2)*alpha];
    
    % Top-left building: [0, intersection_size(1)*alpha] x [intersection_size(2)*(1-alpha), intersection_size(2)]
    tl_poly_x = [0, intersection_size(1)*alpha, intersection_size(1)*alpha, 0];
    tl_poly_y = [intersection_size(2)*(1-alpha), intersection_size(2)*(1-alpha), intersection_size(2), intersection_size(2)];
    
    % Top-right building: [intersection_size(1)*(1-alpha), intersection_size(1)] x [intersection_size(2)*(1-alpha), intersection_size(2)]
    tr_poly_x = [intersection_size(1)*(1-alpha), intersection_size(1), intersection_size(1), intersection_size(1)*(1-alpha)];
    tr_poly_y = [intersection_size(2)*(1-alpha), intersection_size(2)*(1-alpha), intersection_size(2), intersection_size(2)];
    
    % Store the building polygons in a cell array
    building_polys = {struct('x', bl_poly_x, 'y', bl_poly_y), ...
                      struct('x', br_poly_x, 'y', br_poly_y), ...
                      struct('x', tl_poly_x, 'y', tl_poly_y), ...
                      struct('x', tr_poly_x, 'y', tr_poly_y)};
    
    % Loop over each grid point (absolute coordinates)
    for j = 1:num_x
        for i = 1:num_y
            p = [x_points(j); y_points(i)];
            blocked = false;
            for k = 1:length(building_polys)
                poly = building_polys{k};
                % If the grid point lies inside a building footprint, mark as blocked.
                if inpolygon(p(1), p(2), poly.x, poly.y)
                    blocked = true;
                    break;
                end
                % If the line segment from tx_pos to p intersects any building edge, mark as blocked.
                if check_line_segment_intersection(tx_pos, p, poly.x, poly.y)
                    blocked = true;
                    break;
                end
            end
            mask(i, j) = ~blocked;
        end
    end
    
    mask = logical(mask);
end


function intersect = check_line_segment_intersection(p1, p2, poly_x, poly_y)
    % Ensure we use only the first two elements (x,y)
    p1 = p1(1:2);
    p2 = p2(1:2);
    intersect = false;
    num_vertices = length(poly_x);
    for k = 1:num_vertices
        k_next = mod(k, num_vertices) + 1;
        if lineseg_intersect(p1, p2, [poly_x(k); poly_y(k)], [poly_x(k_next); poly_y(k_next)])
            intersect = true;
            return;
        end
    end
end

function flag = lineseg_intersect(p1, p2, p3, p4)
    % 2D line segment intersection test.
    d1 = p2 - p1;
    d2 = p4 - p3;
    denominator = d1(1)*d2(2) - d1(2)*d2(1);
    if abs(denominator) < 1e-10
        flag = false; % segments are parallel or collinear
        return;
    end
    diff = p3 - p1;
    t = (diff(1)*d2(2) - diff(2)*d2(1)) / denominator;
    u = (diff(1)*d1(2) - diff(2)*d1(1)) / denominator;
    flag = (t >= 0 && t <= 1 && u >= 0 && u <= 1);
end
