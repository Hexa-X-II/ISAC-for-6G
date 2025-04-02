function valid_walls = plot_intersection_with_reflections(params, FoV_size, alpha)
% plot_intersection_with_reflections
%   Plots a 3D urban intersection visualization showing the intersection,
%   buildings, Tx and Rx positions, and for each wall in params.walls computes
%   and overlays the expected reflection point and corresponding ray.
%   Each reflection is annotated with its wall number.
%
%   The function uses a tolerance (tol_bound) when checking if the computed
%   reflection point lies within the wall segment boundaries. This is useful
%   for walls that are degenerate in one coordinate (e.g. vertical walls with
%   constant z) so that slight deviations (e.g. due to numerical effects) are
%   accepted.
%
%   If the computed reflection point is nearly equal to the receiver, it is
%   skipped.
%
% Inputs:
%   params   - simulation parameters structure containing:
%                pTx: [3x1] transmitter position,
%                pRx: [3x1] receiver position,
%                walls: structure array with fields:
%                         p     : representative point on the wall,
%                         n     : unit normal vector,
%                         start : starting point of the wall segment,
%                         end   : ending point of the wall segment.
%   FoV_size - [width, depth, building height] of the intersection.
%   alpha    - building boundary fraction (e.g., 0.3).
%
% Output:
%   valid_walls - vector of wall indices for which the computed reflection
%                 point is valid.
%
% Example usage:
%   valid_walls = plot_intersection_with_reflections(params, FoV_size, alpha);
%   disp('Valid reflection walls:');
%   disp(valid_walls);

% Tolerance for checking boundaries (in meters)
tol_bound = 0.2;

% Create figure
figure('Name', '3D Urban Intersection with Reflections', 'Position', [100, 100, 800, 600]);
hold on;
grid on;
axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('3D Urban Intersection Visualization with Reflection Points');

%% Plot Intersection Boundaries (Road Edges)
plot3([0 FoV_size(1)], [0 0], [0 0], 'k-', 'LineWidth', 2);       % Bottom edge
plot3([FoV_size(1) FoV_size(1)], [0 FoV_size(2)], [0 0], 'k-', 'LineWidth', 2);  % Right edge
plot3([FoV_size(1) 0], [FoV_size(2) FoV_size(2)], [0 0], 'k-', 'LineWidth', 2);  % Top edge
plot3([0 0], [FoV_size(2) 0], [0 0], 'k-', 'LineWidth', 2);       % Left edge

%% Plot Building Patches (Assuming Buildings Occupy the Four Corners)
ix_start = FoV_size(1) * alpha;
ix_end   = FoV_size(1) * (1 - alpha);
iy_start = FoV_size(2) * alpha;
iy_end   = FoV_size(2) * (1 - alpha);
building_height = FoV_size(3);

draw_building_patch([0, ix_start], [0, iy_start], building_height);
draw_building_patch([ix_end, FoV_size(1)], [0, iy_start], building_height);
draw_building_patch([0, ix_start], [iy_end, FoV_size(2)], building_height);
draw_building_patch([ix_end, FoV_size(1)], [iy_end, FoV_size(2)], building_height);

%% Plot Transmitter and Receiver
pTx = params.pTx;
pRx = params.pRx;
plot3(pTx(1), pTx(2), pTx(3), 'b^', 'MarkerSize', 10, 'LineWidth', 2); % Tx (blue triangle)
plot3(pRx(1), pRx(2), pRx(3), 'gs', 'MarkerSize', 10, 'LineWidth', 2); % Rx (green square)

%% Process All Walls and Compute Reflection Points
valid_walls = [];  % initialize list of valid wall indices
tol = 1e-6;        % tolerance for comparing points

for idx = 1:length(params.walls)
    wall = params.walls(idx);
    
    % Use wall.p as the representative point on the wall.
    p0 = wall.p;
    
    % Compute the wall's boundaries from its start and end points.
    bounds = [min(wall.start(1), wall.end(1)), max(wall.start(1), wall.end(1));
              min(wall.start(2), wall.end(2)), max(wall.start(2), wall.end(2));
              min(wall.start(3), wall.end(3)), max(wall.start(3), wall.end(3))];
    
    % Compute the mirror image of Tx with respect to the wall plane.
    pTx_mirror = pTx - 2 * dot(pTx - p0, wall.n) * wall.n;
    
    % Parameterize the line from pTx_mirror to Rx:
    % p(t) = pTx_mirror + t*(pRx - pTx_mirror)
    denom = dot(pRx - pTx_mirror, wall.n);
    if abs(denom) > eps
        t = dot(p0 - pTx_mirror, wall.n) / denom;
    else
        t = 0;
    end
    
    % Compute the reflection point on the wall plane.
    p_reflect = pTx_mirror + t * (pRx - pTx_mirror);
    
    % If the computed reflection point is (nearly) equal to the receiver, skip this wall.
    if norm(p_reflect - pRx) < tol
        fprintf('Wall %d: Reflection point equals receiver, skipping.\n', idx);
        continue;
    end
    
    % Check whether the reflection point falls within the wall segment bounds,
    % allowing a tolerance (tol_bound).
    isValid = (p_reflect(1) >= bounds(1,1)-tol_bound) && (p_reflect(1) <= bounds(1,2)+tol_bound) && ...
              (p_reflect(2) >= bounds(2,1)-tol_bound) && (p_reflect(2) <= bounds(2,2)+tol_bound) && ...
              (p_reflect(3) >= bounds(3,1)-tol_bound) && (p_reflect(3) <= bounds(3,2)+tol_bound);
    
    if isValid
        markerStyle = 'ko';  % valid reflection: black circle
        lineColor   = 'k';
        valid_walls(end+1) = idx; %#ok<AGROW>
    else
        markerStyle = 'ro';  % invalid (extrapolated): red circle
        lineColor   = 'r';
    end
    
    % Plot the reflection point.
    plot3(p_reflect(1), p_reflect(2), p_reflect(3), markerStyle, 'MarkerSize', 8, 'LineWidth', 2);
    
    % Plot the reflecting ray segments.
    line([pTx(1) p_reflect(1)], [pTx(2) p_reflect(2)], [pTx(3) p_reflect(3)], ...
         'LineStyle', '--', 'Color', lineColor, 'LineWidth', 2);
    line([p_reflect(1) pRx(1)], [p_reflect(2) pRx(2)], [p_reflect(3) pRx(3)], ...
         'LineStyle', '--', 'Color', lineColor, 'LineWidth', 2);
    
    % Annotate the reflection point with its wall number.
    text(p_reflect(1), p_reflect(2), p_reflect(3), sprintf(' Wall %d', idx), ...
         'FontSize', 10, 'Color', lineColor);
end

view(3);
hold off;
end

%% Helper Function: Draw a 3D Building Patch
function draw_building_patch(x_range, y_range, height)
% draw_building_patch draws a 3D rectangular building patch as a semi-transparent prism.
%
% Inputs:
%   x_range - two-element vector [xmin, xmax]
%   y_range - two-element vector [ymin, ymax]
%   height  - height of the building

    x = [x_range(1), x_range(2)];
    y = [y_range(1), y_range(2)];
    z = [0, height];

    [X, Y, Z] = ndgrid(x, y, z);
    vertices = [X(:), Y(:), Z(:)];

    faces = [...
        1 2 4 3;  % bottom face
        5 6 8 7;  % top face
        1 2 6 5;  % side face
        2 4 8 6;  % side face
        4 3 7 8;  % side face
        3 1 5 7]; % side face

    patch('Vertices', vertices, 'Faces', faces, ...
          'FaceColor', 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
