function plot_building(footprint, height, color, alpha_val, label)
    numCorners = size(footprint, 1);
    vertices = [footprint, zeros(numCorners, 1); footprint, repmat(height, numCorners, 1)];
    for i = 1:numCorners
        next_i = mod(i, numCorners) + 1;
        wall_face = [i, next_i, next_i+numCorners, i+numCorners];
        patch('Vertices', vertices, 'Faces', wall_face, 'FaceColor', color, 'FaceAlpha', alpha_val, 'EdgeColor', 'k', 'LineWidth', 1);
    end
    roof_face = numCorners+1:2*numCorners;
    patch('Vertices', vertices, 'Faces', roof_face, 'FaceColor', color, 'FaceAlpha', alpha_val, 'EdgeColor', 'k', 'LineWidth', 1);
    center = mean(vertices(roof_face, :));
    text(center(1), center(2), center(3), label, 'HorizontalAlignment', 'center', 'FontSize', 10);
end