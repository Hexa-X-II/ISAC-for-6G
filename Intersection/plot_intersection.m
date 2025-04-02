function plot_intersection(params, intersection_size, alpha)
    figure('Color', 'white', 'Position', [100 100 1000 800]);
    ax = gca;
    ax.Box = 'on';
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.3;
    hold on;
    
    building_color = [0.8, 0.2, 0.2];
    ground_color = [0.9, 0.9, 0.9];
    road_color = [0.7, 0.7, 0.7];
    wall_alpha = 0.6;
    
    ix_start = intersection_size(1) * alpha;
    ix_end   = intersection_size(1) * (1 - alpha);
    iy_start = intersection_size(2) * alpha;
    iy_end   = intersection_size(2) * (1 - alpha);
    
    floor_corners = [0 0 0;
                     intersection_size(1) 0 0;
                     intersection_size(1) intersection_size(2) 0;
                     0 intersection_size(2) 0];
    patch('Vertices', floor_corners, 'Faces', [1 2 3 4], 'FaceColor', ground_color, 'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 1);
    
    road_corners = [ix_start iy_start 0.01;
                    ix_end iy_start 0.01;
                    ix_end iy_end 0.01;
                    ix_start iy_end 0.01];
    patch('Vertices', road_corners, 'Faces', [1 2 3 4], 'FaceColor', road_color, 'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 1);
    
    BL = [0 0; ix_start 0; ix_start iy_start; 0 iy_start];
    plot_building(BL, intersection_size(3), building_color, wall_alpha, 'Bottom-Left');
    BR = [ix_end 0; intersection_size(1) 0; intersection_size(1) iy_start; ix_end iy_start];
    plot_building(BR, intersection_size(3), building_color, wall_alpha, 'Bottom-Right');
    TL = [0 iy_end; ix_start iy_end; ix_start intersection_size(2); 0 intersection_size(2)];
    plot_building(TL, intersection_size(3), building_color, wall_alpha, 'Top-Left');
    TR = [ix_end iy_end; intersection_size(1) iy_end; intersection_size(1) intersection_size(2); ix_end intersection_size(2)];
    plot_building(TR, intersection_size(3), building_color, wall_alpha, 'Top-Right');
    
    hTx = plot3(params.pTx(1), params.pTx(2), params.pTx(3), '^', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    hRx = plot3(params.pRx(1), params.pRx(2), params.pRx(3), 's', 'MarkerSize', 15, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    text(params.pTx(1), params.pTx(2), params.pTx(3)+1, sprintf('Tx (%.1f, %.1f, %.1f)', params.pTx), 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(params.pRx(1), params.pRx(2), params.pRx(3)+1, sprintf('Rx (%.1f, %.1f, %.1f)', params.pRx), 'FontSize', 10, 'HorizontalAlignment', 'center');
    
    plot3([params.pTx(1) params.pTx(1)], [params.pTx(2) params.pTx(2)], [0 params.pTx(3)], 'r:', 'LineWidth', 1);
    plot3([params.pRx(1) params.pRx(1)], [params.pRx(2) params.pRx(2)], [0 params.pRx(3)], 'b:', 'LineWidth', 1);
    plot3([params.pTx(1) params.pRx(1)], [params.pTx(2) params.pRx(2)], [params.pTx(3) params.pRx(3)], 'k--', 'LineWidth', 1);
    
    xlabel('X Position (m)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y Position (m)', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('Z Position (m)', 'FontSize', 12, 'FontWeight', 'bold');
    title('3D Urban Intersection Visualization', 'FontSize', 14, 'FontWeight', 'bold');
    legend([hTx, hRx], 'Location', 'best', 'FontSize', 10);
    
    view(45, 30);
    z_max = max([params.pTx(3), params.pRx(3)]) + 2;
    axis([0 intersection_size(1) 0 intersection_size(2) 0 z_max]);
    grid on;
    set(gca, 'FontSize', 10);
    hold off;
end