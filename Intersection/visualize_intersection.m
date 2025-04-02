function main()
    % Clear workspace, close figures, and clear command window
    clear; close all; clc;
    
    % Define intersection dimensions
    intersection_size = [50, 50, 5];  % [width, depth, building height]
    
    % Initialize simulation parameters (which will use intersection_size to set boundaries)
    params = initialize_params(intersection_size);
    
    % Define grid parameters for 2D evaluation (x-y plane)
    grid_size = 5;  % grid_size x grid_size grid
    % Define alpha and compute grid over the intersection area (absolute coordinates)
    alpha = 0.3;
    % Define intersection boundaries explicitly
    ix_start = intersection_size(1) * alpha;  % e.g. 15 for 50 m width
    ix_end   = intersection_size(1) * (1 - alpha);  % e.g. 35 for 50 m width
    iy_start = intersection_size(2) * alpha;  % e.g. 15 for 50 m depth
    iy_end   = intersection_size(2) * (1 - alpha);  % e.g. 35 for 50 m depth

    % Define grid parameters for 2D evaluation (x-y plane)
    x_points = linspace(ix_start, ix_end, grid_size);
    y_points = linspace(iy_start, iy_end, grid_size);
    
    % For demonstration, we visualize the 3D intersection with buildings:
    visualizeintersection = true;
    if visualizeintersection
        figure('Name','3D Intersection Visualization','Position',[100 100 800 600]);
        plot_intersection(params, intersection_size, alpha);
    end
    
    %% --- Compute and Plot the LOS Mask ---
    % To use intersection_create_mask, we first convert our grid points from
    % absolute coordinates to local coordinates (centered at the intersection center).
    center = [intersection_size(1)/2, intersection_size(2)/2];
    % Local grid coordinates:
    x_local = x_points - center(1);
    y_local = y_points - center(2);
    
    % Convert Tx position to local coordinates (only x and y components)
    tx_local = params.pTx;
    tx_local(1:2) = tx_local(1:2) - center';
    
    % Define building corner positions relative to the intersection center.
    % Using the simulation boundaries:
    %   Bottom-left building touches the intersection at (ix_start,iy_start)
    %   Top-right building touches it at (ix_end,iy_end), where:
    ix_start = intersection_size(1)*alpha;  % 15 for intersection width 50
    ix_end   = intersection_size(1)*(1-alpha); % 35 for intersection width 50
    iy_start = intersection_size(2)*alpha;  % 15 for intersection depth 50
    iy_end   = intersection_size(2)*(1-alpha); % 35 for intersection depth 50
    % In local coordinates (subtracting the center [25,25]):
    %   bottom-left becomes [15-25,15-25] = [-10,-10]
    %   top-right becomes [35-25,35-25] = [10,10]
    % Similarly, we set the other two corners.
    tx_height = tx_local(3); % Use the same height as the Tx (and assumed target)
    corner_tr = [ ix_end - center(1); iy_end - center(2); tx_height ];  % [10; 10; tx_height]
    corner_tl = [ ix_start - center(1); iy_end - center(2); tx_height ];  % [15-25; 35-25] = [-10; 10]
    corner_bl = [ ix_start - center(1); iy_start - center(2); tx_height ];  % [-10; -10]
    corner_br = [ ix_end - center(1); iy_start - center(2); tx_height ];    % [10; -10]
    
    % Compute the LOS mask over the local grid
    mask = intersection_create_mask(x_local, y_local, tx_local, corner_tr, corner_tl, corner_bl, corner_br);
    
    % Plot the mask
    figure('Name','Intersection LOS Mask','Position',[100 100 600 500]);
    % Use imagesc to display the mask. Note: The mask matrix rows correspond to y.
    imagesc(x_local, y_local, mask);
    colormap(gray);
    colorbar;
    xlabel('x (m) (Local Coordinates)');
    ylabel('y (m) (Local Coordinates)');
    title('Intersection LOS Mask (true = LOS available)');
    axis equal;
    axis xy;

    % Loop over grid points (x-y positions)
    for i = 1:length(x_points)
        for j = 1:length(y_points)
            % Update target position (fix z at intersection mid-height)
            % clc;
            tic
            params.pT = [x_points(i); y_points(j); intersection_size(3)/2];
            
            % Calculate complete Fisher Information Matrix (FIM) for the current target position
            [FIM, ~, flag] = calculate_complete_fim_3d(params);
            
            % Calculate position error bounds 
            % (Assumed to return [TPEB, CRB_xT, CRB_yT, CRB_zT])
            [TPEB, CRB_xT, CRB_yT, CRB_zT, flag] = calculate_position_error_bound_3d(FIM, params, flag);
            
            % Store the computed bounds and flag
            TPEB_grid(i,j)   = TPEB;
            CRB_xT_grid(i,j) = CRB_xT;
            CRB_yT_grid(i,j) = CRB_yT;
            CRB_zT_grid(i,j) = CRB_zT;
            flag_grid(i,j)   = flag;
            toc
        end
    end

    % Modify TPEB for visualization: cap large values
    fixed_value = sqrt((intersection_size(1)/2)^2 + (intersection_size(2)/2)^2 + (intersection_size(3)/2)^2);
    modified_TPEB = TPEB_grid;
    modified_TPEB(modified_TPEB > fixed_value) = fixed_value;
    modified_TPEB(flag_grid == 2) = fixed_value;

    %% Plot TPEB Heatmap
    figure;
    valid_TPEB_flag0 = modified_TPEB(flag_grid == 0);
    minTick_TPEB = floor(log10(min(valid_TPEB_flag0)));
    maxTick_TPEB = ceil(log10(max(valid_TPEB_flag0)));
    hHeat = imagesc(x_points, y_points, log10(modified_TPEB)');
    %caxis([minTick_TPEB, maxTick_TPEB]);
    set(gca, 'YDir', 'normal');
    set(hHeat, 'AlphaData', double(TPEB_grid' ~= fixed_value));
    cb = colorbar;
    title('Target Position Error Bound (TPEB)');
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    axis equal tight;
    cb.Ticks = minTick_TPEB:1:maxTick_TPEB;
    cb.Label.String = 'Error (m)';
    cb.TickLabels = arrayfun(@(x) sprintf('10^{%d}', x), cb.Ticks, 'UniformOutput', false);
    hold on;
    overlay_flag_rectangles(x_points, y_points, flag_grid');
    
    % Overlay transmitter, receiver, and intersection boundaries
    hTx = plot(params.pTx(1), params.pTx(2), 'r^', 'MarkerSize', 10);
    hRx = plot(params.pRx(1), params.pRx(2), 'rs', 'MarkerSize', 10);
    plot([0 intersection_size(1)], [0 0], 'k-', 'LineWidth', 2);
    plot([intersection_size(1) intersection_size(1)], [0 intersection_size(2)], 'k-', 'LineWidth', 2);
    plot([intersection_size(1) 0], [intersection_size(2) intersection_size(2)], 'k-', 'LineWidth', 2);
    plot([0 0], [intersection_size(2) 0], 'k-', 'LineWidth', 2);
    
    % Create dummy patches for overlay flags in the legend
    hOutsideFoV = patch('Visible','on', 'FaceColor','w', 'EdgeColor','none', 'DisplayName','Outside FoV');
    hUnresolved = patch('Visible','on', 'FaceColor','k', 'EdgeColor','none', 'DisplayName','Unresolved Position');
    legend([hTx, hRx, hOutsideFoV, hUnresolved], ...
           {'Tx', 'Rx', 'Outside FoV', 'Unresolved Position'}, 'Location', 'northeast');
    set(gca, 'FontSize', 14);
    set(gcf, 'Color', 'white');
    saveas(gcf, [num2str(params.fc/1e9) 'GHz_Heatmap.pdf']);
    saveas(gcf, [num2str(params.fc/1e9) 'GHz_Heatmap.fig']);
    
    %% Plot TPEB CDF
    % Use only grid cells where flag is 0 or 2
    cdf_mask = (flag_grid == 0) | (flag_grid == 2);
    cdf_data = modified_TPEB(cdf_mask);
    
    figure;
    plot_cdf2(cdf_data, 'Target PEB (m)');
    set(gca, 'FontSize', 14);
    set(gcf, 'Color', 'white');
    saveas(gcf, [num2str(params.fc/1e9) 'GHz_CDF.pdf']);
    saveas(gcf, [num2str(params.fc/1e9) 'GHz_CDF.fig']);
    save([num2str(params.fc/1e9) 'GHz_workspace']);
    
    %% Display Statistics (only for cells with flag 0 or 2)
    display_statistics('Target PEB', TPEB_grid(cdf_mask));
    display_statistics('X Position RMSE', sqrt(CRB_xT_grid(cdf_mask)));
    display_statistics('Y Position RMSE', sqrt(CRB_yT_grid(cdf_mask)));
    display_statistics('Z Position RMSE', sqrt(CRB_zT_grid(cdf_mask)));

end

%% Helper Functions

function params = initialize_params(intersection_size)
    % Initialize the simulation parameters.
    params.c = 299792458;
    params.fc = 100e9;
    
    switch params.fc
        case 10e9
            params.delta_f = 120e3;
            params.Ntx_x = 2;
            params.Ntx_y = 2;
            params.Nrx_x = 2;
            params.Nrx_y = 2;
        case 60e9
            params.delta_f = 960e3;
            params.Ntx_x = 4;
            params.Ntx_y = 4;
            params.Nrx_x = 4;
            params.Nrx_y = 4;
        case 100e9
            params.delta_f = 1920e3;
            params.Ntx_x = 8;
            params.Ntx_y = 4;
            params.Nrx_x = 4;
            params.Nrx_y = 4;
    end
    params.K = 10;
    params.wavelength = params.c / params.fc;
    params.BW = params.delta_f * params.K;
    
    % Convert Es from dBm to Watts and normalize by bandwidth
    params.Es = 60;    
    params.Es = 10^(params.Es/10)*1e-3;
    params.Es = params.Es / params.BW;
    
    params.kB = 1.380649e-23;
    params.T0 = 290;
    params.N0 = params.kB * params.T0;
    params.NNF = 10^(7/10);
    
    params.Ntx = params.Ntx_x * params.Ntx_y;
    params.Nrx = params.Nrx_x * params.Nrx_y;
    
    % (Assuming generate_Q_matrix is defined elsewhere)
    params.QTx = generate_Q_matrix(params.Ntx_x, params.Ntx_y, params.wavelength);
    params.QRx = generate_Q_matrix(params.Nrx_x, params.Nrx_y, params.wavelength);
    params.N = params.Ntx * params.Nrx;
    
    % Positions for Tx, Rx, and Target (these remain unchanged)
    params.pTx = [25; 1; intersection_size(3)/2];
    params.pRx = [25; 49; intersection_size(3)/2];
    params.pT  = [2; 2; intersection_size(3)/2];
    
    % Euler angles (in radians)
    params.TxTheta = 0;
    params.TxPhi   = 0;
    params.TxPsi   = 0;
    params.RxTheta = 0;
    params.RxPhi   = pi;
    params.RxPsi   = 0;
    params.RTx = euler_rotation_matrix(params.TxTheta, params.TxPhi, params.TxPsi);
    params.RRx = euler_rotation_matrix(params.RxTheta, params.RxPhi, params.RxPsi);
    
    % Field of view parameters
    params.FoV_tx = pi - 0.0001;
    params.FoV_rx = pi - 0.0001;
    params.FoV_tx_elevation = pi - 0.0001;
    params.FoV_rx_elevation = pi - 0.0001;
    
    % --- Define Intersection Boundaries based on intersection size ---
    alpha = 0.3;
    ix_start = intersection_size(1)*alpha;         % e.g. 15 for 50 m width
    ix_end   = intersection_size(1)*(1 - alpha);     % e.g. 35 for 50 m width
    iy_start = intersection_size(2)*alpha;           % e.g. 15 for 50 m depth
    iy_end   = intersection_size(2)*(1 - alpha);       % e.g. 35 for 50 m depth
    
    % --- Define Walls (each with a point, normal, start and end) ---
    % We assume walls are defined at mid-height (z = intersection_size(3)/2).
    z_mid = intersection_size(3)/2;
    
    % 1. Ground (floor) -- for reflection (not used for building color)
    params.walls(1).p = [intersection_size(1)/2; intersection_size(2)/2; 0];
    params.walls(1).n = [0; 0; 1];  % Upward normal
    params.walls(1).start = [0; 0; 0];
    params.walls(1).end   = [intersection_size(1); 0; 0];
    
    % 2. Bottom-Left Building (occupies x: [0, ix_start], y: [0, iy_start])
    %    - East wall (faces into intersection)
    params.walls(2).p = [ix_start; iy_start/2; z_mid];
    params.walls(2).n = [1; 0; 0];
    params.walls(2).start = [ix_start; 0; z_mid];
    params.walls(2).end   = [ix_start; iy_start; z_mid];
    %    - North wall
    params.walls(3).p = [ix_start/2; iy_start; z_mid];
    params.walls(3).n = [0; 1; 0];
    params.walls(3).start = [0; iy_start; z_mid];
    params.walls(3).end   = [ix_start; iy_start; z_mid];
    
    % 3. Bottom-Right Building (occupies x: [ix_end, intersection_width], y: [0, iy_start])
    %    - West wall
    params.walls(4).p = [ix_end; iy_start/2; z_mid];
    params.walls(4).n = [-1; 0; 0];
    params.walls(4).start = [ix_end; 0; z_mid];
    params.walls(4).end   = [ix_end; iy_start; z_mid];
    %    - North wall
    params.walls(5).p = [(ix_end + intersection_size(1))/2; iy_start; z_mid];
    params.walls(5).n = [0; 1; 0];
    params.walls(5).start = [ix_end; iy_start; z_mid];
    params.walls(5).end   = [intersection_size(1); iy_start; z_mid];
    
    % 4. Top-Left Building (occupies x: [0, ix_start], y: [iy_end, intersection_depth])
    %    - East wall
    params.walls(6).p = [ix_start; (iy_end + intersection_size(2))/2; z_mid];
    params.walls(6).n = [1; 0; 0];
    params.walls(6).start = [ix_start; iy_end; z_mid];
    params.walls(6).end   = [ix_start; intersection_size(2); z_mid];
    %    - South wall
    params.walls(7).p = [ix_start/2; iy_end; z_mid];
    params.walls(7).n = [0; -1; 0];
    params.walls(7).start = [0; iy_end; z_mid];
    params.walls(7).end   = [ix_start; iy_end; z_mid];
    
    % 5. Top-Right Building (occupies x: [ix_end, intersection_width], y: [iy_end, intersection_depth])
    %    - West wall
    params.walls(8).p = [ix_end; (iy_end + intersection_size(2))/2; z_mid];
    params.walls(8).n = [-1; 0; 0];
    params.walls(8).start = [ix_end; iy_end; z_mid];
    params.walls(8).end   = [ix_end; intersection_size(2); z_mid];
    %    - South wall
    params.walls(9).p = [(ix_end+intersection_size(1))/2; iy_end; z_mid];
    params.walls(9).n = [0; -1; 0];
    params.walls(9).start = [ix_end; iy_end; z_mid];
    params.walls(9).end   = [intersection_size(1); iy_end; z_mid];
    
    params.sigma_RCS = 1.0;
    params.epsilon_r = 5.0;
    params.phase_offset = pi/3;
    params.target_phase_offset = pi/8;
end

function R = euler_rotation_matrix(theta, phi, psi)
    % Compute the Euler rotation matrix.
    Rx = [1, 0, 0;
          0, cos(theta), -sin(theta);
          0, sin(theta), cos(theta)];
    Rz = [cos(phi), -sin(phi), 0;
          sin(phi), cos(phi), 0;
          0, 0, 1];
    Ry = [cos(psi), 0, sin(psi);
          0, 1, 0;
          -sin(psi), 0, cos(psi)];
    R = Ry * Rx * Rz;
end

function plot_intersection(params, intersection_size, alpha)
    % Visualize the 3D intersection scenario with building volumes.
    wall_color = [1, 0, 0];   % red color for buildings/walls
    wall_alpha = 0.5;
    
    % Compute intersection boundaries
    ix_start = intersection_size(1) * alpha;
    ix_end   = intersection_size(1) * (1 - alpha);
    iy_start = intersection_size(2) * alpha;
    iy_end   = intersection_size(2) * (1 - alpha);
    
    hold on;
    
    % Plot the ground (floor)
    floor_corners = [0 0 0;
                     intersection_size(1) 0 0;
                     intersection_size(1) intersection_size(2) 0;
                     0 intersection_size(2) 0];
    patch('Vertices', floor_corners, 'Faces', [1 2 3 4], ...
          'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
      
    % Draw dashed lines to indicate road (intersection) boundaries
    plot([ix_start ix_start], [0 intersection_size(2)], 'k--', 'LineWidth', 2);
    plot([ix_end ix_end], [0 intersection_size(2)], 'k--', 'LineWidth', 2);
    plot([0 intersection_size(1)], [iy_start iy_start], 'k--', 'LineWidth', 2);
    plot([0 intersection_size(1)], [iy_end iy_end], 'k--', 'LineWidth', 2);

    % --- Plot the mask corners ---
    % Use the same boundaries as defined earlier in the function:
    ix_start = intersection_size(1) * alpha;
    ix_end   = intersection_size(1) * (1 - alpha);
    iy_start = intersection_size(2) * alpha;
    iy_end   = intersection_size(2) * (1 - alpha);
    % Use the same height as the Tx (which equals intersection_size(3)/2)
    z_corner = intersection_size(3)/2;

    % Define the corners in global coordinates
    corner_bl = [ix_start, iy_start, z_corner];
    corner_tl = [ix_start, iy_end,   z_corner];
    corner_tr = [ix_end,   iy_end,   z_corner];
    corner_br = [ix_end,   iy_start, z_corner];
    
    % Plot the corners with markers (here using a magenta pentagram)
    plot3(corner_bl(1), corner_bl(2), corner_bl(3), 'mp', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
    plot3(corner_tl(1), corner_tl(2), corner_tl(3), 'mp', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
    plot3(corner_tr(1), corner_tr(2), corner_tr(3), 'mp', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
    plot3(corner_br(1), corner_br(2), corner_br(3), 'mp', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
    
    % Optionally add text labels near each marker
    text(corner_bl(1), corner_bl(2), corner_bl(3), ' BL', 'Color','m','FontSize',10);
    text(corner_tl(1), corner_tl(2), corner_tl(3), ' TL', 'Color','m','FontSize',10);
    text(corner_tr(1), corner_tr(2), corner_tr(3), ' TR', 'Color','m','FontSize',10);
    text(corner_br(1), corner_br(2), corner_br(3), ' BR', 'Color','m','FontSize',10);
    
    % --- Plot building volumes (filled footprints) ---
    % Bottom-Left Building: occupies x in [0,ix_start], y in [0,iy_start]
    BL = [0 0; ix_start 0; ix_start iy_start; 0 iy_start];
    plot_building(BL, intersection_size(3), wall_color, wall_alpha);
    
    % Bottom-Right Building: occupies x in [ix_end, intersection_size(1)], y in [0,iy_start]
    BR = [ix_end 0; intersection_size(1) 0; intersection_size(1) iy_start; ix_end iy_start];
    plot_building(BR, intersection_size(3), wall_color, wall_alpha);
    
    % Top-Left Building: occupies x in [0,ix_start], y in [iy_end, intersection_size(2)]
    TL = [0 iy_end; ix_start iy_end; ix_start intersection_size(2); 0 intersection_size(2)];
    plot_building(TL, intersection_size(3), wall_color, wall_alpha);
    
    % Top-Right Building: occupies x in [ix_end, intersection_size(1)], y in [iy_end, intersection_size(2)]
    TR = [ix_end iy_end; intersection_size(1) iy_end; intersection_size(1) intersection_size(2); ix_end intersection_size(2)];
    plot_building(TR, intersection_size(3), wall_color, wall_alpha);
    
    % --- Plot walls as vertical planes using their start/end definitions ---
    for k = 2:length(params.walls)
        S = params.walls(k).start;  % defined at mid-height
        E = params.walls(k).end;
        vertices = [S(1), S(2), 0;
                    E(1), E(2), 0;
                    E(1), E(2), intersection_size(3);
                    S(1), S(2), intersection_size(3)];
        patch('Vertices', vertices, 'Faces', [1 2 3 4], ...
              'FaceColor', wall_color, 'FaceAlpha', wall_alpha, 'EdgeColor', 'k');
    end
    
    % Plot Tx, Rx, and Target positions
    hTx = plot3(params.pTx(1), params.pTx(2), params.pTx(3), 'r^', ...
                'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Tx');
    hRx = plot3(params.pRx(1), params.pRx(2), params.pRx(3), 'bs', ...
                'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Rx');
    hT  = plot3(params.pT(1), params.pT(2), params.pT(3), 'ko', ...
                'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Target');
    
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('3D Intersection Visualization');
    legend([hTx, hRx, hT], 'Tx', 'Rx', 'Target', 'Location', 'best');
    axis([0 intersection_size(1) 0 intersection_size(2) 0 intersection_size(3)]);
    grid on;
    view(3);
    hold off;
end

function plot_building(footprint, height, color, alpha_val)
    % Given a 2D footprint (Nx2 matrix) and a height, plot the building roof (filled area)
    numCorners = size(footprint,1);
    roofVertices = [footprint, repmat(height, numCorners, 1)];
    patch('Vertices', roofVertices, 'Faces', 1:numCorners, ...
          'FaceColor', color, 'FaceAlpha', alpha_val, 'EdgeColor', 'k');
end

function overlay_flag_rectangles(x_points, y_points, flag_grid)
    % Overlay white or black rectangles on the heatmap based on flag_grid.
    grid_size = size(flag_grid, 1);
    dx = x_points(2) - x_points(1);
    dy = y_points(2) - y_points(1);
    
    for i = 1:grid_size
        for j = 1:grid_size
            if flag_grid(i,j) == 1
                rectangle('Position', [x_points(j)-dx/2, y_points(i)-dy/2, dx, dy], ...
                          'FaceColor', 'w', 'EdgeColor', 'none');
            elseif flag_grid(i,j) == 2
                rectangle('Position', [x_points(j)-dx/2, y_points(i)-dy/2, dx, dy], ...
                          'FaceColor', 'k', 'EdgeColor', 'none');
            end
        end
    end
end

function plot_cdf2(data, label)
    % Plot the CDF of the input data and mark specific percentiles.
    sorted_data = sort(data);
    cdf = (1:length(sorted_data)) / length(sorted_data);
    
    plot(sorted_data, cdf, 'LineWidth', 2);
    grid on;
    title(['CDF of ' label]);
    xlabel(label);
    ylabel('Cumulative Probability');
    hold on;
    
    percentiles = [50, 90, 95];
    colors = {'r--', 'g--', 'b--'};
    for i = 1:length(percentiles)
        p_val = percentile(data, percentiles(i));
        p_frac = percentiles(i) / 100;
        plot([p_val p_val], [0 p_frac], colors{i}, 'LineWidth', 1.5);
        plot([min(sorted_data) p_val], [p_frac p_frac], colors{i}, 'LineWidth', 1.5);
    end
    legend({'CDF', '50th percentile','','90th percentile','','95th percentile',''}, 'Location', 'best');
    xlim([min(sorted_data), max(sorted_data)])
end

function p_val = percentile(data, p)
    % Compute the p-th percentile of the data.
    sorted_data = sort(data);
    n = length(data);
    rank = 1 + (p/100) * (n - 1);
    rank_floor = floor(rank);
    rank_ceil = ceil(rank);
    
    if rank_floor == rank_ceil
        p_val = sorted_data(rank_floor);
    else
        weight_ceil = rank - rank_floor;
        p_val = (1 - weight_ceil) * sorted_data(rank_floor) + weight_ceil * sorted_data(rank_ceil);
    end
end

function display_statistics(label, data)
    % Display mean, median, percentiles, and standard deviation for the given data.
    fprintf('\n%s Statistics:\n', label);
    fprintf('Mean: %.4f\n', mean(data));
    fprintf('Median (50th percentile): %.4f\n', percentile(data, 50));
    fprintf('90th percentile: %.4f\n', percentile(data, 90));
    fprintf('95th percentile: %.4f\n', percentile(data, 95));
    fprintf('Standard Deviation: %.4f\n', std(data));
end

