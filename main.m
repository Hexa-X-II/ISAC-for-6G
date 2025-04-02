function main()
    % Clear workspace, close figures, and clear command window
    clear; close all; clc;
    
    % Initialize simulation parameters
    params = initialize_params();
    
    % Define grid parameters for 2D evaluation (x-y plane)
    grid_size = 25;         % Grid will be grid_size x grid_size
    room_size = [50, 50, 5];  % Room dimensions: [width, depth, height] in meters
    margin = 0.05;           % 10% margin from walls
    
    % Calculate x and y grid points (fix z at mid-height of the room)
    x_margin = room_size(1) * margin;
    y_margin = room_size(2) * margin;
    x_points = linspace(x_margin, room_size(1) - x_margin, grid_size);
    y_points = linspace(y_margin, room_size(2) - y_margin, grid_size);
    
    % Preallocate result matrices for the 2D grid
    TPEB_grid    = zeros(grid_size, grid_size);
    CRB_xT_grid  = zeros(grid_size, grid_size);
    CRB_yT_grid  = zeros(grid_size, grid_size);
    CRB_zT_grid  = zeros(grid_size, grid_size);
    flag_grid    = zeros(grid_size, grid_size);  % flag: 0 = normal, 1 = white overlay, 2 = black overlay
    
    % Toggle room visualization (set true to view the 3D room)
    visualizeRoom = false;
    if visualizeRoom
        figure('Name','3D Room Visualization','Position',[100 100 800 600]);
        plot_room(params);
    end
    
    % Loop over grid points (x-y positions)
    for i = 1:grid_size
        for j = 1:grid_size
            % Update target position (fix z at room mid-height)
            % clc;
            tic
            params.pT = [x_points(i); y_points(j); room_size(3)/2];
            
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
    fixed_value = sqrt((room_size(1)/2)^2 + (room_size(2)/2)^2 + (room_size(3)/2)^2);
    modified_TPEB = TPEB_grid;
    modified_TPEB(modified_TPEB > fixed_value) = fixed_value;
    modified_TPEB(flag_grid == 2) = fixed_value;

    %% Plot TPEB Heatmap
    figure;
    valid_TPEB_flag0 = modified_TPEB(flag_grid == 0);
    minTick_TPEB = floor(log10(min(valid_TPEB_flag0)));
    maxTick_TPEB = ceil(log10(max(valid_TPEB_flag0)));
    hHeat = imagesc(x_points, y_points, log10(modified_TPEB)', [minTick_TPEB, maxTick_TPEB]);
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
    
    % Overlay transmitter, receiver, and room boundaries
    hTx = plot(params.pTx(1), params.pTx(2), 'r^', 'MarkerSize', 10);
    hRx = plot(params.pRx(1), params.pRx(2), 'rs', 'MarkerSize', 10);
    plot([0 room_size(1)], [0 0], 'k-', 'LineWidth', 2);
    plot([room_size(1) room_size(1)], [0 room_size(2)], 'k-', 'LineWidth', 2);
    plot([room_size(1) 0], [room_size(2) room_size(2)], 'k-', 'LineWidth', 2);
    plot([0 0], [room_size(2) 0], 'k-', 'LineWidth', 2);
    
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

%% Helper Function Definitions

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

function params = initialize_params()
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
    params.K = 792;
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
    
    params.QTx = generate_Q_matrix(params.Ntx_x, params.Ntx_y, params.wavelength);
    params.QRx = generate_Q_matrix(params.Nrx_x, params.Nrx_y, params.wavelength);
    params.N = params.Ntx * params.Nrx;
    
    % Positions for Tx, Rx, and Target
    params.pTx = [25; 1; 2.5];
    params.pRx = [25; 49; 2.5];
    params.pT = [2; 2; 2.5];
    
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
    
    % Wall definitions
    params.walls(1).p = [25; 25; 0];
    params.walls(1).n = [0; 0; 1];
    params.walls(2).p = [50; 25; 2.5];
    params.walls(2).n = [-1; 0; 0];
    params.walls(3).p = [25; 25; 5];
    params.walls(3).n = [0; 0; -1];
    params.walls(4).p = [0; 25; 2.5];
    params.walls(4).n = [1; 0; 0];
    params.walls(5).p = [25; 50; 2.5];
    params.walls(5).n = [0; -1; 0];
    params.walls(6).p = [25; 0; 2.5];
    params.walls(6).n = [0; 1; 0];
    
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

function plot_room(params)
    % Visualize a 3D room with floor, ceiling, and walls.
    room_size = [5, 5, 5];  % Example room dimensions for visualization
    
    floor_corners = [0 0 0; room_size(1) 0 0; room_size(1) room_size(2) 0; 0 room_size(2) 0];
    ceiling_corners = [0 0 room_size(3); room_size(1) 0 room_size(3); room_size(1) room_size(2) room_size(3); 0 room_size(2) room_size(3)];
    
    hold on;
    patch('Vertices', floor_corners, 'Faces', [1 2 3 4], 'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch('Vertices', ceiling_corners, 'Faces', [1 2 3 4], 'FaceColor', 'magenta', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    patch('Vertices', [0 0 0; 0 room_size(2) 0; 0 room_size(2) room_size(3); 0 0 room_size(3)],...
          'Faces', [1 2 3 4], 'FaceColor', 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch('Vertices', [room_size(1) 0 0; room_size(1) room_size(2) 0; room_size(1) room_size(2) room_size(3); room_size(1) 0 room_size(3)],...
          'Faces', [1 2 3 4], 'FaceColor', 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch('Vertices', [0 0 0; room_size(1) 0 0; room_size(1) 0 room_size(3); 0 0 room_size(3)],...
          'Faces', [1 2 3 4], 'FaceColor', 'green', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch('Vertices', [0 room_size(2) 0; room_size(1) room_size(2) 0; room_size(1) room_size(2) room_size(3); 0 room_size(2) room_size(3)],...
          'Faces', [1 2 3 4], 'FaceColor', 'green', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Plot Tx, Rx, and Target positions
    hTx = plot3(params.pTx(1), params.pTx(2), params.pTx(3), 'r^', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Tx');
    hRx = plot3(params.pRx(1), params.pRx(2), params.pRx(3), 'bs', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Rx');
    hT  = plot3(params.pT(1), params.pT(2), params.pT(3), 'ko', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Target');
    
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    axis equal; grid on; view(3);
    title('3D Room Visualization');
    legend([hTx, hRx, hT], 'Location', 'best');
    hold off;
end