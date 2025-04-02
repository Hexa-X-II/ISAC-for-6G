

%% 1) Housekeeping
clear; close all; clc;

%% 2) Intersection Parameters
FoV_size = [50, 50, 5];   % [width, depth, building height]
alpha    = 0.3;           % Building boundary fraction
grid_size = 100;           % Grid resolution

% Initialize simulation parameters (must have initialize_params.m)
params = initialize_params(FoV_size);

%% 3) Generate Grid and Building Coordinates
x_full = linspace(0, FoV_size(1), grid_size);
y_full = linspace(0, FoV_size(2), grid_size);
[X_full, Y_full] = meshgrid(x_full, y_full);

% Building corners
ix_start = FoV_size(1)*alpha;  
ix_end   = FoV_size(1)*(1-alpha);
iy_start = FoV_size(2)*alpha;
iy_end   = FoV_size(2)*(1-alpha);

% Identify building patches in corners
in_building = ( (X_full <= ix_start & Y_full <= iy_start) | ...  % bottom-left
                (X_full >= ix_end   & Y_full <= iy_start) | ...  % bottom-right
                (X_full <= ix_start & Y_full >= iy_end)   | ...  % top-left
                (X_full >= ix_end   & Y_full >= iy_end));        % top-right
diary('log.txt');  % Start logging output to log.txt

% Points outside the buildings
outside_building_mask = ~in_building;
X_valid = X_full(outside_building_mask);
Y_valid = Y_full(outside_building_mask);

% Sorted unique x,y positions (for FIM grid)
valid_pts = [X_valid'; Y_valid'];
x_positions = sort(unique(valid_pts(1, :)));
y_positions = sort(unique(valid_pts(2, :)));



%% 4) Compute LOS Mask
% intersection_create_mask must be defined
mask = intersection_create_mask(x_positions, y_positions, params.pTx(1:2), FoV_size, alpha);

%% 5) Plot the Intersection LOS Mask + Buildings
figure('Name','Intersection LOS Mask','Position',[100 100 600 500]);

% (a) Draw imagesc first
imagesc(x_positions, y_positions, mask);
colormap(gray);
colorbar;
axis equal; axis xy;
xlabel('X (m)'); ylabel('Y (m)');
title('Intersection LOS Mask (white = LOS available)');
hold on;

% (b) Draw building patches as solid red rectangles
% Bottom-left
patch([0, ix_start, ix_start, 0], [0, 0, iy_start, iy_start], ...
      'r', 'FaceColor','r', 'FaceAlpha',1, 'EdgeColor','none');
% Bottom-right
patch([ix_end, FoV_size(1), FoV_size(1), ix_end], [0, 0, iy_start, iy_start], ...
      'r', 'FaceColor','r', 'FaceAlpha',1, 'EdgeColor','none');
% Top-left
patch([0, ix_start, ix_start, 0], [iy_end, iy_end, FoV_size(2), FoV_size(2)], ...
      'r', 'FaceColor','r', 'FaceAlpha',1, 'EdgeColor','none');
% Top-right
%patch([ix_end, FoV_size(1), FoV_size(1), ix_end], [iy_end, iy_end, FoV_size(2), iy_end], ...
%      'r', 'FaceColor','r', 'FaceAlpha',1, 'EdgeColor','none');
patch([ix_end, FoV_size(1), FoV_size(1), ix_end], [iy_end, iy_end, FoV_size(2), FoV_size(2)], ...
      'r', 'FaceColor','r', 'FaceAlpha',1, 'EdgeColor','none');

% (c) Plot Tx and Rx
plot(params.pTx(1), params.pTx(2), 'b^', 'MarkerSize',10, 'LineWidth',2); % Tx: blue triangle
plot(params.pRx(1), params.pRx(2), 'gs', 'MarkerSize',10, 'LineWidth',2); % Rx: green square

% (d) Create custom legend handles (outside the axes)
% Dummy handle for Building: square marker
dummyTx       = plot(NaN, NaN, 'b^', 'MarkerSize',10, 'LineWidth',2);
dummyRx       = plot(NaN, NaN, 'gs', 'MarkerSize',10, 'LineWidth',2);
dummyBuilding = plot(NaN, NaN, 'rs', 'MarkerFaceColor','r','MarkerSize',10);
lgd = legend([dummyTx, dummyRx, dummyBuilding], ...
             {'Transmitter','Receiver','Building'}, ...
             'Location','eastoutside');  % <--- places legend to the right
set(lgd, 'Box','on');  % optional: box around legend
hold off;

%% 6) FIM Computation Over the Grid
% Identify valid walls for reflection
valid_walls = [];
for idx = 1:length(params.walls)
    if is_valid_reflection(params.pTx, params.pRx, params.walls(idx))
         valid_walls = [valid_walls, idx]; 
    end
end

figure('Name','3D Intersection Visualization','Position',[100 100 800 600]);
    % This function plots the 3D intersection along with reflection points and rays
valid_walls = plot_intersection_with_reflections(params, FoV_size, alpha);

fixed_value = sqrt((FoV_size(1)/2)^2 + (FoV_size(2)/2)^2 + (FoV_size(3)/2)^2);

% Preallocate
TPEB_grid   = zeros(length(x_positions), length(y_positions));
CRB_xT_grid = zeros(length(x_positions), length(y_positions));
CRB_yT_grid = zeros(length(x_positions), length(y_positions));
CRB_zT_grid = zeros(length(x_positions), length(y_positions));
flag_grid   = zeros(length(x_positions), length(y_positions));


for i = 1:length(x_positions)
    for j = 1:length(y_positions)
        if mask(j,i) == 1   % LOS available
            params.pT = [x_positions(i); y_positions(j); 1];  % user equipment pos (z=1)
            
            % Compute FIM and TPEB
            [FIM, ~, flag] = calculate_complete_fim_3d(params, valid_walls);
            [TPEB, CRB_xT, CRB_yT, CRB_zT, flag] = ...
                calculate_position_error_bound_3d(FIM, params, flag);
            
            TPEB_grid(i,j)   = TPEB;
            CRB_xT_grid(i,j) = CRB_xT;
            CRB_yT_grid(i,j) = CRB_yT;
            CRB_zT_grid(i,j) = CRB_zT;
            flag_grid(i,j)   = flag;
        else
            % blocked
            TPEB_grid(i,j)   = fixed_value;
            CRB_xT_grid(i,j) = fixed_value;
            CRB_yT_grid(i,j) = fixed_value;
            CRB_zT_grid(i,j) = fixed_value;
            flag_grid(i,j)   = 1;
        end
    end
end

%% Plot TPEB Heatmap
    %% --- Modify TPEB for visualization ---
    modified_TPEB = TPEB_grid;
    modified_TPEB(modified_TPEB > fixed_value) = fixed_value;
    modified_TPEB(flag_grid == 2) = fixed_value;

    figure;
    valid_TPEB_flag0 = modified_TPEB(flag_grid == 0);
    minTick_TPEB = floor(log10(min(valid_TPEB_flag0)));
    maxTick_TPEB = ceil(log10(max(valid_TPEB_flag0)));
    hHeat = imagesc(x_positions, y_positions, log10(modified_TPEB)');
    if minTick_TPEB >= maxTick_TPEB
        minTick_TPEB = minTick_TPEB - 1;
        maxTick_TPEB = maxTick_TPEB + 1;
    end
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
    
    % Overlay flag rectangles
    overlay_flag_rectangles(x_positions, y_positions, flag_grid');
    
    
    patch([0, ix_start, ix_start, 0], [0, 0, iy_start, iy_start], 'red', 'FaceAlpha', 1, 'EdgeColor', 'none');
    patch([ix_end, FoV_size(1), FoV_size(1), ix_end], [0, 0, iy_start, iy_start], 'red', 'FaceAlpha', 1, 'EdgeColor', 'none');
    patch([0, ix_start, ix_start, 0], [iy_end, iy_end, FoV_size(2), FoV_size(2)], 'red', 'FaceAlpha', 1, 'EdgeColor', 'none');
    %patch([ix_end, FoV_size(1), FoV_size(1), ix_end], [iy_end, iy_end, FoV_size(2), iy_end], 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([ix_end, FoV_size(1), FoV_size(1), ix_end], [iy_end, iy_end, FoV_size(2), FoV_size(2)], ...
      'red', 'FaceColor','r', 'FaceAlpha',1, 'EdgeColor','none');

    xlim([0 FoV_size(1)]);
    ylim([0 FoV_size(2)]);
    
    hTx = plot(params.pTx(1), params.pTx(2), 'r^', 'MarkerSize', 10);
    hRx = plot(params.pRx(1), params.pRx(2), 'rs', 'MarkerSize', 10);
    plot([0 FoV_size(1)], [0 0], 'k-', 'LineWidth', 2);
    plot([FoV_size(1) FoV_size(1)], [0 FoV_size(2)], 'k-', 'LineWidth', 2);
    plot([FoV_size(1) 0], [FoV_size(2) FoV_size(2)], 'k-', 'LineWidth', 2);
    plot([0 0], [FoV_size(2) 0], 'k-', 'LineWidth', 2);
    
    hOutsideFoV = patch('Visible','on', 'FaceColor','w', 'EdgeColor','none', 'DisplayName','Outside FoV');
    hUnresolved = patch('Visible','on', 'FaceColor','k', 'EdgeColor','none', 'DisplayName','Unresolved Position');
    hBuilding = patch('Visible','on', 'FaceColor','red', 'EdgeColor','none', 'DisplayName','Building');
    legend([hTx, hRx, hBuilding, hOutsideFoV, hUnresolved], ...
           {'Tx', 'Rx', 'Building', 'Outside FoV', 'Unresolved Position'}, 'Location', 'northeast');
    set(gca, 'FontSize', 14);
    set(gcf, 'Color', 'white');
    
    saveas(gcf, [num2str(params.fc/1e9) 'GHz_Heatmap.pdf']);
    saveas(gcf, [num2str(params.fc/1e9) 'GHz_Heatmap.fig']);

%% 8) Plot TPEB CDF
cdf_mask = (flag_grid == 0) | (flag_grid == 2);
cdf_data = modified_TPEB(cdf_mask);

figure('Name','TPEB CDF','Position',[300 100 600 400]);
plot_cdf(cdf_data, 'Target PEB (m)');  % must have plot_cdf.m
set(gca, 'FontSize',14);
set(gcf, 'Color','white');
saveas(gcf, [num2str(params.fc/1e9) 'GHz_CDF.pdf']);
saveas(gcf, [num2str(params.fc/1e9) 'GHz_CDF.fig']);

%% 9) Save Workspace
save([num2str(params.fc/1e9) 'GHz_workspace.mat']);

%% 10) Display Statistics
display_statistics('Target PEB', TPEB_grid(cdf_mask));  % must have display_statistics.m
display_statistics('X Position RMSE', sqrt(CRB_xT_grid(cdf_mask)));
display_statistics('Y Position RMSE', sqrt(CRB_yT_grid(cdf_mask)));
display_statistics('Z Position RMSE', sqrt(CRB_zT_grid(cdf_mask)));

diary off;  % Stop logging output
