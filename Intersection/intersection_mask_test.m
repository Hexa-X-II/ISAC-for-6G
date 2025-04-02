
clear all;
close all;

% Grid
x_points = linspace(-50, 50, 500);
y_points = linspace(-50, 50, 500);

% Tx position
tx_pos = [2.0; 30.0; 1.0];

% Corners
corner_tr = [20.0; 20.0; 1.0]; % top right
corner_tl = [-25.0; 25.0; 1.0]; % top left
corner_bl = [-20.0; -20.0; 1.0]; % bottom left
corner_br = [20.0; -20.0; 1.0]; % bottom right

mask = intersection_create_mask( ...
    x_points, ...
    y_points, ...
    tx_pos, ...
    corner_tr, ...
    corner_tl, ...
    corner_bl, ...
    corner_br);

figure();
s = pcolor(x_points, y_points, cast(mask, "uint8")); hold on
scatter(tx_pos(1), tx_pos(2), "o", "blue");
scatter(corner_tr(1), corner_tr(2), "x", "red");
scatter(corner_tl(1), corner_tl(2), "x", "yellow");
scatter(corner_bl(1), corner_bl(2), "x", "green");
scatter(corner_br(1), corner_br(2), "x", "magenta");
s.EdgeColor = "none";
colormap("gray");
axis equal
legend("Mask", "Tx", "Top Right", "Top Left", "Bottom Left", "Bottom Right", "Location", "bestoutside");

% Example of how to use the mask in our case
img = NaN(length(y_points), length(x_points));
for j = 1:length(x_points)
    for i = 1:length(y_points)
        if mask(i, j)
            img(i, j) = x_points(j).^2 + y_points(i).^2; % arbitrary function (should be CRB)
        end
    end
end

figure();
s2 = pcolor(x_points, y_points, img);
s2.EdgeColor = "none";
axis equal