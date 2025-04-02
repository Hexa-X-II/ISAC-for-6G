figure;
% Calculate min and max ticks for TPEB grid separately
% minTick_TPEB = floor(log10(min(TPEB_grid(:))));
% maxTick_TPEB = ceil(log10(max(TPEB_grid(TPEB_grid<max(TPEB_grid)))));
h1 = imagesc(x_points, y_points, dd');%, [minTick_TPEB, maxTick_TPEB]);
set(gca, 'YDir', 'normal');
% Mark outside-FoV (default) cells as black
set(h1, 'AlphaData', double(dd'>=4.7950));
% set(h1, 'AlphaData', double(dd'<4.7950 & dd'>4.66));
cb1 = colorbar;
title('Total distance');
xlabel('X Position (m)');
ylabel('Y Position (m)');
axis equal tight;
% cb1.Ticks = minTick_TPEB:1:maxTick_TPEB;
cb1.Label.String = 'd1+d2 (m)';
% cb1.TickLabels = arrayfun(@(x) sprintf('10^{%d}', x), cb1.Ticks, 'UniformOutput', false);

% Overlay transmitter, receiver, and room boundaries
hold on;
hTxCommon=plot(params.pTx(1), params.pTx(2), 'r^', 'MarkerSize', 10);
hRxCommon=plot(params.pRx(1), params.pRx(2), 'rs', 'MarkerSize', 10);
plot([0 5], [0 0], 'k-', 'LineWidth', 2);
plot([5 5], [0 5], 'k-', 'LineWidth', 2);
plot([5 0], [5 5], 'k-', 'LineWidth', 2);
plot([0 0], [5 0], 'k-', 'LineWidth', 2);
legend([hTxCommon, hRxCommon], {'Tx','Rx'}, 'Location', 'northeast');  % Place legend at the top right corner
set(gcf, 'Color', 'white');
set(gca,'fontsize', 14)