function plot_cdf(data, label)
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
    sorted_data = sort(data(:));
    N = length(sorted_data);
    
    if N == 0
        error('Data array is empty.');
    end
    
    % Ensure percentile value p is between 0 and 1
    p = max(min(p, 1), 0);
    
    rank = p * (N - 1) + 1;
    rank_floor = floor(rank);
    rank_ceil = ceil(rank);
    
    % Clamp indices to valid range
    rank_floor = max(rank_floor, 1);
    rank_ceil = min(rank_ceil, N);
    
    weight_ceil = rank - rank_floor;
    p_val = (1 - weight_ceil) * sorted_data(rank_floor) + weight_ceil * sorted_data(rank_ceil);
end