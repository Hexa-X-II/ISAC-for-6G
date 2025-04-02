function display_statistics(label, data)
    fprintf('\n%s Statistics:\n', label);
    fprintf('Mean: %.4f\n', mean(data));
    fprintf('Median (50th percentile): %.4f\n', prctile(data, 50));
    fprintf('90th percentile: %.4f\n', prctile(data, 90));
    fprintf('95th percentile: %.4f\n', prctile(data, 95));
    fprintf('Standard Deviation: %.4f\n', std(data));
end

