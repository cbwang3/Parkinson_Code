%% peak detection - energy or x 
%Energy calculation (e.g. for peak detection)
x_plot = matrix(:, 2);
[peaks, peakLocInds] = findpeaks(x_plot, 'minPeakHeight', 2, 'minPeakDistance', 30);
time_stamps = matrix(:, 1);
figure;
plot(matrix(:, 1), x_plot);
hold on;
peakLocs = time_stamps(peakLocInds);
plot(peakLocs, peaks, 'r.');
