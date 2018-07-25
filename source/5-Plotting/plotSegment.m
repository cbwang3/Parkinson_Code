% Plot segments of Main Parkinson's
function [figurenumber] = plotSegment(figure_counter, gait, matrix, segmentStarts, segmentEnds)
figure(figure_counter);
plot(matrix(:, 1), matrix(:, gait)); hold on;
for s = 1:length(segmentStarts)
    plot([segmentStarts(s), segmentStarts(s)], [-10 10], 'm');
    
end
for e = 1:length(segmentEnds)
    if segmentEnds(e) <= matrix(end, 1)
        plot([segmentEnds(e), segmentEnds(e)], [-10 10], 'm');
    else
        segmentEnds = segmentEnds(1:e-1);
    end
end
figurenumber = figure_counter + 1;