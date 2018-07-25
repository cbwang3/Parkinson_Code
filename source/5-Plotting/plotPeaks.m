function [figurenumber] = plotPeaks(matrix,gait,peakLocs,positivePeaks,figure_counter,label,subtitle,x_axis,y_axis,neg_peakLocs,neg_peaks,pd_counter,nonpd_counter)
if subtitle(7) == 'A'
        label = ['PD ',label];
        figure(figure_counter); set(gcf, 'name', label); 
        subplot(2, 5, pd_counter);
else
        label = ['non-PD ',label];
        figure(figure_counter+1); set(gcf, 'name', label); 
        subplot(2, 3, nonpd_counter);
end
title(subtitle);
xlabel(x_axis); ylabel(y_axis);
hold on;
plot(matrix(:,1),matrix(:, gait),'r');
plot(peakLocs, positivePeaks, 'k.');

figurenumber = figure_counter + 2;