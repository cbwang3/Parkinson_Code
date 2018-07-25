% Plots data for transformed features
function [figurenumber] = plotData(subject,x_value,y_value,figure_counter,label,subtitle,x_axis,y_axis,pd_counter,nonpd_counter)
chartitle = char(subtitle);
if chartitle(7) == 'A'
        label = ['PD ',label];
        figure(figure_counter); set(gcf, 'name', label); 
        subplot(2, 5, pd_counter);
else
        label = ['non-PD ',label];
        figure(figure_counter+1); set(gcf, 'name', label);
        subplot(2, 3, nonpd_counter);
        
end
plot(x_value,y_value);
title(subtitle);
xlabel(x_axis); ylabel(y_axis);
figurenumber = figure_counter + 2;
