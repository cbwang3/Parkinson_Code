% Plot Raw and Filtered signals for PD and non-PD
function [figure_counter] = plotSignals(subject_name,subject,gaits_info,gait,figurenumber,dtitle,subtitle,matrix,pd,non_pd)
if subject_name(7) == 'A'
    dtitle = ['PD ',dtitle];
    figure(figurenumber); set(gcf, 'name', dtitle);
    ha(subject) = subplot(2, 5, pd);
    plot(matrix(:,1),matrix(:, gait),'r');  %matrix(:, 1),matrix(:, 2),'b',
    title(subtitle);
else
    dtitle = ['non-PD ',dtitle];
    figure(figurenumber+1); set(gcf, 'name', dtitle);
    ha(subject) = subplot(2, 3, non_pd);
    plot(matrix(:,1),matrix(:, gait),'r');  %matrix(:, 1),matrix(:, 2),'b'
    title(subtitle);
    
end
linkaxes(ha);
figure_counter = figurenumber + 2;