% Plot spectrograms
function [figure_counter] = plotSpectrograms(subject_name,gait,figurenumber,dtitle,subtitle,matrix,pd,non_pd)
frequencyLimits = [0 100]/pi; % Normalized frequency (*pi rad/sample)
    leakage = 0.2;
    overlapPercent = 50;
    
    if subject_name(7) == 'A'
        dtitle = ['PD',dtitle];
        figure(figurenumber); set(gcf, 'name', dtitle); 
        subplot(2, 5, pd);
        pspectrum(matrix(:, gait), 'spectrogram','FrequencyLimits',frequencyLimits, ...
        'Leakage',leakage, 'OverlapPercent',overlapPercent);
        title(subtitle);
    else
        dtitle = ['non-PD ',dtitle];
        figure(figurenumber+1); set(gcf, 'name', dtitle); 
        subplot(2, 3, non_pd);
        pspectrum(matrix(:, gait), 'spectrogram','FrequencyLimits',frequencyLimits, ...
        'Leakage',leakage, 'OverlapPercent',overlapPercent);
        title(subtitle);
       
    end
    figure_counter = figurenumber + 2;
end
