%parse data
%use init_processing along with parseData function

%% main_parkinsons
%Nx7 matrices, saved by ID
all_subjects = ["001A", "002A","004A", "010A", "115A", "118A", "120A", "215A", "218A", "220A",  "031B", "079B", "111B", "211B", "121B", "221B"]';

%% loading data, only have to use once
for subject = 1:length(all_subjects)
    load(strcat('kav',all_subjects(subject),'_acc.mat'));
    load(strcat('kav',all_subjects(subject),'_gyro.mat'));
    load(strcat('kav',all_subjects(subject),'_orien.mat'));
    Tentries_acc = table2array(Tentries_acc); Tentries_gyro = table2array(Tentries_gyro); Tentries_orien = table2array(Tentries_orien);
    matrix = interpolate_data(Tentries_acc(1:end-1, :), Tentries_gyro(1:end-1, :), Tentries_orien(1:end-1, :));
    save(strcat('kav', all_subjects(subject), '_main.mat'), 'matrix');
end
%figure;

%%

clf;
for subject = 1:length(all_subjects)
    id = char(all_subjects(subject));
    load(strcat('kav',id,'_main.mat'));

% apply low pass filter to smooth data
    sfq = 100; %sampling frequency in Hz
    cfq =10; %cutoff frequency in Hz
    low_cutoff = cfq/(sfq/2); %high cutoff -> more inclusive 
    [b,a] = butter(1,low_cutoff, 'low'); %filter params for cutoff of 0.2
    data_acc_sm = zeros(size(matrix));
    data_acc_sm(:,2:end) = filter(b,a,matrix(: ,2:end)); %applying filter to all accel data
    %plot parkinsons vs non parkinsons 
%     if id(4) == 'A'
%         figure(1); set(gcf, 'name', 'PD Raw and Filtered Signals');
%         subplot(2, 5, subject);
%         plot(matrix(:,1),data_acc_sm(:, 2),'r');  %matrix(:, 1),matrix(:, 2),'b',      
%         title(strcat('kav',all_subjects(subject)));
%     else
%         figure(2); set(gcf, 'name', 'non-PD Raw and Filtered Signals');
%         subplot(2, 3, subject-10);
%         plot(matrix(:,1),data_acc_sm(:, 2),'r');  %matrix(:, 1),matrix(:, 2),'b'  
%         title(strcat('kav',all_subjects(subject)));
%     end
    matrix(:, 2:end) = data_acc_sm(:, 2:end);
%{
%% spectrograms
    %plot spectrograms
    frequencyLimits = [0 100]/pi; %Normalized frequency (*pi rad/sample)
    leakage = 0.2;
    overlapPercent = 50;
    
    if id(4) == 'A'
        figure(3); set(gcf, 'name', 'PD Spectrogram'); 
        subplot(2, 5, subject);
        pspectrum(data_acc_sm(:, 2), 'spectrogram','FrequencyLimits',frequencyLimits, ...
        'Leakage',leakage, 'OverlapPercent',overlapPercent);
        title('Spectrogram');
    else
        figure(4); set(gcf, 'name', 'non-PD Spectrogram'); 
        subplot(2, 3, subject-10);
        pspectrum(data_acc_sm(:, 2), 'spectrogram','FrequencyLimits',frequencyLimits, ...
        'Leakage',leakage, 'OverlapPercent',overlapPercent);
        title('Spectrogram');
    end
%% fourier transform
 % fourier transform 
    Y = fft(matrix(:, 2));
    L = size(matrix(:, 2), 1);
    %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
    P2 = abs(Y); P1 = P2(1:L/2+1);
    %Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.
    f = (1000/15)*[0:(L/2)]/L;
    
    if id(4) == 'A'
        figure(5); set(gcf, 'name', 'PD Spectrum'); subplot(2, 5, subject);
        plot(f,P1);
        title('PD Single-Sided Amplitude Spectrum of X(t)');
        xlabel('f (Hz)'); ylabel('|P1(f)|');
    else
        figure(6); set(gcf, 'name', 'non-PD Spectrum'); subplot(2, 5, subject);
        subplot(2, 3, subject-10);
        plot(f,P1);
        title('non-PD Single-Sided Amplitude Spectrum of X(t)');
        xlabel('f (Hz)'); ylabel('|P1(f)|');
    end
%}
%% wavelet transformation - for easier detection of peaks
    wt = modwt(matrix(:, 2));
%     figure;
% %     analyze individual wavelets from the decomposition
%     for a  = 1:length(wt(:, 1)) 
%         subplot(6, 2, a); %m x n plot
%         plot(wt(a, :));
%     end
    % take out columns 7-11 for reconstruction
    wtrec = zeros(size(wt));
    wtrec(5:10, :) = wt(5:10, :);
    modified_signal = imodwt(wtrec);    
% 
    if id(4) == 'A'
        figure(7); set(gcf, 'name', 'Reconstructed Wavelet');
        subplot(2, 5, subject);
        plot(matrix(:, 1), modified_signal); hold on;       
        title(strcat('kav',all_subjects(subject)));
    else
        figure(8); set(gcf, 'name', 'Reconstructed Wavelet');
        subplot(2, 3, subject-10);
        plot(matrix(:, 1), modified_signal); hold on;
        title(strcat('kav',all_subjects(subject)));
    end

%% peak detection - energy or x 
    maxpeak = max(modified_signal);
    for peakheight = maxpeak:-0.1:0
        x_plot = matrix(:, 2);
        [peaks, peakLocInds] = findpeaks(x_plot, 'minPeakHeight', peakheight, 'minPeakDistance', 30);
        if length(peaks)>=20
            length(peaks)
            break
        end
    end
    
    time_stamps = matrix(:, 1);
    figure;
    plot(matrix(:, 1), x_plot);
    hold on;
    peakLocs = time_stamps(peakLocInds);
    plot(peakLocs, peaks, 'r.');

    %detect positive peaks
    [peaks, peakLocInds] = findpeaks(modified_signal, 'minPeakHeight', 0.6, 'minPeakDistance', 30);
    time_stamps = matrix(:, 1);
    peakLocs = time_stamps(peakLocInds);
    
    %detect negative peaks
    [neg_peaks, neg_peakLocInds] = findpeaks(-modified_signal, 'minPeakHeight', 0.6, 'minPeakDistance', 30);
    neg_peakLocs = time_stamps(neg_peakLocInds);
    neg_peaks = -neg_peaks;
    
    if id(4) == 'A'
        figure(7); 
        hold on;
        subplot(2, 5, subject);
        plot(peakLocs, peaks, 'r.');  
        hold on; plot(neg_peakLocs, neg_peaks, 'k.');
    else
        figure(8); 
        hold on;
        subplot(2, 3, subject-10);
        plot(peakLocs, peaks, 'r.'); 
        hold on; plot(neg_peakLocs, neg_peaks, 'k.');
    end

%% segmentation
% loop to calculate the distances of the negative peaks from each positve
%peak
%timestamps: peakLocs, neg_peakLocs
    num_A = 0; num_B = 0;
    sum_A = 0; sum_B = 0;
    for ts = 1: length(peakLocs)
        %detect if there is a negative peak to the left
        closest_lneg_ind = max(find(neg_peakLocs < peakLocs(ts)));
        if ~isempty(closest_lneg_ind)
            sum_A = sum_A + (peakLocs(ts) - neg_peakLocs(closest_lneg_ind));
            num_A= num_A + 1;
        end
        %detect if there is negative peak to the right
        closest_rneg_ind = min(find(neg_peakLocs > peakLocs(ts)));
        if ~isempty(closest_rneg_ind)
            sum_B = sum_B + (neg_peakLocs(closest_rneg_ind) - peakLocs(ts));
            num_B = num_B + 1;
        end
    end

    %take the average of the distances; in ms 
    segment_A = sum_A/num_A;
    segment_B = sum_B/num_B;

    segmentStarts = peakLocs - segment_A;
    segmentEnds = peakLocs + segment_B;
    %segmentStartings you have the segment indices.

    %plot wavelet signals with the segments 
    if id(4) == 'A'
        figure(7); set(gcf, 'name', 'Reconstructed Wavelet');
        subplot(2, 5, subject);
        plot(matrix(:, 1), modified_signal); hold on;  
        for s = 1:length(segmentStarts)
            plot([segmentStarts(s), segmentStarts(s)], [-5 5], 'm');
        end
        for e = 1:length(segmentEnds)
            if segmentEnds(e) <= matrix(end, 1)
                plot([segmentEnds(e), segmentEnds(e)], [-5 5], 'k');
            end
        end
        title(strcat('kav',all_subjects(subject)));
    else
        figure(8); set(gcf, 'name', 'Reconstructed Wavelet');
        subplot(2, 3, subject-10);
        plot(matrix(:, 1), modified_signal); hold on;
        for s = 1:length(segmentStarts)
            plot([segmentStarts(s), segmentStarts(s)], [-5 5], 'm');
        end
        for e = 1:length(segmentEnds)
            if segmentEnds(e) <= matrix(end, 1)
                plot([segmentEnds(e), segmentEnds(e)], [-5 5], 'k');
            else
                segmentEnds = segmentEnds(1:e-1);
            end
        end
        title(strcat('kav',all_subjects(subject)));
    end

    
%     figure;
%     % plot(matrix(:, 1), modified_signal); hold on;
% 
%     plot(matrix(:, 1), matrix(:, 2)); hold on;
%     for s = 1:length(segmentStarts)
%         plot([segmentStarts(s), segmentStarts(s)], [-5 6], 'm');
%     end
% 
%     for e = 1:length(segmentEnds)
%         if segmentEnds(e) <= matrix(end, 1)
%             plot([segmentEnds(e), segmentEnds(e)], [-5 6], 'k');
%         end
%     end
end

%%
wt = modwt(matrix(:, 2));
figure;
%     analyze individual wavelets from the decomposition
for a  = 1:length(wt(:, 1)) 
    subplot(6, 2, a); %m x n plot
    plot(wt(a, :));
end
%take out columns 7-11 for reconstruction
wtrec = zeros(size(wt));
wtrec(5:10,:) = wt(5:10, :);
modified_signal = imodwt(wtrec);  

figure;
% plot(matrix(:, 1), matrix(:, 2));
% hold on;
plot(matrix(:, 1), modified_signal);




    %% peak detection - energy or x 
% %Energy calculation (e.g. for peak detection)
% energy_acc = matrix(:,2).^2 + matrix(:,3).^2 + matrix(:,4).^2;
% %plot energy
% 
%     if id(4) == 'A'
%         figure(7); set(gcf, 'name', 'PD Energy'); subplot(2, 5, subject);subplot(2, 5, subject);
%         plot(matrix(:, 1), energy_acc);
%         title('Energy, Acceleration g^2');
%     else
%         figure(8); set(gcf, 'name', 'non-PD Energy'); 
%         subplot(2, 3, subject-10);
%         plot(matrix(:, 1), energy_acc);
%         title('Energy, Acceleration g^2');
%     end


%data segmentation
%run for loop to determine segment A and segment B
%segmentation - try finding negative peaks for start and stop of each
%window 

% plot(matrix(:, 1), matrix(:, 3), matrix(:, 1), matrix(:,4));

%%
plot(data_acc_nts(:, 1), data_acc_nts_b(:, 2)); plot(data_acc_nts(:, 1), data_acc_nts_b(:, 3));
title("Filtered Accelerometer Data w/ Wavelet");
xlabel("m/s^2");
ylabel("Time Stamp (ms)");
legend('x', 'y', 'z');
hold off;



% if id(4) == 'A'
%     figure(5); subplot(2, 5, subject);
%     plot(matrix(:, 1),modified_signal);
% else
%     figure(6); subplot(2, 3, subject-10);
%     plot(matrix(:, 1),modified_signal);
% end


