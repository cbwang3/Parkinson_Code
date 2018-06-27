%% main_parkinsons
all_subjects = ["001A", "002A","004A", "010A", "115A", "118A", "120A", "215A", "218A", "220A", "031B", "079B", "111B", "211B", "121B", "211B"]';

for subject = 1:2
    load(strcat('kav',all_subjects(subject),'_acc.mat'));
    load(strcat('kav',all_subjects(subject),'_gyro.mat'));
    load(strcat('kav',all_subjects(subject),'_orien.mat'));
    Tentries_acc = table2array(Tentries_acc); Tentries_gyro = table2array(Tentries_gyro); Tentries_orien = table2array(Tentries_orien);
    matrix = interpolate_data(Tentries_acc(1:end-1, :), Tentries_gyro(1:end-1, :), Tentries_orien(1:end-1, :));
    save(strcat('kav', all_subjects(subject), '_main.mat'), 'matrix');
end
%figure;
    
%% apply low pass filter to smooth data


% plot(data_acc(:, 1) - data_acc(1, 1), data_acc(:, 2));
% hold on;
cutoff = 1/100; %high cutoff -> more inclusive 
[b,a] = butter(1,cutoff); %filter params for cutoff of 0.2
data_acc_sm = zeros(size(data_acc));
data_acc_sm(:,2:4) = filter(b,a,data_acc(:,2:4)); %applying filter to all accel data


%subplot(2, 8, subject);
plot(data_acc(:, 1) - data_acc(1, 1), data_acc_sm(:, 2));
%     plotter.plotSignalBig(data_acc_sm(:,2),'Filtered','Sample','Acceleration [g]');

%end

%%
PD = all_subjects(1:10, :);
non_PD = all_subjects(11:end, :);

for subject = 1:length(PD)
    %accelerometer only
    load(strcat('kav',PD(subject),'_acc.mat'));
    data_acc = table2array(Tentries_orien);
%     load(strcat('kav',extractBefore(patients(patient), 4),'_gyro.mat'));
%     load(strcat('kav',extractBefore(patients(patient), 4),'_orien.mat')); 
%     data_gyro = table2array(Tentries_gyro);
    data_acc = data_acc(1:end-1, :); %remove time stamp
    
    %plot raw signal
    figure(1);
    subplot(2, 5, subject);
    plot(data_acc(:, 1) - data_acc(1,1), data_acc(:, 2)); hold on;
    %apply filter to take out noise
    cutoff = 1/100; %high cutoff -> more inclusive 
    [b,a] = butter(1,cutoff); %filter params for cutoff of 0.2
    data_acc_sm = zeros(size(data_acc));
    data_acc_sm(:,2:4) = filter(b,a,data_acc(:,2:4)); %applying filter to all accel data
    %plot filtered signal
    plot(data_acc(:, 1) - data_acc(1,1), data_acc_sm(:, 2));
%     %plot spectrograms
%     figure(2);
%     subplot(2, 5, subject);
%     frequencyLimits = [0 100]/pi; %Normalized frequency (*pi rad/sample)
%     leakage = 0.2;
%     overlapPercent = 50;
%     
%     pspectrum(data_acc(:, 2), 'spectrogram','FrequencyLimits',frequencyLimits, ...
%     'Leakage',leakage, 'OverlapPercent',overlapPercent);
%     title('Spectrogram');
% 
%     % fourier transform 
%     Y = fft(data_acc(:, 2));
%     L = size(data_acc(:, 2), 1);
%     %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
%     P2 = abs(Y);
%     P1 = P2(1:L/2+1);
%     %Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.
% 
%     f = (1000/15)*[0:(L/2)]/L;
%     figure(2);
%     subplot(2, 5, subject);
%     plot(f,P1);
%     title('Single-Sided Amplitude Spectrum of X(t)')
%     xlabel('f (Hz)')
%     ylabel('|P1(f)|')
end

for subject = 1:length(non_PD)
    load(strcat('kav',non_PD(subject),'_acc.mat'));
%     load(strcat('kav',extractBefore(patients(patient), 4),'_gyro.mat'));
%     load(strcat('kav',extractBefore(patients(patient), 4),'_orien.mat'));
    %load(strcat('kav',patients(patient),'.mat'))
    data_acc = table2array(Tentries_acc);
%     data_gyro = table2array(Tentries_gyro);
    data_acc = data_acc(1:end-1, :);
    %plot raw signal
    figure(2);
    subplot(2, 3, subject);
    plot(data_acc(:, 1) - data_acc(1,1), data_acc(:, 2)); hold on;
    %apply filter to take out noise
    cutoff = 1/100; %high cutoff -> more inclusive 
    [b,a] = butter(1,cutoff); %filter params for cutoff of 0.2
    data_acc_sm = zeros(size(data_acc));
    data_acc_sm(:,2:4) = filter(b,a,data_acc(:,2:4)); %applying filter to all accel data
    %plot filtered signal
    plot(data_acc(:, 1) - data_acc(1,1), data_acc_sm(:, 2));

%     %plot spectrograms
%     figure(3);
%     subplot(2, 3, subject);
%     frequencyLimits = [0 100]/pi; %Normalized frequency (*pi rad/sample)
%     leakage = 0.2;
%     overlapPercent = 50;
%     
%     pspectrum(data_acc(:, 2), 'spectrogram','FrequencyLimits',frequencyLimits, ...
%     'Leakage',leakage, 'OverlapPercent',overlapPercent);
%     title('Spectrogram');
% 
%     % fourier transform 
%     Y = fft(data_acc(:, 2));
%     L = size(data_acc(:, 2), 1);
%     %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
%     P2 = abs(Y);
%     P1 = P2(1:L/2+1);
%     %Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.
% 
%     f = (1000/15)*[0:(L/2)]/L;
%     figure(4);
%     subplot(2, 3, subject);
%     plot(f,P1);
%     title('Single-Sided Amplitude Spectrum of X(t)')
%     xlabel('f (Hz)')
%     ylabel('|P1(f)|')
end

%% interpolation (ignore this for now)
%uses interpolate_data function
%this is so all the timestamps for a patient are the same 
matrix2 = interpolate_data(data_acc, data_gyro, data_orien);

%% plot data using the plot_data function
%can use signal analyzer for now
%% 

plotter = Plotter_p();

plotter.plotSpectrogram(data_acc(:, 2), "Raw Acc Data 221B", "Sample", "Frequency");

% [sp, fp, tp] = pspectrum(matrix(:, 2), matrix(:, 1), 'spectrogram');
% figure;
% mesh(tp, fp, sp);
% % view(-15, 60);
% xlabel('Time (ms)');
% ylabel('Frequency (Hz)');

%% fourier transform 


Y = fft(data_acc(:, 2));
L = size(data_acc(:, 2), 1);
%Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2 = abs(Y);
P1 = P2(1:L/2+1);
%Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.

f = (1000/15)*(0:(L/2))/L;
figure;
plot(f,P1);
title('Single-Sided Amplitude Spectrum of X(t) 221B')
xlabel('f (Hz)')
ylabel('|P1(f)|')
hold on;


%% Filter 
%try three different filters:

%butterworth %design a lowpass filter, experiment with different types of 
%cutoff frequencies
Wn = 0.2; %cutoff frequency
[B, A] = butter(1, Wn);
data_acc_nts_b = data_acc_nts(:, 2:4);
data_acc_nts_b = filter(B,A,data_acc_nts_b); %applying filter to all accel data

figure(4);
plot(data_acc_nts(:, 1), data_acc_nts_b(:, 1));hold on;
plot(data_acc_nts(:, 1), data_acc_nts_b(:, 2)); plot(data_acc_nts(:, 1), data_acc_nts_b(:, 3));
title("Filtered Accelerometer Data");
xlabel("m/s^2");
ylabel("Time Stamp (ms)");
legend('x', 'y', 'z');
hold off;

%fourier

%% wavelet transformation
%
wt = modwt(kav001A_acc(:, 2));
%analyze individual wavelets from the decomposition
for a  = 1:length(wt(:, 1))
    subplot(6, 2, a); %m x n plot
    plot(wt(a, :));
end

%take out columns 7-11 for reconstruction
wtrec = zeros(size(wt));
wtrec(1:5, :) = wt(1:5, :);
modified_signal = imodwt(wtrec);
figure(5);
plot(kav001A_acc(:, 2)); hold on;
plot(modified_signal);



%%
plot(data_acc_nts(:, 1), data_acc_nts_b(:, 2)); plot(data_acc_nts(:, 1), data_acc_nts_b(:, 3));
title("Filtered Accelerometer Data w/ Wavelet");
xlabel("m/s^2");
ylabel("Time Stamp (ms)");
legend('x', 'y', 'z');
hold off;