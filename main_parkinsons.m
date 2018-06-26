%% main_parkinsons
%file to load data, perform machine learning 
%don't need to use this code once the files are saved
kav001A = struct('Tentries_acc', Tentries_acc, 'Tentries_gyro', Tentries_gyro,... 
'Tentries_orien', Tentries_orien);

save('kav001A', 'kav001A');
%% 
PD = ["001A", "002A","004A", "010A", "115A", "118A", "120A", "215A", "218A", "220A"]';
non_PD = ["031B", "079B", "111B", "211B", "121B", "211B"]';
for subject = 1:length(PD)
    load(strcat('kav',PD(subject),'_acc.mat'));
%     load(strcat('kav',extractBefore(patients(patient), 4),'_gyro.mat'));
%     load(strcat('kav',extractBefore(patients(patient), 4),'_orien.mat'));
    %load(strcat('kav',patients(patient),'.mat'))
    data_acc = table2array(Tentries_acc);
%     data_gyro = table2array(Tentries_gyro);

    %plot spectrograms
    figure(1);
    subplot(2, 5, subject);
    frequencyLimits = [0 100]/pi; %Normalized frequency (*pi rad/sample)
    leakage = 0.2;
    overlapPercent = 50;
    
    pspectrum(data_acc(:, 2), 'spectrogram','FrequencyLimits',frequencyLimits, ...
    'Leakage',leakage, 'OverlapPercent',overlapPercent);
    title('Spectrogram');

    % fourier transform 
    Y = fft(data_acc(:, 2));
    L = size(data_acc(:, 2), 1);
    %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
    P2 = abs(Y);
    P1 = P2(1:L/2+1);
    %Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.

    f = (1000/15)*[0:(L/2)]/L;
    figure(2);
    subplot(2, 5, subject);
    plot(f,P1);
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
end

for subject = 1:length(non_PD)
    load(strcat('kav',non_PD(subject),'_acc.mat'));
%     load(strcat('kav',extractBefore(patients(patient), 4),'_gyro.mat'));
%     load(strcat('kav',extractBefore(patients(patient), 4),'_orien.mat'));
    %load(strcat('kav',patients(patient),'.mat'))
    data_acc = table2array(Tentries_acc);
%     data_gyro = table2array(Tentries_gyro);

    %plot spectrograms
    figure(3);
    subplot(2, 3, subject);
    frequencyLimits = [0 100]/pi; %Normalized frequency (*pi rad/sample)
    leakage = 0.2;
    overlapPercent = 50;
    
    pspectrum(data_acc(:, 2), 'spectrogram','FrequencyLimits',frequencyLimits, ...
    'Leakage',leakage, 'OverlapPercent',overlapPercent);
    title('Spectrogram');

    % fourier transform 
    Y = fft(data_acc(:, 2));
    L = size(data_acc(:, 2), 1);
    %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
    P2 = abs(Y);
    P1 = P2(1:L/2+1);
    %Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.

    f = (1000/15)*[0:(L/2)]/L;
    figure(4);
    subplot(2, 3, subject);
    plot(f,P1);
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
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