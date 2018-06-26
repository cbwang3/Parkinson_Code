%% main_parkinsons
%file to load data, perform machine learning 
%don't need to use this code once the files are saved
kav001A = struct('Tentries_acc', Tentries_acc, 'Tentries_gyro', Tentries_gyro,... 
'Tentries_orien', Tentries_orien);

save('kav001A', 'kav001A');

%% loading data
%each table has the time stamp as the last row, so we just plot
%the data from the first entry to end-1.
%for each file loaded; change name when extracting tables from saved mat
%structures
data_acc = table2array(kav031B.Tentries_acc);
data_acc = data_acc(1:(end-1),:); %%%
data_gyro = table2array(kav031B.Tentries_gyro);
data_gyro = data_gyro(1:(end-1), :);
data_orien = table2array(kav031B.Tentries_orien);
data_orien = data_orien(1:(end-1), :);

%reformat time stamp

%% interpolation
%uses interpolate_data function
%this is so all the timestamps for a patient are the same 
matrix2 = interpolate_data(data_acc, data_gyro, data_orien);

%% plot data using the plot_data function
%can use signal analyzer for now
%% 

plotter = Plotter_p();

plotter.plotSpectrogram(matrix1(:, 2), "Raw Acc Data", "Sample", "Frequency");

% [sp, fp, tp] = pspectrum(matrix(:, 2), matrix(:, 1), 'spectrogram');
% figure;
% mesh(tp, fp, sp);
% % view(-15, 60);
% xlabel('Time (ms)');
% ylabel('Frequency (Hz)');
% 
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