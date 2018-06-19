%% main_parkinsons
%file to load data, perform machine learning 
kav001A = struct('Tentries_acc', Tentries_acc, 'Tentries_gyro', Tentries_gyro,... 
'Tentries_orien', Tentries_orien);

save('kav001A', 'kav001A');

%% loading and plotting data
%each table has the time stamp as the last row, so we just plot
%the data from the first entry to end-1.

%% Load data from mat files 
%1) have to open mat files
%2) separate the data from the structure

data_acc = table2array(kav001A.Tentries_acc);
kav001A_acc = data_acc(1:(end-1),:); %%%
data_gyro = table2array(kav001A.Tentries_gyro);
kav001A_gyro = data_gyro(1:(end-1), :);
data_orien = table2array(kav001A.Tentries_orien);
kav001A_orien = data_orien(1:(end-1), :);


%% plot data
%first column is a time stamp
figure;
plot(data_acc_nts(:, 1), data_acc_nts(:, 2));
% hold on;
% plot(data_acc_nts(:, 1), data_acc_nts(:, 3)); plot(data_acc_nts(:, 1), data_acc_nts(:, 4));
title("Raw Accelerometer Data");
xlabel("m/s^2");
ylabel("Time Stamp (ms)");
legend('x', 'y', 'z');
hold off;

%% 
figure(2);
plot(data_gyro_nts(:, 1), data_gyro_nts(:, 2));hold on;
plot(data_gyro_nts(:, 1), data_gyro_nts(:, 3)); plot(data_gyro_nts(:, 1), data_gyro_nts(:, 4));
title("Raw Gyroscope Data");
xlabel("rad/sec");
ylabel("Time Stamp (ms)");
legend('x', 'y', 'z');
hold off;

figure(3);
plot(data_orien_nts(:, 1), data_orien_nts(:, 2));hold on;
plot(data_orien_nts(:, 1), data_orien_nts(:, 3)); plot(data_orien_nts(:, 1), data_orien_nts(:, 4));
title("Raw Orientation Data");
xlabel("");
ylabel("Time Stamp (ms)");
legend('x', 'y', 'z');
hold off;

%% plot spectrogram 
plotter = Plotter_p(); %create Plotter object, in folder 5-Plotting

% Explore Frequency Domain
%our data is in m/s^2 and deg/sec, how to make spectogram
plotter.plotSpectrogram(kav001A_acc(:,2),'Spectrogram','Sample','Frequency [Hz]');
plotter.plotSpectrogram(kav031B_acc(:,2),'Spectrogram','Sample','Frequency [Hz]');

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
wtrec(1:6, :) = wt(1:6, :);

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