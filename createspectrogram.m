%% main_parkinsons
%file to load data, perform machine learning 

columnpatients = ["001A",'002A','004A', '010A', '031B', '079B', '111B', '115A', '118A', '120A', '121B', '211B', '215A', '218A', '220A']
patients = columnpatients'

for patient = 1:length(patients)
load(strcat('kav',extractBefore(patients(patient), 4),'_acc.mat'))
load(strcat('kav',extractBefore(patients(patient), 4),'_gyro.mat'))
load(strcat('kav',extractBefore(patients(patient), 4),'_orien.mat'))
%load(strcat('kav',patients(patient),'.mat'))
data_acc = table2array(Tentries_acc)
data_gyro = table2array(Tentries_gyro)
plotter = Plotter(); %create Plotter object, in folder 5-Plotting
plotter.plotSpectrogram(data_acc(:,2), strcat("Raw Acc Data ",patients(patient), "Sample", "Frequency"))

plotter = Plotter(); %create Plotter object, in folder 5-Plotting
plotter.plotSpectrogram(data_gyro(:,2), strcat("Raw Gyro Data ",patients(patient), "Sample", "Frequency"))

%struct('Tentries_acc', Tentries_acc, 'Tentries_gyro', Tentries_gyro,... 
%'Tentries_orien', Tentries_orien);

%save(strcat('kav',patients(patient)),strcat('kav',patients(patient)))

end

%% loading and plotting data
%each table has the time stamp as the last row, so we just plot
%the data from the first entry to end-1.

%% Load data from mat files 
%1) have to open mat files
%2) separate the data from the structure

nonPDpatients = dir('*B.mat'); 
for patientnumber = 1:length(nonPDpatients) 
    load(nonPDpatients(patientnumber))
    data_acc = table2array(kav111B.Tentries_acc)
    kav031B_acc = data_acc(1:(end-1),:)
end 

%data_gyro = table2array(kav031B.Tentries_gyro);
%kav001A_gyro = data_gyro(1:(end-1), :);

%data_orien = table2array(kav031B.Tentries_orien);
%kav001A_orien = data_orien(1:(end-1), :);

%% plot spectrogram 
plotter = Plotter_p(); %create Plotter object, in folder 5-Plotting
plotter.plotSpectrogram(data_acc(:,2), "Raw Acc Data 004A", "Sample", "Frequency");

% Explore Frequency Domain
%our data is in m/s^2 and deg/sec, how to make spectogram
%plotter.plotSpectrogram(kav001A_acc(:,2),'Spectrogram','Sample','Frequency [Hz]');
%plotter.plotSpectrogram(kav031B_acc(:,2),'Spectrogram','Sample','Frequency [Hz]');
