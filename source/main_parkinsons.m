%
% ________________________________________________________  ____  ______
% \______   \______   \_   _____/\______   \_____  \   _  \/_   |/  __  \
%  |     ___/|       _/|    __)_  |     ___//  ____/  /_\  \|   |>      <
%  |    |    |    |   \|        \ |    |   /       \  \_/   \   /   --   \
%  |____|    |____|_  /_______  / |____|   \_______ \_____  /___\______  /
%                   \/        \/                   \/     \/           \/
% Code for detecting Parkinson's Disease in patient subjects using gait data
% Supervised by Juan Haladjian and Professor Bernd Bruegge
% Cassia Wang
% Justin Zhu
% @ Technical University of Munich
% with Collaboration from University of Tampere


%% Manual Variables to be optimized by Unsupervised Learning
numFeatures = 15;
minPeakDistance = 30;
percent_test = .20; %20% test, 80% training


%% Add Paths for Helper Functions: Extraction, Selection, & Classification
addpath(genpath('0-Parsing'));
addpath(genpath('1-Preprocessing'));
addpath(genpath('2-FeatureExtraction'));
addpath(genpath('3-FeatureSelection'));
addpath(genpath('4-Classification'));
addpath(genpath('5-Plotting'));
addpath(genpath('mat_files'));
addpath(genpath('metadata'));
currentdirectory = pwd;
filedirectory = '/mat_files/';
figuredirectory = '/figures/';

%% Read Data Inputs
patient_number = 13;

patient_labels = readtable("patient_labels_kaveli.xlsx");
patient_id = convertCharsToStrings(table2array(patient_labels(:,1)));
patient_id = patient_id(1:patient_number);
patient_group = convertCharsToStrings(table2array(patient_labels(:,2)));
patient_group = patient_group(1:patient_number);

all_subjects = patient_id + patient_group;

all_subjects(5) = "kav111B";
all_subjects = [all_subjects; "kav211B"];
all_subjects(6) = "kav115A";
all_subjects = [all_subjects; "kav215A"];
all_subjects(7) = "kav118A";
all_subjects = [all_subjects; "kav218A"];
all_subjects(8) = "kav120A";
all_subjects = [all_subjects; "kav220A"];
all_subjects(9) = "kav121B";
all_subjects = [all_subjects; "kav221B"];

all_subjects = all_subjects(all_subjects~="kav027B"&all_subjects~="kav028B");
patient_number = size(all_subjects);



% Load Accelerometer, Gyroscope, and Orientation Data into cleanly parsed matrix
for subject = 1:length(all_subjects)
    load(strcat(all_subjects(subject),'_acc.mat'));
    load(strcat(all_subjects(subject),'_gyro.mat'));
    load(strcat(all_subjects(subject),'_orien.mat'));
    
    Tentries_acc = table2array(Tentries_acc);
    Tentries_gyro = table2array(Tentries_gyro);
    Tentries_orien = table2array(Tentries_orien);
    
    matrix = interpolate_data(Tentries_acc(1:end-1, :), Tentries_gyro(1:end-1, :), Tentries_orien(1:end-1, :));
    save(strcat(currentdirectory,filedirectory,all_subjects(subject), '_main.mat'), 'matrix');
end

%% Initializing variables for feature extraction
featureExtractor = FeatureExtractor();
nFeatures = featureExtractor.nFeatures;

PD_class = zeros(1,nFeatures);
nonPD_class = zeros(1, nFeatures);

% allTable contains indicies of each subjects' strides
num_strides = 1;
strides_per_subj = [1];
% Counts the number of figures
figure_counter = 1;

% matrix is Nx10 displayed like this for each kav patient subject:

gait_key = ["timestamp";"(x-acc)";"(y-acc)";"(z-acc)";"(x-gyro)";"y_gyro";"(z-gyro)";"(x-orien)";"(y-orien)";"(z-orien)"];
timestamp = 1;
x_acc = 2;
y_acc = 3;
z_acc = 4;
x_gyro = 5;
y_gyro = 6;
z_gyro = 7;
x_orien = 8;
y_orien = 9;
z_orien = 10;

% Specify what gait data to analyze
gaits = [x_acc];
gaits_size = numel(gaits);

% Create a counter for the number of patients with PD and non-PD
pd_counter = 0;
nonpd_counter = 0;

% Iterate over all subjects and gait data, transforming and plotting features
for subject = 1:length(all_subjects)
    figure_counter = 1;
    % load subject files
    subject_name = char(all_subjects(subject));
    subject_filename = strcat(subject_name,'_main.mat');
    load(subject_filename);
    
    if subject_name(7) == 'A'
        pd_counter = pd_counter + 1;
    else
        nonpd_counter = nonpd_counter + 1;
    end
    
    
    % Iterate for specified gait data
    for gait = 1:gaits_size
        
        gait_info = char(gait_key(gaits(gait)));
        gait_to_plot = gaits(gait);
        subtitle = [subject_name,' ',gait_info];
        x_value = matrix(:, gaits(gait));
        id = subject_name(7);
        % Smooth using Butterworth low-pass filter
        matrix = lowPassFilter(matrix)
        time_stamps = matrix(:, 1);
        
        % Plot PD Raw and Filtered Signals
        figure_counter = plotSignals(subject_name,subject,gait_info,gaits(gait),figure_counter,'Raw and Filtered Signals',subtitle,matrix,pd_counter,nonpd_counter);
        
        % Plot spectrograms
        figure_counter = plotSpectrograms(subject_name,gaits(gait),figure_counter,'Spectrogram',subtitle,matrix,pd_counter,nonpd_counter);
        
        %% Fourier transform
        Y = fft(matrix(:, gaits(gait)));
        L = size(matrix(:, gaits(gait), 1));
        % Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
        P2 = abs(Y); P1 = P2(1:L/2+1);
        % Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.
        f = (1000/15)*(0:(L/2))/L(1);
        
        % Plot Single-Sided Amplitude Spectrum of X(t)
        [figure_counter] = plotData(subject,f,P1,figure_counter,'Single-Sided Amplitude Spectrum',subtitle,'f (Hz)','|P1(f)|', pd_counter, nonpd_counter);
        
        %% Wavelet transformation - for easier detection of peaks
        
        wt = modwt(matrix(:, gait_to_plot));
        wtrec = zeros(size(wt));
        wtrec(5:10, :) = wt(5:10, :);
        modified_signal = imodwt(wtrec);
        
        [figure_counter] = plotData(subject,time_stamps,modified_signal,figure_counter,'Reconstructed Wavelet',subtitle,'time (s)','|P1(f)|',pd_counter,nonpd_counter);
        
        %% Peak detection
        % Use peak detection on the modified wavelet signal
        maxpeak = max(modified_signal);
        
        % Detect positive peaks
        % Optimization - detect twenty peaks (or more) with known 30 ms
        % minimimum peak distance
        
        for peakheight = maxpeak:-0.1:0
            [peaks, peakLocInds] = findpeaks(modified_signal, 'minPeakHeight', peakheight, 'minPeakDistance', minPeakDistance);
            if length(peaks)>=20
                length(peaks);
                break
            end
        end
        
        peakLocs = time_stamps(peakLocInds); %peakLocs are actual time stamps where peaks occurred
        signalpeaks = x_value(peakLocInds); % signalpeaks are actual x_values where peaks occured
        % Detect negative peaks
        [neg_peaks, neg_peakLocInds] = findpeaks(-modified_signal, 'minPeakHeight', 0.6, 'minPeakDistance', 30);
        neg_peakLocs = time_stamps(neg_peakLocInds);
        neg_peaks = -neg_peaks;
        
        % Mark peaks
        
        figure_counter = plotPeaks(matrix,...
            gait_to_plot,...
            peakLocs,...
            signalpeaks,...
            figure_counter,...
            'Marked Peaks',...
            subtitle,...
            '',...
            '',...
            neg_peakLocs,...
            neg_peaks,...
            pd_counter,...
            nonpd_counter);
    end
    
    %% Segmentation
    % loop to calculate the distances of the negative peaks from each positve
    % peak
    num_A = 0; num_B = 0;
    sum_A = 0; sum_B = 0; %running sum of differences of time stamp indices
    for ts = 1: length(peakLocInds)
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
    
    %take the average of the distances; by index of timestamp
    segment_A = floor(sum_A/(10*num_A)) * 10;
    segment_B = floor(sum_B/(10*num_B)) * 10;
    
    segmentStarts = peakLocs - segment_A; %indices
    %     length_of_segments(1 ,subject) = length(segmentStarts);
    segmentEnds = peakLocs + segment_B;
    %     length_of_segments(2, subject) = length(segmentEnds);
    %segmentStartings you have the segment indices.
    
    if subject == 2
        figure_counter = plotSegment(figure_counter, gait, matrix, segmentStarts, segmentEnds);
        
    end
    
    %plot original filtered signals with the segments
    %     if id(4) == 'A'
    %         figure(7); hold on;
    %         subplot(2, 5, subject);
    %
    %         %plot(matrix(:, 1), modified_signal); hold on;
    %         for s = 1:length(segmentStarts)
    %             plot([segmentStarts(s), segmentStarts(s)], [-10 10], 'm');
    %         end
    %         for e = 1:length(segmentEnds)
    %             if segmentEnds(e) <= matrix(end, 1)
    %                 plot([segmentEnds(e), segmentEnds(e)], [-10 10], 'k');
    %             else
    %                 segmentEnds = segmentEnds(1:e-1);
    %             end
    %         end
    %         title(strcat('kav',all_subjects(subject)));
    %     else
    %         figure(8); hold on;
    %         subplot(2, 3, subject-10);
    %         %plot(matrix(:, 1), modified_signal); hold on;
    %         for s = 1:length(segmentStarts)
    %             plot([segmentStarts(s), segmentStarts(s)], [-10 10], 'm');
    %         end
    %         for e = 1:length(segmentEnds)
    %             if segmentEnds(e) <= matrix(end,1)
    %                 plot([segmentEnds(e), segmentEnds(e)], [-10 10], 'k');
    %             else
    %                 segmentEnds = segmentEnds(1:e-1);
    %             end
    %         end
    %         title(strcat('kav',all_subjects(subject)));
    %     end
    
    %% Feature Extraction
    %min, max of segment
    %keep an array of data sample we are currently workig with
    fc1 = []; fc2 = [];%keeps count of #strides/patient
    
    for i = 1 : length(segmentEnds)
        start_ts = segmentStarts(i);
        end_ts = segmentEnds(i);
        
        if end_ts > matrix(end, 1)
            continue
        end
        start_idx = find(matrix(:, 1) == start_ts);
        end_idx = find(matrix(:, 1) == end_ts);
        segment = matrix(start_idx:end_idx, 1:end);
        
        [featureVector] = featureExtractor.extractFeaturesForSegment(segment); %%can't get class to work
        if id == 'A'
            fc1(i,:) = featureVector;
        else
            fc2(i, :) = featureVector;
        end
        
        num_strides = num_strides + 1;
    end
    
    strides_per_subj = [strides_per_subj num_strides];
    PD_class = [PD_class; fc1];
    nonPD_class = [nonPD_class; fc2];
    %Note: you could save the features using the save command for faster
    %loading next time
    
    
end

%% Merge feature tables
numFeatures = 15;
PD_class(:,nFeatures+1) = 1;%label every feature vector with a 1 to indicate a parkinson's subject
nonPD_class(:, nFeatures+1) = 0; %label every feature vector with a 0 to indicate a non-pd subject


%%
allTable = array2table([PD_class; nonPD_class]);
allTable.Properties.VariableNames = [featureExtractor.featureNames, 'label'];
allTable = allTable(~any(ismissing(allTable),2),:);

%% Normalize features
dataNormalizer = DataNormalizer();
dataNormalizer.fit(allTable);
allTable = dataNormalizer.normalize(allTable);

%% Feature Selection
featureSelector = FeatureSelector();
bestFeatures = featureSelector.findBestFeatures(allTable,numFeatures);
allTable = featureSelector.selectFeatures(allTable,bestFeatures);

%% take out samples for testing data
allArray = table2array(allTable);
%% 1st option: randomized. To use this, must take code out from the for loop to train and test the classifier

test_PD_samples = randperm(size(PD_class, 1),floor(percent_test*size(PD_class, 1)));

test_PD_class = allArray([test_PD_samples], :);
PD_class_indices = [1: size(PD_class, 1)];
train_PD_class= allArray(setdiff(PD_class_indices, test_PD_samples), :);

test_nonPD_samples = randperm(size(nonPD_class, 1),floor(percent_test*size(nonPD_class, 1)));

test_nonPD_class = allArray(size(PD_class, 1) + [test_nonPD_samples], :);
nonPD_class_indices = [1: size(nonPD_class, 1)];
train_nonPD_class= allArray(size(PD_class, 1) + setdiff(nonPD_class_indices, test_nonPD_samples), :);

trainTable = array2table([train_PD_class; train_nonPD_class]);
testTable = array2table([test_PD_class; test_nonPD_class]);

%% 2nd option: leave one subject out validation (going through all subjects).
for xtimes = 1:length(all_subjects)
    %% leave one out validation
    %choose a patient's test to leave out from strides_per_subj
    patient_i = xtimes;
    %leave strides from stride_per_subj(patient_i) to stride_per_subj(patient_i
    %+ 1)
    start_test = strides_per_subj(patient_i);
    end_test = strides_per_subj(patient_i + 1);
    test_samples = allArray(start_test+1:end_test, :);
    
    if start_test ~= 1 && end_test ~= length(allArray)
        train_samples = allArray([1:start_test-1, end_test+1:end], :);
    elseif start_test == 1
        train_samples = allArray(end_test+1:end, :);
    elseif end_test == length(allArray)
        train_samples = allArray(1:start_test-1, :);
    end
    
    trainTable = array2table(train_samples); testTable= array2table(test_samples);
    %% Train Classifier
    %this classifiers uses a predefined algorithm (SVM) with a polynomial
    %kernel
    
    %do own SVM training
    % svm_class_model = fitcsvm(trainTable(:, 1:end-1), trainTable(:, end));
    
    trainer = Trainer();
    % [beta_vals, bias]
    trainer.train(trainTable);
    %svm is sign(beta_vals*matrix + bias);
    
    %% Test Classifier
    
    % [label,score] = predict(svm_class_model,testTable(:, 1:end-1));
    testTable.Properties.VariableNames = trainTable.Properties.VariableNames;
    labels = trainer.test(testTable);
    
    shouldBeLabels = table2array(testTable(:,end));
    
    %compare labels to shouldBeLabels to get you accuracy
    
    %to test further algorithms, open the Classification Learner Tool in the
    %Matlab-Toolbox and select the variable 'table'
    
    
    %% Plot Results
    plotter = Plotter();
    confusionMatrix = confusionmat(shouldBeLabels,labels);
    plotter.plotConfusionMatrix(confusionMatrix,["non-Parkinsons","Parkinsons"]);
    figure_counter = figure_counter + 1;
end

%% Save and close all figures generated
figlist=findobj('type','figure');
[num_figs,~] = size(figlist);
for figure=1:num_figs
    titlenum = num_figs-figure+1;
    savedindex = num2str(titlenum);
    savedfigure = figlist(figure);
    saveas(savedfigure,[currentdirectory,figuredirectory,'fig_',savedindex]);
end
close all;
