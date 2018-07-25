%parse data
%use init_processing along with parseData function

%% main_parkinsons
%Nx7 matrices, saved by ID
all_subjects = ["001A", "002A", "004A", "010A", "115A", "118A", "120A", "215A", "218A", "220A",  "031B", "079B", "111B", "211B", "121B", "221B"]';

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

%initialize 
%length_of_segments = zeros(2, length(all_subjects));
%feature classes 
featureExtractor = FeatureExtractor();
nFeatures = featureExtractor.nFeatures;
PD_class = zeros(1,nFeatures);
nonPD_class = zeros(1, nFeatures);
num_strides = 1;
strides_per_subj = [1]; %keeps indicies of each subjects' strides in allTable
for subject = 1:length(all_subjects)
    id = char(all_subjects(subject));
    load(strcat('kav',id,'_main.mat'));
    
%     if id == "031B"
%         
%     
% apply low pass filter to smooth data
    sfq = 100; %sampling frequency in Hz
    cfq =10; %cutoff frequency in Hz
    low_cutoff = cfq/(sfq/2); %high cutoff -> more inclusive 
    [b,a] = butter(1,low_cutoff, 'low'); %filter params for cutoff of 0.2
    data_acc_sm = zeros(size(matrix));
    data_acc_sm(:,2:end) = filter(b,a,matrix(: ,2:end)); %applying filter to all accel data
    %plot parkinsons vs non parkinsons 
   
    if id(4) == 'A'
        figure(1); set(gcf, 'name', 'PD Raw and Filtered Signals');
        ha(subject) = subplot(2, 5, subject);
        plot(matrix(:,1),data_acc_sm(:, 3),'r');  %matrix(:, 1),matrix(:, 2),'b',      
       
        title(strcat('kav',all_subjects(subject)));
    else
        figure(2); set(gcf, 'name', 'non-PD Raw and Filtered Signals');
        ha(subject) = subplot(2, 3, subject-10);
        plot(matrix(:,1),data_acc_sm(:, 3),'r');  %matrix(:, 1),matrix(:, 2),'b'  
        title(strcat('kav',all_subjects(subject)));
        
%         plot(matrix(:,1), data_acc_sm(:, 4), 'm');
    end
    linkaxes(ha);
    matrix(:, 2:end) = data_acc_sm(:, 2:end);

%% spectrograms
%{
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
%{
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
%     if id(4) == 'A'
%         figure(7); set(gcf, 'name', 'Reconstructed Wavelet');
%         subplot(2, 5, subject);
%         plot(matrix(:, 1), modified_signal); hold on;       
%         title(strcat('kav',all_subjects(subject)));
%     else
%         figure(8); set(gcf, 'name', 'Reconstructed Wavelet');
%         subplot(2, 3, subject-10);
%         plot(matrix(:, 1), modified_signal); hold on;
%         title(strcat('kav',all_subjects(subject)));
%     end

%% peak detection - energy or x 

    %use peak detection on the modified wavelet signal 
    maxpeak = max(modified_signal);
    %detect positive peaks
    %optimization- detect twenty peaks (or more)
    for peakheight = maxpeak:-0.1:0
        %x_plot = matrix(:, 2);
        [peaks, peakLocInds] = findpeaks(modified_signal, 'minPeakHeight', peakheight, 'minPeakDistance', 30);
        if length(peaks)>=20
            length(peaks);
            break
        end
    end
    time_stamps = matrix(:, 1);
    peakLocs = time_stamps(peakLocInds); %peakLocs are actual time stamps where peaks occurred
    
    %detect negative peaks
    [neg_peaks, neg_peakLocInds] = findpeaks(-modified_signal, 'minPeakHeight', 0.6, 'minPeakDistance', 30);
    neg_peakLocs = time_stamps(neg_peakLocInds);
    neg_peaks = -neg_peaks;
    
%     %mark peaks on the wavelet signals 
%     if id(4) == 'A'
%         figure(7); 
%         hold on;
%         subplot(2, 5, subject);
%         plot(peakLocs, peaks, 'r.');  
%         hold on; plot(neg_peakLocs, neg_peaks, 'k.');
%     else
%         figure(8); 
%         hold on;
%         subplot(2, 3, subject-10);
%         plot(peakLocs, peaks, 'r.'); 
%         hold on; plot(neg_peakLocs, neg_peaks, 'k.');
%     end
%% segmentation
% loop to calculate the distances of the negative peaks from each positve
%peak
%timestamps: peakLocs, neg_peakLocs
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
       figure;
       plot(matrix(:, 1), matrix(:, 2)); hold on;  
       for s = 1:length(segmentStarts)
            plot([segmentStarts(s), segmentStarts(s)], [-10 10], 'm');
            
       end
       for e = 1:length(segmentEnds)
            if segmentEnds(e) <= matrix(end, 1)
                plot([segmentEnds(e), segmentEnds(e)], [-10 10], 'm');
            else
                segmentEnds = segmentEnds(1:e-1);
            end
       end
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
        if id(4) == 'A'
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
    %}
end

%% Merge feature tables
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
nFeatures = 10;
featureSelector = FeatureSelector(); 
bestFeatures = featureSelector.findBestFeatures(allTable,nFeatures);
allTable = featureSelector.selectFeatures(allTable,bestFeatures);

%% take out samples for testing data
allArray = table2array(allTable);
%% 1st option: randomized. To use this, must take code out from the for loop to train and test the classifier

percent_test = .20; %20% test, 80% training
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
    [beta_vals, bias] = trainer.trainSVM(trainTable);
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
end



%% IGNORE THE FOLLOWING SECTIONS
%{
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

%% Feature Extraction

featureExtractor = FeatureExtractor();
nFeatures = featureExtractor.nFeatures;
PD_class = zeros(length(segmentStarts),nFeatures);

%min, max of segment
for i = 1 : length(segmentStarts)
    start_ts = segmentStarts(i);
    end_ts = segmentEnds(i);
    start_idx = find(matrix(:, 1) == start_ts);
    end_idx = find(matrix(:, 1) == end_ts);
    segment = matrix(start_idx:end_idx, 1:end);
    [featureVector] = featureExtractor.extractFeaturesForSegment(segment); %%can't get class to work
    PD_class(i,:) = featureVector;
end


%% segmentation
% loop to calculate the distances of the negative peaks from each positve
%peak
%timestamps: peakLocs, neg_peakLocs
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
   
%}
