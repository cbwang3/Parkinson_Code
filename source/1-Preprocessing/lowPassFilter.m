% Low-pass filter to smooth data

function [smoothed_data]  = lowPassFilter(matrix)
    % Sampling frequency in Hz
    sfq = 100;
    % cutoff frequency in Hz
    cfq =10; 
    
    low_cutoff = cfq/(sfq/2); % high cutoff -> more inclusive 
    [b,a] = butter(1,low_cutoff, 'low'); % filter params for cutoff of 0.2
    smoothed_data = zeros(size(matrix));
    smoothed_data(:,2:end) = filter(b,a,matrix(: ,2:end)); %applying filter to all accel data
    smoothed_data(:,1) = matrix(:,1);
    