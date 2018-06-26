%function for interpolating our data
function matrix  = interpolate_data(acc, gyro, orien)
%takes in the accelerometer, gyroscope, and orientation data for a 
%subject and synchronizes the data. Matrix holds readjusted timestamp 
%(measurement for each ms) with the corresponding accel and gyro data.

%this is so all the timestamps for a patient are the same 
% find maximum from time stamps

%starting time stamp
max_start_ts = max([acc(1, 1) gyro(1, 1) orien(1, 1)]);

inds_st_acc = find(acc(:, 1) > max_start_ts); inds_st_gyro = find(gyro(:, 1) > max_start_ts);
inds_st_orien = find(orien(:, 1) > max_start_ts);

%ending time stamp
min_stop_ts = min([acc(end, 1) gyro(end, 1) orien(end, 1)]);
inds_end_acc = find(acc(:, 1) < min_stop_ts); inds_end_gyro = find(gyro(:, 1) < min_stop_ts);
inds_end_orien = find(orien(:, 1) < min_stop_ts);

%cut the signal from each index
new_data_acc = acc(inds_st_acc(1)-1:inds_end_acc(end)+1, :);
new_data_gyro = gyro(inds_st_gyro(1)-1:inds_end_gyro(end)+1, :);
new_data_orien = orien(inds_st_orien(1)-1:inds_end_orien(end)+1, :);

%interpolate on all axes for acc and gyro
tq_acc = new_data_acc(1, 1):1: new_data_acc(end,1); %new timestamp
vq_acc = interp1(new_data_acc(:, 1), new_data_acc(:, 2:4), tq_acc); %interpolate for all directions
start = find(tq_acc == max_start_ts);
finish = find(tq_acc == min_stop_ts);
% tq_acc = tq_acc(start:end)';
vq_acc = vq_acc(start:finish, :);

tq_gyro = new_data_gyro(1, 1):1: new_data_gyro(end,1); %new timestamp
vq_gyro = interp1(new_data_gyro(:, 1), new_data_gyro(:, 2:4), tq_gyro);
start = find(tq_gyro == max_start_ts);
finish = find(tq_gyro == min_stop_ts);
vq_gyro = vq_gyro(start:finish, :);

% create matrix
matrix = zeros(length(vq_acc), 7);
matrix(:, 1) = [0:length(vq_acc)-1]'; %time stamp
matrix(:, 2:4) = vq_acc;
matrix(:, 5:7) = vq_gyro;