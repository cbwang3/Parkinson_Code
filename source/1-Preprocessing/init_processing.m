%% parsing Parkinson's daa
addpath(genpath('13_subjects_2'));

%shared variables
num_var = 4;
for i = 24:27
    chari = int2str(i);
    filenames_acc = ['kaveli4_acc_201805', chari, '*'];
    filenames_gyro = ['kaveli4_gyro_201805', chari, '*'];
    filenames_orien = ['kaveli4_orien_201805', chari, '*'];
    directories_acc = dir(filenames_acc);
    directories_gyro = dir(filenames_gyro);
    directories_orien = dir(filenames_orien);
    all_acc = [];
    all_gyro = [];
    all_orien = [];
    acc_length = length(directories_acc);
    gyro_length = length(directories_gyro);
    orien_length = length(directories_orien);
    
    for filename = 1:acc_length
        %% acceleration
        symbol = "%s";
        indexeddirectory = directories_acc(filename);
        fid = fopen(indexeddirectory.name);
        M = textscan(fid, symbol, ...
            'Delimiter', '{}', ...
            'CollectOutput', true );
        M = M{1,1};
        M = M(3:2:end-2);
        
        %Mtable = cell2mat(M);
        % newdata = table2array(M);
        % time_stamp_acc = M(end-2);
        
        %2) parse out the data for each sample. Sampling starts at 0
        %stats = Mtable(2:2:end);
        entries_acc = parseData_cell(M, num_var); %array
        
        Tentries_acc = array2table(entries_acc, 'VariableNames', {'at', 'x', 'y', 'z'});
        % Create the accelerometer array, generates a table for one file,
        % Tentries_acc
        % Add table of acc_array = [acc-array Tentries_acc]
        % this is included in for loop
        all_acc = [all_acc; Tentries_acc];
        disp(filename)
        disp(filename / acc_length)
        disp("acc")
    end
    savedfile = ['kaveli4_acc_201805', chari];
    save(savedfile, 'all_acc');
    %end
    
    
    %% orientation data
    for filename = 1:orien_length
    %imported data with {} delimiters
    %1) import data
    indexeddirectory = directories_orien(filename);
        fid = fopen(indexeddirectory.name);
        M = textscan(fid, "%s", ...
            'Delimiter', '{}', ...
            'CollectOutput', true );
        M = M{1,1};
        M = M(3:2:end-2);
        
        entries_orien = parseData_cell(M, num_var); %array
        
        Tentries_orien = array2table(entries_orien, 'VariableNames', {'at', 'x', 'y', 'z'});
       
        all_orien = [all_orien; Tentries_orien];
        disp(filename)
        disp(filename / orien_length)
        disp("orien")
    end
    savedfile = ['kaveli4_orien_201805', chari];
    save(savedfile, 'all_orien');
    
    %% gyroscope
    
    for filename = 1:gyro_length
    %imported data with {} delimiters
    %1) import data
    indexeddirectory = directories_gyro(filename);
        fid = fopen(indexeddirectory.name);
        M = textscan(fid, "%s", ...
            'Delimiter', '{}', ...
            'CollectOutput', true );
        M = M{1,1};
        M = M(3:2:end-2);
        
        entries_gyro = parseData_cell(M, num_var); %array
        
        Tentries_gyro = array2table(entries_gyro, 'VariableNames', {'at', 'x', 'y', 'z'});
       
        all_gyro = [all_gyro; Tentries_gyro];
        disp(filename)
        disp(filename / gyro_length)
        disp("gyro")
    end
    savedfile = ['kaveli4_gyro_201805', chari];
    save(savedfile, 'all_gyro');
    
end
    
