%% init_processing function to create data tables

function entries = parseData(data, num_var)
%2d array for holding at, x, y, z
%data: 
addpath(genpath('13_subjects_2'));
entries = zeros(length(data), num_var);

for a = 1:length(data)
    C = strsplit(data(a), ',');
    row = [];
    end_entry = data(a);
    if a == length(data) && end_entry(1:3)== "Snd" 
        break
        
        %last entry 
%         disp(C);
%         D = strsplit(C(4),' ');
%         D = strjoin(strsplit(D(6), ':'), '');
%         D = str2double(D);
%         entries(end,1) = D;
    else
        for i = 1: length(C)
        D = strsplit(C(i), ':');
        D = str2double(D);
        row = [row, D(end)];
        end
        entries(a, :) = row;
    end
    
   
end
%add column for m?
