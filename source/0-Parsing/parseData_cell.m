%% init_processing function to create data tables
function entries = parseData_cell(data, num_var)
%data is class cell
%2d array for holding at, x, y, z
%data: 
entries = zeros(length(data), num_var);

for a = 1:length(data)
    C = strsplit(data{a}, ',');
    row = [];
    end_entry = data{a};
    if a == length(data) && end_entry(1:3) == "Snd" %last entry 
        break
    else
        for i = 1: length(C)
            D = strsplit(C{i}, ':');
            D = D{end};
            row = [row, D];
        end
    end
    entries(a, :) = row;
   
   
end
%add column for m?
