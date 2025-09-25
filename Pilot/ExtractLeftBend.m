% Define the folder path
folderPath = 'N:\VLOpt\Pilot\S211381\LB\';

% Get a list of all text files in the folder
fileList = dir(fullfile(folderPath, '*.txt'));

% Loop through each file
for k = 1:length(fileList)
    % Get the full path of the current file
    FileName = fullfile(folderPath, fileList(k).name);
    TF = contains(fileList(k).name, "Shortened");
    if TF ==1
        continue
    end
    % Read the file
    TempFile = readtable(FileName);
    
    column_data = TempFile.My;
    sign_changes = sign(column_data);
    
    % Find the indexes where the sign changes (crossing zero threshold)
    cross_zero_indexes = find(diff(sign_changes) ~= 0);
    
    % Display the indexes
    disp(['Indexes where values cross the zero threshold for file: ', fileList(k).name]);
    disp(cross_zero_indexes);
    
    keep_index = true(size(cross_zero_indexes));
    
    % Iterate through the list to compare each index with the others
    for i = 1:length(cross_zero_indexes)
        if keep_index(i) % Only process if this index is not already marked for deletion
            % Compare current index with subsequent indexes
            for j = i+1:length(cross_zero_indexes)
                if abs(cross_zero_indexes(j) - cross_zero_indexes(i)) <= 100
                    % Mark the subsequent index for deletion
                    keep_index(j) = false;
                end
            end
        end
    end
    
    % Keep only the indexes that are not marked for deletion
    filtered_indexes = cross_zero_indexes(keep_index);
    
    % Display the filtered indexes
    disp(['Filtered indexes for file: ', fileList(k).name]);
    disp(filtered_indexes);
    
    % Select the rows based on filtered indexes (adjust this as needed)
    FinalTable = TempFile(filtered_indexes(5,1):filtered_indexes(6,1), :);
    
    % Define the output file name
    [~, name, ~] = fileparts(FileName);
    txtFileName = fullfile(folderPath, [name, '_LeftBendShortened.txt']);
    
    % Write the output to a new text file
    writetable(FinalTable, txtFileName, 'WriteVariableNames', true, 'Delimiter', '\t');
end
