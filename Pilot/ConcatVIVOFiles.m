inputDir = 'N:\MartineSingleSegment_S24\C210577_L2L3\10000steps_Gait'; 

% Get a list of all files in the input directory
Files = dir(fullfile(inputDir, '*.txt'));

FinalTable = table();
% Check if any files are found
if isempty(Files)
    disp('No TXT files found in the input directory.');
else

    % Loop through each file
    for k = 1:length(Files)
        A = contains(fileList(k).name, "Shortened");
        B = contains(fileList(k).name, "Resected");
    
        if A == 1 && B ==1
            
        end
        % Get the full path of the current file
        FileName = fullfile(inputDir, Files(k).name);
    
        TempFile = readtable(FileName);
    
        FinalTable = [FinalTable; TempFile];
    end    

end

txtFileName = "ConcatFile.txt";

writetable(FinalTable, txtFileName, 'WriteVariableNames', true, 'Delimiter', '\t');
