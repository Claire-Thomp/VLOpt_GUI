FileName = 'N:\VLOpt\Pilot\DispControl_Quad1_0_00001.txt';


TempFile = readtable(FileName);

% column_data = TempFile.Mx;
% sign_changes = sign(column_data);
% 
% % Find the indexes where the sign changes (crossing zero threshold)
% cross_zero_indexes = find(diff(sign_changes) ~= 0);
% 
% % Display the indexes
% disp('Indexes where values cross the zero threshold:');
% disp(cross_zero_indexes);
% 
% 
% keep_index = true(size(cross_zero_indexes));
% 
% % Iterate through the list to compare each index with the others
% for i = 1:length(cross_zero_indexes)
%     if keep_index(i) % Only process if this index is not already marked for deletion
%         % Compare current index with subsequent indexes
%         for j = i+1:length(cross_zero_indexes)
%             if abs(cross_zero_indexes(j) - cross_zero_indexes(i)) <= 100
%                 % Mark the subsequent index for deletion
%                 keep_index(j) = false;
%             end
%         end
%     end
% end
% 
% % Keep only the indexes that are not marked for deletion
% filtered_indexes = cross_zero_indexes(keep_index);
% 
% % Display the filtered indexes
% disp('Filtered indexes:');
% disp(filtered_indexes);

%FinalTable = TempFile(filtered_indexes(2,1):filtered_indexes(4,1), :);

FinalTable = TempFile(8000:12000, :);
txtFileName = "ShortenedAxialRotQuad1.txt";

writetable(FinalTable, txtFileName, 'WriteVariableNames', true, 'Delimiter', '\t');