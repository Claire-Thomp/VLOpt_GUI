function [RefLength] = GetRefLength(LigLengthTable,outputTable)
    

% for i = 1:width(LigLengthTable)
%     RefLength(i,1) = LigLengthTable(1,i);
% end

ProxCoords = table2array(outputTable(:, 1:3));
DistCoords = table2array(outputTable(:, 4:6));
for i = 1:height(outputTable)
    RefLength(i,1) = sqrt(((DistCoords(i,1))-(ProxCoords(i,1)))^2 + ((DistCoords(i,2))-(ProxCoords(i,2)))^2 + ((DistCoords(i,3))-(ProxCoords(i,3)))^2);
end

% RefLength = table2array(RefLength);