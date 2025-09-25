function [outputTable] = makeOutput(outputTable, x)
    
    outputTable = [outputTable, array2table(zeros(height(outputTable), 2))];
    
    K_opt = x(1);
    RefStrain_opt = x(2);
    consts = [1,2.34,1.35,1.72];
    
    K1 = K_opt*consts(1);
    K2 = K_opt*consts(2);
    K3 = K_opt*consts(3);
    K4 = K_opt*consts(4);


    

    % for i = 1:height(outputTable)
    %     if i < (height(outputTable)/4)
    %         outputTable(i,7) = table(K1);
    %     elseif i > (height(outputTable)/4) && i < (height(outputTable)/2)
    %         outputTable(i,7) = table(K2);
    %     elseif i > (height(outputTable)/2) && i < (3*(height(outputTable))/2)
    %         outputTable(i,7) = table(K3);
    %     else 
    %         outputTable(i,7) = table(K4);
    %     end
    % end


    outputTable(:,8) = RefStrain_opt;

  %Used when divided into quadrants 
   % repeatingValue1 = K_opt*consts(1);
   % repeatingValue2 = K_opt*consts(2);
   % repeatingValue3 = K_opt*consts(3);
   % repeatingValue4 = K_opt*consts(4);
   % numRows = height(outputTable);
   % matrixData1 = repmat({repeatingValue1}, numRows, 1);
   % matrixData2 = repmat({repeatingValue2}, numRows, 1);
   % matrixData3 = repmat({repeatingValue3}, numRows, 1);
   % matrixData4 = repmat({repeatingValue4}, numRows, 1);
   % matrixData5 = repmat({RefStrain_opt}, numRows,1);
   % myTable1 = cell2table(matrixData1, 'VariableNames', {'K_1'});
   % myTable2 = cell2table(matrixData2, 'VariableNames', {'K_2'});
   % myTable3 = cell2table(matrixData3, 'VariableNames', {'K_3'});
   % myTable4 = cell2table(matrixData4, 'VariableNames', {'K_4'});
   % myTable5 = cell2table(matrixData5, 'VariableNames', {'RefStrain'});
   % 
   % insertIdx1 = 7;
   % leftTable1 = outputTable(:, 1:insertIdx1-1);
   % rightTable1 = outputTable(:, insertIdx1:end);
   % outputTable = [leftTable1, myTable1, rightTable1];
   % 
   % insertIdx2 = 14;
   % leftTable2 = outputTable(:, 1:insertIdx2-1);
   % rightTable2 = outputTable(:, insertIdx2:end);
   % outputTable = [leftTable2, myTable2, rightTable2];
   % 
   % insertIdx3 = 21;
   % leftTable3 = outputTable(:, 1:insertIdx3-1);
   % rightTable3 = outputTable(:, insertIdx3:end);
   % outputTable = [leftTable3, myTable3, rightTable3];
   % 
   % insertIdx4 = 28;
   % leftTable4 = outputTable(:, 1:insertIdx4-1);
   % rightTable4 = outputTable(:, insertIdx4:end);
   % outputTable = [myTable5,leftTable4, myTable4, rightTable4];


end