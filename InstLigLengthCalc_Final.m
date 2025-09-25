function [LigLengthTable, trans_results] = InstLigLengthCalc_Final(dataIn, nrows, outputTable, flipped)
    
    lig_prox_coord = table2array(outputTable(:, 1:3));
    lig_dist_coord = table2array(outputTable(:, 4:6));

    colTwo = [-0.2574, -0.0222, -28.7181, -0.5773, -0.3933, 0.7908]';
    colThree = [0, 0, 0, 0, 0, 0]';
    vdat = [];
    for i = 1:(height(dataIn))
        dataTemp = dataIn(i, 1:6);
        dataTempT = dataTemp';
        dataTempT = [dataTempT, colTwo, colThree];
        vdat = [vdat; dataTempT];
    end

    numCols = size(lig_prox_coord, 1);      % Total number of ligaments
    numEllipses = 4;
    numPoints = numCols / (numEllipses * 2);  % Original + flipped per ellipse

    ligamentNames = {};
    
        
    for i = 1:numCols
        
        if flipped ~= 0 
            iseven = rem(i, 2) == 0;
            if iseven == 0
            flipTag = 'O';
            else
            flipTag = 'F';
            end
        else 
           flipTag = 'O';
        end 
    
        if i <= numCols / 4
            ellipseNum = ceil(1);
        elseif i > numCols / 4 && i <= numCols / 2
            ellipseNum = ceil(2);
        elseif i > numCols / 2 && i <= 3*numCols / 4
            ellipseNum = ceil(3);
        else
            ellipseNum = ceil(4);
        end
            
        name = sprintf('Lig_E%d_%s_%d', ellipseNum, flipTag,i);
        ligamentNames{end+1} = name;
    end

    LigLengthTable = table('Size', [nrows, numCols], ...
                       'VariableTypes', repmat({'double'}, 1, numCols), ...
                       'VariableNames', ligamentNames);

    
    rowIdx = 1;
    for i = 1:6:height(vdat)
        pose = vdat(i:i+5, 1);
        Roffset = vdat(4:6, 2) + [20 0 0]';
        T = [0 0 0];
        [T_fem, T_tib] = VIVO_Positions_rotmatrix(pose(4), pose(5), pose(6), pose(2), pose(1), pose(3), Roffset, T, 'r');
    
        fem_trans = T_fem(1:3, 1);
        fem_rot = T_fem(1:3, 2:4) - fem_trans;
        fem_transformation = [fem_rot fem_trans];
    
        tib_trans = T_tib(1:3, 1);
        tib_rot = T_tib(1:3, 2:4) - tib_trans;
        tib_transformation = [tib_rot tib_trans];
        
        [flexion_transformation, gimbal_transformation, platen_transformation] = ActuatorOffsets(vdat, fem_transformation, tib_transformation);
    
        trans_results{1} = fem_transformation;
        trans_results{2} = tib_transformation;
        trans_results{3} = flexion_transformation;
        trans_results{4} = gimbal_transformation;
        trans_results{5} = platen_transformation;
        
        % Compute current ligament insertion positions
        F = ((trans_results{1}(1:3, 1:3)) * lig_prox_coord' + (trans_results{1}(1:3, 4)))';
        T = ((trans_results{2}(1:3, 1:3)) * lig_dist_coord' + (trans_results{2}(1:3, 4)))';
        
        % Compute ligament vectors and lengths
        lig_vect = F - T;
        lig_lengths = vecnorm(lig_vect, 2, 2); % Compute lengths along rows
        % Store ligament lengths in the table
        LigLengthTable(rowIdx, :) = array2table(lig_lengths', 'VariableNames', LigLengthTable.Properties.VariableNames);
        rowIdx = rowIdx + 1;
    end
    
end
