
    tic
    file1 = "N:\VLOptTest3\ShortenedintactTension.txt"; % Intact file
    file2 = "N:\VLOptTest3\ShortenedcutTension.txt"; % Cut file

    % Constants
    OuterAPLength = 35.9; % mm
    OuterMLLength = 49.75; % mm
    InnerAPLength = 26.76; % mm
    InnerMLLength = 39.25; % mm
    disc_height = 8.34; % mm
    OuterLamenaThickness = 0.55; % mm
    InnerLamenaThickness = 0.55; % mm
    % RefLengthAnterior = 7.7; % mm
    % RefLengthLateral = 7.7;  
    % RefLengthPosterior = 7.7;
    quad = 1;
    numLigs = 12;
   
    
    % Read input data
    intact = readtable(file1);
    cut = readtable(file2);
    order = 4;
    framelen = 2999;
    intact.Fz = sgolayfilt(intact.Fz,order,framelen);
    cut.Fz = sgolayfilt(cut.Fz,order,framelen);

    plot(intact.Fz);
    hold on
    plot(cut.Fz);
    
    % Preprocessing
    [intact, cutLig] = FileLengthDiscrep(intact, cut);
    nrows = height(intact);
    
    % Get ligament coordinates
    [lig_prox_coords, lig_dist_coords] = GetLigCoordinates(OuterAPLength, OuterMLLength, ...
        InnerAPLength, InnerMLLength, disc_height, OuterLamenaThickness, InnerLamenaThickness, numLigs);
    
    [curr_prox_coords_f,curr_dist_coords_f, curr_prox_coords1,  curr_prox_coords2,  curr_prox_coords3,...  
        curr_prox_coords4,  curr_dist_coords1,  curr_dist_coords2,  curr_dist_coords3,  curr_dist_coords4] = GetCurrCoords(lig_prox_coords, lig_dist_coords,quad, numLigs);
    
    prox_vert = [curr_prox_coords1;curr_prox_coords2;curr_prox_coords3;curr_prox_coords4];
    dist_vert = [curr_dist_coords1;curr_dist_coords2;curr_dist_coords3;curr_dist_coords4];

    % Create the output table for ligament coordinates
    [outputTable] = InitializeOutputTable(curr_prox_coords_f,curr_dist_coords_f);
    
    WorkingLig = LigamentResultantForces(intact, cutLig);
    
    % Compute ligament lengths for all quadrants
    LigLengthTable = [];
    
    curr_prox_coords_cell = {curr_prox_coords1, curr_prox_coords2, curr_prox_coords3, curr_prox_coords4}; 
    curr_dist_coords_cell = {curr_dist_coords1, curr_dist_coords2, curr_dist_coords3, curr_dist_coords4}; 
    
    for i = 1:4
        LigLengthTable = [LigLengthTable, InstLigLengthCalc_Final(intact, nrows, curr_prox_coords_cell{i}, curr_dist_coords_cell{i}, char('a' + i - 1))];
    end 
    RefLength = GetRefLength(LigLengthTable, outputTable, quad);
    
    % Pad LigLengthTable to match WorkingLig size
    numRowsToPad = height(WorkingLig) - height(LigLengthTable);
    if numRowsToPad > 0
        rowsToRepeat = LigLengthTable(end-numRowsToPad+1:end, :);
        LigLengthTable = [LigLengthTable; repmat(rowsToRepeat, numRowsToPad, 1)];
    end
    
    [x_opt] = StiffandStrainOptimization_ga(OuterAPLength, InnerAPLength,quad, RefLength, LigLengthTable, WorkingLig,prox_vert, dist_vert);
    outputTable = makeOutput(outputTable, x_opt);
%% 
    % Sec1 = [outputTable(1:5, 2:8), outputTable(1:5, 1)];
    % Sec1.Properties.VariableNames = {'ProxX', 'ProxY', 'ProxZ', 'DistX', 'DistY', 'DistZ','K','RefStrain'};
    % Sec2 = [outputTable(1:5, 9:15),outputTable(1:5, 1)];
    % Sec2.Properties.VariableNames = {'ProxX', 'ProxY', 'ProxZ', 'DistX', 'DistY', 'DistZ','K','RefStrain'};
    % Sec3 = [outputTable(1:5, 16:22),outputTable(1:5, 1)];
    % Sec3.Properties.VariableNames = {'ProxX', 'ProxY', 'ProxZ', 'DistX', 'DistY', 'DistZ','K','RefStrain'};
    % Sec4 = [outputTable(1:5, 23:29),outputTable(1:5, 1)];
    % Sec4.Properties.VariableNames = {'ProxX', 'ProxY', 'ProxZ', 'DistX', 'DistY', 'DistZ','K','RefStrain'};
    % FinalOutput = [Sec1; Sec2; Sec3; Sec4];
    toc
%end
