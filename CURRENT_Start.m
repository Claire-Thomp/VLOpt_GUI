% File Org
%file1 = "N:\VLOptCURRENT\InputData\D_Sham1DOFMzIntact.txt"; % Intact file
%file2 = "N:\VLOptCURRENT\InputData\D_Sham1DOFMzCut.txt"; % Cut file
file1 = "/Users/clairethompson/Documents/MATLAB/VLOptCURRENT/InputData/D_ShamDataMxZmmIntact.txt"; % Intact file
file2 = "/Users/clairethompson/Documents/MATLAB/VLOptCURRENT/InputData/D_ShamDataMxZmmCut.txt"; % Cut file

intact = readtable(file1,"VariableNamingRule","preserve");
cut = readtable(file2,"VariableNamingRule","preserve");
nrows = height(intact);
if height(intact) ~= height(cut)
    disp("Note files have different lengths")
    [intact, cut] = FileLengthDiscrep(intact, cut);
        nrows = height(intact);
end

dataIn = [intact.Xmm, intact.Ymm, intact.Zmm, intact.Xdeg, intact.Ydeg, intact.Zdeg];
%FILTER IF NEEDED
% order = 4;
% framelen = 2499;
% intact.Fz = sgolayfilt(intact.Fz,order,framelen);
% cut.Fz = sgolayfilt(cut.Fz,order,framelen);

%VISUALIZE IF NEEDED
% figure
% plot(intact.Mx);
% hold on 
% plot(cut.Mx);
% hold off 
% figure 
% plot(intact.Xdeg);
% hold on 
% plot(cut.Xdeg);
%% Coordinates 

%global params
OuterAPLength = 33.8; % mm
OuterMLLength = 45.8; % mm
InnerAPLength = 26.76; % mm
InnerMLLength = 39.25; % mm
disc_height = 9.4; % mm
OuterLamenaThickness = 0.55; % mm
InnerLamenaThickness = 0.55; % mm
LamenaThickness = 0.55;
numLigs = 12;
flipped = 0; %0 = none, 1 = angled, 2 = straight
coupled = 3; % 0 = no, 1 = yes
numEllipses = 2;
RefPose = [0,0,0,0,0,0]; %1x6 array with the corresponding reference coordinates


%Simple Straight
% ProxCoords = [0,OuterAPLength,disc_height/2;0,-OuterAPLength,disc_height/2;0,OuterMLLength,disc_height/2;0,-OuterMLLength,disc_height/2];
% DistCoords = [0,OuterAPLength,-disc_height/2;0,-OuterAPLength,-disc_height/2;0,OuterMLLength,-disc_height/2;0,-OuterMLLength,-disc_height/2];

if coupled == 0
%Made from ellipse 
    [ProxCoords, DistCoords, Prox_flipped, Dist_flipped] = GetLigCoordinates(OuterAPLength, OuterMLLength, ...
        InnerAPLength, InnerMLLength, disc_height, OuterLamenaThickness, InnerLamenaThickness, numLigs, flipped, numEllipses, coupled);
elseif coupled ==1 
    [ProxCoords, DistCoords] = GetLigCoordinates(OuterAPLength, OuterMLLength, ...
        InnerAPLength, InnerMLLength, disc_height, OuterLamenaThickness, InnerLamenaThickness, numLigs, flipped, numEllipses, coupled);
else 
    [ProxCoords, DistCoords] = GetLigCoordinates(OuterAPLength, OuterMLLength, InnerAPLength, InnerMLLength, disc_height, ...
    LamenaThickness, numLigs, 45);
end

%if sectioning the disc into quadrants 
% quad = 1;
%[curr_prox_coords_f, curr_dist_coords_f] = GetCurrCoords(lig_prox_coords, lig_dist_coords, quad, numLigs);

[outputTable] = InitializeOutputTable(ProxCoords,DistCoords);
%% Ligament Forces and Lengths

WorkingLig = LigamentResultantForces(intact, cut);

[LigLengthTable, trans_results] = InstLigLengthCalc_Final(dataIn, nrows,outputTable, flipped);
RefLength = GetRefLength(LigLengthTable, outputTable);


numRowsToPad = height(WorkingLig) - height(LigLengthTable);
    if numRowsToPad > 0
        rowsToRepeat = LigLengthTable(end-numRowsToPad+1:end, :);
        LigLengthTable = [LigLengthTable; repmat(rowsToRepeat, numRowsToPad, 1)];
    end
%% 

%[x_opt, FinalOutput] = StiffandStrainOptimization_ga( outputTable, quad, RefLength, LigLengthTable, WorkingLig, ProxCoords, DistCoords);
%[x_opt, FinalOutput] = MultiObj_StiffandStrainOptimization_ga( outputTable, quad, RefLength, LigLengthTable, WorkingLig, ProxCoords, DistCoords);
[x_opt, FinalOutput] = New_StiffandStrainOptimization_ga(RefPose, outputTable, dataIn,trans_results, WorkingLig);


% Display results
fprintf('Optimal K: %.4f\n', x_opt(1));
fprintf('Optimal RefStrain: %.4f\n', x_opt(2));

