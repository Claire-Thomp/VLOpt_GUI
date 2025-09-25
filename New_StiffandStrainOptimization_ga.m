function [x_opt, FinalOutput] = New_StiffandStrainOptimization_ga(refPose, outputTable, dataIn, trans_results, WorkingLig)

%Variable Preparation
ligModel = table2array(outputTable(:, 1:6));
numLigs = length(ligModel);   
testPoses = num2cell(dataIn);
numTestPoses  = length(dataIn);
Ttrans = trans_results{2};
Ftrans = trans_results{1};


%adjust to interanl standard of VIVO 
ligModelN = cell(size(ligModel,1),1);
ligcoordscale = [0.001 0.001 0.001 0.001 0.001 0.001];
for lig = 1:numLigs
    ligRow = ligModel(lig,:);
    ligModelN{lig} = ligRow .* ligcoordscale;
end

testPosesN = cell(size(dataIn, 1),1);
coordScale = [0.001 0.001 0.001 pi/180 pi/180 pi/180];
coordOffset = [0 0 0 0 pi/2 0];
for pose = 1:numTestPoses
    testPosesN{pose} = (testPoses{pose} .* coordScale) + coordOffset;
end

refPoseN = (refPose .* coordScale) + coordOffset;

%remap distal coords of ligs to T-ref coords 
refCell = num2cell(refPoseN);
tZfRef = Ttrans*refCell{:};
for lig = 1:numLigs
    dlig = [1; (ligModelN{lig,4:6})'];%to col vec 
    dlig = tZfRef * dlig;
    ligModelN{lig,4:6} = (dlig(2:4))'; %back to row
end

%output arrays
LigFinal = cell(numTestPoses, numLigs);
NetFinal = cell(numTestPoses,1);

%rotational sign correction (alpha, beta are negative --> gamma is negative
%when chirality is right) 
rotsign = [-1 -1 1];


ParamsTable = table(size(height(outputTable),1),'VariableNames',"Strain Const");
for i = 1:height(outputTable) 
    quadrants = string(outputTable{:,8});
    ventral = contains(quadrants(i),"Q1") || contains(quadrants(i),"Q4");
    dorsal = contains(quadrants(i),"Q2") || contains(quadrants(i),"Q3");
    %Ventral interior stiffness has a const of 1, these ratios are
    %built from the Hozapfel paper 
     if outputTable.EllipseID(i) == 1 &&  ventral == 1 
            const = 2.34; %Ventral Exterior
     elseif outputTable.EllipseID(i) == 1 && dorsal == 1
            const = 1.72/0.99; %Dorsal Exterior
     elseif outputTable.EllipseID(i) == 2 && dorsal == 1
            const = 1/0.99; %Dorsal Interior
     else
            const = 1; %Ventral Interior
     end

     ParamsTable(i,1) = table(const);

end

     outputTable= [outputTable,ParamsTable];

% Optimization
    lb = [0, 0]; % Lower bounds: [K, RefStrain]
    ub = [100000, 100]; % Upper bounds: [K, RefStrain]

    % Define the objective function
    objectiveFcn = @(params) computeObjective1(params(1), params(2), umTestPoses,numLigs,testPosesN,Ttrans,ligModelN,rotsign,NetFinal, LigFinal, ParamsTable);

    % Options for ga
    options = optimoptions('ga', ...
        'PopulationSize', 100, ...
        'MaxGenerations', 10000);

    % Run the optimization
    [x_opt] = ga(objectiveFcn, 2, [], [], [], [], lb, ub, [], options);

    % After optimization
    K_opt = x_opt(1);
    RefStrain_opt = x_opt(2);

    % Compute RowSum using the optimized parameters
    [~, RowSum, MxExp, Results] = computeObjective1(K_opt, RefStrain_opt, ...  %change here
        umTestPoses,numLigs,testPosesN,Ttrans,ligModelN,rotsign,NetFinal, LigFinal, ParamsTable);

    % Plot RowSum against WorkingLig
    figure;
    hold on;

    % Plot RowSum
    plot(MxExp, '-s', 'LineWidth', 1.5, 'DisplayName', 'Mx Experimental'); %change here
    hold on
    % Plot DeltaMx
    plot(RowSum, '-s', 'LineWidth', 1.5, 'DisplayName', 'Mx Simulated');

    % Add labels, legend, and grid
    xlabel('Index');
    ylabel('Values');
    title('Comparison of Simulated and Experimental Data');
    legend('Location', 'best');
    grid on;

    hold off;


% Objective function
function [f1,RowSum, MxExp, Results] = computeObjective1(K, RefStrain, ...  %change here
        numTestPoses,numLigs,testPosesN,Ttrans,ligModelN,rotsign,NetFinal, LigFinal, ParamsTable)
    

    for i = 1:numTestPoses
        tcoords = num2cell(testPosesN{i});
        tZfTC = Ttrans*(tcoords{:}); %transform test coords
        
        %T to Q coord conversion
        [beta gamma] = tcoords{5:6}; %get flexion and abduction angles
        qXt = [cos(gamma)*sin(beta) sin(beta)*sin(gamma)    cos(beta);
                -sin(gamma) cos(gamma)  0;
                0   0   -1];
        NetFinal{i} =zeros(1,6);

        %forces time!
        for ilig = 1:numLigs
            ligModelN{ilig} = [ligModelN{ilig}, K];
            ligModelN{ilig} = [ligModelN{ilig}, RefStrain];
            ligModelN{ilig} = [ligModelN{ilig}, cptL0(ligModelN{ilig})]; %this finds the L0 and adds it to the ligModel
            StiffMultiplier = ParamsTable{ilig};
            tforce = LForce(ligModelN{ilig}, tZfTC, StiffMultiplier); %force and moment in T frame 
            LigFinal{i, ilig} = -1*[(qXt * tforce(1:3)')', (qXt * tforce(4:6)')' .* rotsign]; %not sure about the negative 1 here
            LigFinal{i} = NetFinal{i} + LigFinal{i, ilig}; % accumulate resultant
        end
    end

    function L0 = cptL0(lig)
        % return L0 given ligament insertion site positions at reference and reference strain
        P = lig(1:3);
        D = lig(4:6);
        dlig = P - D;
        lrf = sqrt(dlig * dlig');
        L0 = lrf/(lig(8)+1);  %lig(8) = reference strain
    end
    
    function fT = LForce(lig, tZf, StiffMultiplier)
        % lig is a [1 x 9] specifying the ligament: [Px Py Pz Dx Dy Dz K eRf L0]
        % The distal insertion site has already been statically transformed to the T frame
        % The proximal point is in F and must be converted to T using the dynamic xform passed as tZf.
        P = lig(1:3); % proximal (in F)
        D = lig(4:6); % distal (in T)
        K = lig(7)*StiffMultiplier;
        eRf = lig(8) % unused
        L0 = lig(9);
        
        Pt = tZf * [1; P'];
        Pt = Pt(2:4)'; % proximal insertion site referred to T
        dlig = Pt - D; % vector from distal to proximal
        llig = sqrt( dlig * dlig'); % current length of ligament
        if llig < L0   % detect ligament slack
	        fT = zeros(1,6);
	        return
        end
        nlig = dlig / llig; % unit vector aligned with distal to proximal line
		        %; gives direction of ligament force wrt T
        
        %compute force using Blankevoort model
        e = (llig-L0) / L0;
        esp1 = 0.03;
        if e < 2*esp1
	        flig = ( K * .25 * e * e / esp1) * nlig;
        else
	        flig = ( K * (e-esp1) ) * nlig;
        end
        mlig = cross(D, flig); % moment wrt T origin
        fT = [flig, mlig];
    end

        MxExp = WorkingLig.DeltaMx; %change here
        RowSum = sum(Results(:,4,:), 3);  %change here
        RowSum = RowSum(:);
    % Calculate RMSE
      f1 = rmse(RowSum, MxExp);  %change here
    end

end


