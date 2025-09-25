function [x_opt, FinalOutput] = StiffandStrainOptimization_ga( outputTable, quad, RefLength, LigLengthTable, WorkingLig, ProxCoords, DistCoords)
 % Optimization
    lb = [0, 0]; % Lower bounds: [K, RefStrain]
    ub = [100000, 100]; % Upper bounds: [K, RefStrain]

    % Define the objective function
    objectiveFcn = @(params) computeObjective1(params(1), params(2), RefLength, LigLengthTable, WorkingLig,ProxCoords, DistCoords,quad);

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
        RefLength, LigLengthTable, WorkingLig, ProxCoords, DistCoords,quad);

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

    StiffnessTable = table();
    ligIDs = LigLengthTable.Properties.VariableNames;
    for i = 1:numel(ligIDs)
        
         if i <= numel(ligIDs)/4 || i >= numel(ligIDs)*3/4
                const = contains(ligIDs{i}, {'E1', 'E2'}) * 1 + contains(ligIDs{i}, {'E3', 'E4'}) * 2.34;
         else
                const = contains(ligIDs{i}, {'E1', 'E2'}) * 1.35 + contains(ligIDs{i}, {'E3', 'E4'}) * 1.72;
         end
    
         StiffnessTable(i,1) = table(K_opt*const);
         StiffnessTable(i,2) = table(RefStrain_opt);
    
    end

     FinalOutput= [outputTable,StiffnessTable];

% Objective function
function [f1,RowSum, MxExp, Results] = computeObjective1(K, RefStrain, ...  %change here
    RefLength, LigLengthTable, WorkingLig, curr_prox_coords, curr_dist_coords,quad)

    ligIDs = LigLengthTable.Properties.VariableNames;
    eps1 = 0.03;
    slack_lengths = RefLength ./ (1 + RefStrain / 100);
    FCalc = zeros(height(LigLengthTable), 1);
    Results = zeros(height(LigLengthTable),6,numel(ligIDs));

    for i = 1:numel(ligIDs)
        LigLengths = LigLengthTable{:, ligIDs{i}};
        current_strains = ((LigLengths ./ slack_lengths(i)) - 1);

     if i <= numel(ligIDs)/4 || i >= numel(ligIDs)*3/4
            const = contains(ligIDs{i}, {'E1', 'E2'}) * 1 + contains(ligIDs{i}, {'E3', 'E4'}) * 2.34;
     else
            const = contains(ligIDs{i}, {'E1', 'E2'}) * 1.35 + contains(ligIDs{i}, {'E3', 'E4'}) * 1.72;
     end

     FCalc(:, i) = (K .* const .* (current_strains - eps1) .* (current_strains > (2 * eps1))) .* (current_strains > 0) + ...
                      ((0.25 * K * const) .* ((current_strains .^ 2) ./ eps1) .* (current_strains <= (2 * eps1))) .* (current_strains > 0);

        %forces
        position = [curr_prox_coords(i,:)]-[curr_dist_coords(i,:)];
        unit_position = (position./norm(position));
        Results(:,1,i) = FCalc(:,i).*unit_position(1); %Fx
        Results(:,2,i) = FCalc(:,i).*unit_position(2); %Fy
        Results(:,3,i) = FCalc(:,i).*unit_position(3); %Fz
        ForceVec = [Results(:,1,i),Results(:,2,i),Results(:,3,i)];

        %moments
        radius = [0,0,curr_prox_coords(i,3)] - [curr_prox_coords(i,1),curr_prox_coords(i,2),curr_prox_coords(i,3)];
        unit_radius = radius ./ vecnorm(radius);
        Radius = repmat(unit_radius,height(FCalc),1); %Radius

        for j = 1:height(LigLengthTable)
            Moment = cross(Radius(j,:), ForceVec(j,:));
            Results(j,4,i) = Moment(:,1); %Mx
            Results(j,5,i) = Moment(:,2); %My
            Results(j,6,i) = Moment(:,3); %Mz
        end


    end
        MxExp = WorkingLig.DeltaMx; %change here
        RowSum = sum(Results(:,4,:), 3);  %change here
        RowSum = RowSum(:);
    % Calculate RMSE
        f1 = rmse(RowSum, MxExp);  %change here
    end

    % % Visualization function for K and RefStrain
    % function visualizeObjective(OuterAPLength, InnerAPLength, quad, RefLength, LigLengthTable, WorkingLig)
    %     % Define ranges for K and RefStrain
    %     K_vals = linspace(0, 1000, 50); % Adjust the range and resolution as needed
    %     RefStrain_vals = linspace(0, 100, 50); % Match the adjusted upper bound
    %     [K_grid, RefStrain_grid] = meshgrid(K_vals, RefStrain_vals);
    % 
    %     % Compute objective function values
    %     obj_vals = arrayfun(@(K, RefStrain) computeObjective1(OuterAPLength, InnerAPLength, quad, K, RefStrain, RefLength, LigLengthTable, WorkingLig,curr_prox_coords, curr_dist_coords), K_grid, RefStrain_grid);
    % 
    %     % Plot the objective function
    %     figure;
    %     surf(K_grid, RefStrain_grid, obj_vals, 'EdgeColor', 'none');
    %     colorbar;
    %     title('Objective Function Landscape');
    %     xlabel('K');
    %     ylabel('RefStrain');
    %     zlabel('Objective Value');
    %     view(135, 30); % Adjust the viewing angle for better visualization
    %     grid on;
    % end

end


