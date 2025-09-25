function [x_opt, FinalOutput] = MultiObj_StiffandStrainOptimization_ga( outputTable, quad, RefLength, LigLengthTable, WorkingLig, ProxCoords, DistCoords)
 %% Multi-Objective Optimization
lb = [0, 0]; % Lower bounds: [K, RefStrain]
ub = [100000, 100]; % Upper bounds: [K, RefStrain]

% Define the multi-objective function
objectiveFcn = @(params) computeMultiObjective(params(1), params(2), RefLength, LigLengthTable, WorkingLig, ProxCoords, DistCoords);

% Options for gamultiobj
options = optimoptions('gamultiobj', ...
    'PopulationSize', 50, ...
    'MaxGenerations', 10000);

% Run the multi-objective optimization
[x_opt, f_opt] = gamultiobj(objectiveFcn, 2, [], [], [], [], lb, ub, [], options);

% Extract optimized values (example: choosing the first Pareto-optimal solution)
K_opt = x_opt(1, 1);
RefStrain_opt = x_opt(1, 2);

% Compute RowSum and ForceSum using optimized parameters
[~, ForceSum1, ForceSum2, ForceSum3, ForceSum4, ForceSum5, ForceSum6,...
    FxExp,FyExp,FzExp,MxExp,MyExp,MzExp] = computeMultiObjective(K_opt, RefStrain_opt, RefLength, LigLengthTable, WorkingLig, ProxCoords, DistCoords);
%% Plot Calculated Forces against experimental values
figure;

% Fx
subplot(3,2,1);
plot(FxExp, '-s', 'LineWidth', 1.5, 'DisplayName', 'FxExp');
hold on;
plot(ForceSum1, '-s', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
xlabel('Index');
ylabel('Force');
title('Fx Comparison');
legend('Location', 'best');
grid on;

% Fy
subplot(3,2,2);
plot(FyExp, '-s', 'LineWidth', 1.5, 'DisplayName', 'FyExp');
hold on;
plot(ForceSum2, '-s', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
xlabel('Index');
ylabel('Force');
title('Fy Comparison');
legend('Location', 'best');
grid on;

% Fz
subplot(3,2,3);
plot(FzExp, '-s', 'LineWidth', 1.5, 'DisplayName', 'FzExp');
hold on;
plot(ForceSum3, '-s', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
xlabel('Index');
ylabel('Force');
title('Fz Comparison');
legend('Location', 'best');
grid on;

% Mx
subplot(3,2,4);
plot(MxExp, '-s', 'LineWidth', 1.5, 'DisplayName', 'MxExp');
hold on;
plot(ForceSum4, '-s', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
xlabel('Index');
ylabel('Moment');
title('Mx Comparison');
legend('Location', 'best');
grid on;

% My
subplot(3,2,5);
plot(MyExp, '-s', 'LineWidth', 1.5, 'DisplayName', 'MyExp');
hold on;
plot(ForceSum5, '-s', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
xlabel('Index');
ylabel('Moment');
title('My Comparison');
legend('Location', 'best');
grid on;

% Mz
subplot(3,2,6);
plot(MzExp, '-s', 'LineWidth', 1.5, 'DisplayName', 'MzExp');
hold on;
plot(ForceSum6, '-s', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
xlabel('Index');
ylabel('Moment');
title('Mz Comparison');
legend('Location', 'best');
grid on;

    outputTable = makeOutput(outputTable, x_opt);
    Sec1 = [outputTable(1:5, 2:8), outputTable(1:5, 1)];
    Sec1.Properties.VariableNames = {'ProxX', 'ProxY', 'ProxZ', 'DistX', 'DistY', 'DistZ','K','RefStrain'};
    Sec1.K = Sec1.K*RefLength(length(RefLength)/4);
    Sec2 = [outputTable(1:5, 9:15),outputTable(1:5, 1)];
    Sec2.Properties.VariableNames = {'ProxX', 'ProxY', 'ProxZ', 'DistX', 'DistY', 'DistZ','K','RefStrain'};
    Sec2.K = Sec2.K*RefLength(2*(length(RefLength)/4));
    Sec3 = [outputTable(1:5, 16:22),outputTable(1:5, 1)];
    Sec3.Properties.VariableNames = {'ProxX', 'ProxY', 'ProxZ', 'DistX', 'DistY', 'DistZ','K','RefStrain'};
    Sec3.K = Sec3.K*RefLength(3*(length(RefLength)/4));
    Sec4 = [outputTable(1:5, 23:29),outputTable(1:5, 1)];
    Sec4.Properties.VariableNames = {'ProxX', 'ProxY', 'ProxZ', 'DistX', 'DistY', 'DistZ','K','RefStrain'};
    Sec4.K = Sec4.K*RefLength(length(RefLength));
    FinalOutput = [Sec1; Sec2; Sec3; Sec4];

%% Multi-Objective Function
    function [f, ForceSum1, ForceSum2, ForceSum3, ForceSum4, ForceSum5, ForceSum6,...
        FxExp,FyExp,FzExp,MxExp,MyExp,MzExp] = computeMultiObjective(K, RefStrain, RefLength, LigLengthTable, WorkingLig, curr_prox_coords, curr_dist_coords)
    
        ligIDs = LigLengthTable.Properties.VariableNames;
        eps1 = 0.03;
        slack_lengths = RefLength ./ (1 + RefStrain / 100);
        FCalc = zeros(height(LigLengthTable), 1);
        Results = zeros(height(LigLengthTable), 6, numel(ligIDs));
    
        for i = 1:numel(ligIDs)
            LigLengths = LigLengthTable{:, ligIDs{i}};
            current_strains = ((LigLengths ./ slack_lengths(i)) - 1);
              % if i <= numel(ligIDs)/4 || i >= numel(ligIDs)*3/4
              %           const = contains(ligIDs{i}, {'E1', 'E2'}) * 1 + contains(ligIDs{i}, {'E3', 'E4'}) * 2.34;
              %    else
              %           const = contains(ligIDs{i}, {'E1', 'E2'}) * 1.35 + contains(ligIDs{i}, {'E3', 'E4'}) * 1.72;
              % end
            const = 1;
              FCalc(:, i) = (K .* const .* (current_strains - eps1) .* (current_strains > (2 * eps1))) .* (current_strains > 0) + ...
                          ((0.25 * K * const) .* ((current_strains .^ 2) ./ eps1) .* (current_strains <= (2 * eps1))) .* (current_strains > 0);
    
            % Forces
            position = [curr_prox_coords(i, :)] - [curr_dist_coords(i, :)];
            unit_position = (position ./ norm(position));
            Results(:, 1, i) = FCalc(:, i) .* unit_position(1); % Fx
            Results(:, 2, i) = FCalc(:, i) .* unit_position(2); % Fy
            Results(:, 3, i) = FCalc(:, i) .* unit_position(3); % Fz
            ForceVec = [Results(:, 1, i), Results(:, 2, i), Results(:, 3, i)];
    
            % Moments
            radius = [0, 0, curr_prox_coords(i, 3)] - [curr_prox_coords(i, 1), curr_prox_coords(i, 2), curr_prox_coords(i, 3)];
            unit_radius = radius ./ vecnorm(radius);
            Radius = repmat(unit_radius, height(FCalc), 1);
    
            for j = 1:height(LigLengthTable)
                Moment = cross(Radius(j, :), ForceVec(j, :));
                Results(j, 4, i) = Moment(:, 1); % Mx
                Results(j, 5, i) = Moment(:, 2); % My
                Results(j, 6, i) = Moment(:, 3); % Mz
            end
        end
    
        % Experimental values
        FxExp = WorkingLig.DeltaFx;
        FyExp = WorkingLig.DeltaFy;
        FzExp = WorkingLig.DeltaFz;
        MxExp = WorkingLig.DeltaMx;
        MyExp = WorkingLig.DeltaMy;
        MzExp = WorkingLig.DeltaMz;
    
        % Compute sum of forces and moments
        ForceSum1 = sum(Results(:, 1, :), 3); % Fx sum
        ForceSum1 = ForceSum1(:);
        ForceSum2 = sum(Results(:, 2, :), 3); % Fy sum
        ForceSum2 = ForceSum2(:);
        ForceSum3 = sum(Results(:, 3, :), 3); % Fz sum
        ForceSum3 = ForceSum3(:);
        ForceSum4 = sum(Results(:, 4, :), 3); % Mx sum
        ForceSum4 = ForceSum4(:);
        ForceSum5 = sum(Results(:, 5, :), 3); % My sum
        ForceSum5 = ForceSum5(:);
        ForceSum6 = sum(Results(:, 6, :), 3); % Mz sum
        ForceSum6 = ForceSum6(:);
    
        % Compute RMSE for Mx and Fz
        rmse_Fx = rmse(ForceSum1, FxExp);
        rmse_Fy = rmse(ForceSum2, FyExp);
        rmse_Fz = rmse(ForceSum3, FzExp);
        rmse_Mx = rmse(ForceSum4, MxExp);
        rmse_My = rmse(ForceSum5, MyExp);
        rmse_Mz = rmse(ForceSum6, MzExp);
    
        % Multi-objective function (equally weighted)
        f = [rmse_Mz, rmse_Fz];
    
    end
end

