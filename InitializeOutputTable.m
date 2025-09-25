function [outputTable] = InitializeOutputTable(proxCoords,distCoords)

numPoints = size(proxCoords, 1);
numEllipses = size(proxCoords, 3);

% Preallocate and fill output table
outputTable = table();

outputTable.ProxX = reshape(proxCoords(:, 1, :), [], 1);
outputTable.ProxY = reshape(proxCoords(:, 2, :), [], 1);
outputTable.ProxZ = reshape(proxCoords(:, 3, :), [], 1);

outputTable.DistX = reshape(distCoords(:, 1, :), [], 1);
outputTable.DistY = reshape(distCoords(:, 2, :), [], 1);
outputTable.DistZ = reshape(distCoords(:, 3, :), [], 1);

% Add EllipseID (1 for outermost, up to numEllipses for innermost)
outputTable.EllipseID = repelem((1:numEllipses)', numPoints);
% Compute base quadrant assignment starting from posterior
q = floor(numPoints / 4);  % Points per quadrant

% Create posterior-starting quadrant labels
quadrantPerEllipse = zeros(numPoints, 1);

quadrantPerEllipse(1:q) = 4;           % Posterior Left → now becomes Q4
quadrantPerEllipse(q+1:2*q) = 1;       % Anterior Left   → Q1
quadrantPerEllipse(2*q+1:3*q) = 2;     % Anterior Right  → Q2
quadrantPerEllipse(3*q+1:end) = 3;     % Posterior Right → Q3

% To shift and start from anterior instead:
% Circularly shift the labels by q positions to the left
quadrantPerEllipse = circshift(quadrantPerEllipse, -q);

% Repeat for all ellipses
outputTable.QuadrantID = repmat(quadrantPerEllipse, numEllipses, 1);

% Optional: label as categorical
outputTable.QuadrantID = categorical(outputTable.QuadrantID, 1:4, {'Q1', 'Q2', 'Q3', 'Q4'});
   

   %Used when using quadrants
    % % Start column group index
    % groupIdx = 1;
    % % Loop through the coordinates in groups of 3 (X, Y, Z)
    % while groupIdx <= 4 % Total of 4 groups (ProxX1, ProxX2, etc.)
    %     % Get the current proximal and distal coordinates
    %     % proxCoords = curr_prox_coords_f(:, (groupIdx - 1) * 3 + 1 : groupIdx * 3); % Nx3 matrix
    %     % distCoords = curr_dist_coords_f(:, (groupIdx - 1) * 3 + 1 : groupIdx * 3); % Nx3 matrix
    % 
    %     % Assign each coordinate component to separate columns
    %     outputTable.(['ProxX' num2str(groupIdx)]) = proxCoords(:, 1); % X-coordinates
    %     outputTable.(['ProxY' num2str(groupIdx)]) = proxCoords(:, 2); % Y-coordinates
    %     outputTable.(['ProxZ' num2str(groupIdx)]) = proxCoords(:, 3); % Z-coordinates
    % 
    %     outputTable.(['DistX' num2str(groupIdx)]) = distCoords(:, 1); % X-coordinates
    %     outputTable.(['DistY' num2str(groupIdx)]) = distCoords(:, 2); % Y-coordinates
    %     outputTable.(['DistZ' num2str(groupIdx)]) = distCoords(:, 3); % Z-coordinates
    % 
    %     % Increment group index
    %     groupIdx = groupIdx + 1;
    % end

end