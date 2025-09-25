function [curr_prox_coords_f, curr_dist_coords_f] = GetCurrCoords(lig_prox_coords, lig_dist_coords, quad, numLigs)

    % Number of ligaments per quadrant
    q = (numLigs / 4) - 1;
    q1 = q + 1;

    % Starting indices for each quadrant
    startIdx = [1, q1+1, q1+q+2, q1+2*q+3];
    
    % Initialize
    curr_prox_coords = [];
    curr_dist_coords = [];

    % Loop through all 4 quadrants
    for i = 1:4
        if quad == i || quad > 4  % If selected or default (all)
            idx_start = startIdx(i);
            idx_range = idx_start : idx_start + q;
            for j = 0:3
                offset = j * numLigs;
                curr_prox_coords = [curr_prox_coords; lig_prox_coords(idx_range + offset, :)];
                curr_dist_coords = [curr_dist_coords; lig_dist_coords(idx_range + offset, :)];
            end
        end
    end

    % Zero out small values
    curr_prox_coords(abs(curr_prox_coords) < 0.01) = 0;
    curr_dist_coords(abs(curr_dist_coords) < 0.01) = 0;

    % Extract each quadrant’s coordinates
    curr_prox_coords1 = curr_prox_coords(1:q1, :);
    curr_prox_coords2 = curr_prox_coords(q1+1:q1+1+q, :);
    curr_prox_coords3 = curr_prox_coords(q1+q+2:q1+2*q+2, :);
    curr_prox_coords4 = curr_prox_coords(q1+2*q+3:end, :);

    % Store in cell array for convenience
    curr_prox_coords_f = {curr_prox_coords1, curr_prox_coords2, curr_prox_coords3, curr_prox_coords4};

    % Same for distal coords
    curr_dist_coords1 = curr_dist_coords(1:q1, :);
    curr_dist_coords2 = curr_dist_coords(q1+1:q1+1+q, :);
    curr_dist_coords3 = curr_dist_coords(q1+q+2:q1+2*q+2, :);
    curr_dist_coords4 = curr_dist_coords(q1+2*q+3:end, :);

    curr_dist_coords_f = {curr_dist_coords1, curr_dist_coords2, curr_dist_coords3, curr_dist_coords4};


    % figure
    % scatter(lig_dist_coords(:,1), lig_dist_coords(:,2));
    % hold on
    % scatter(lig_prox_coords(:,1), lig_prox_coords(:,2));
    % figure;
    % hold on;
    % 
    % % Plot proximal ligament coordinates
    % scatter3(lig_prox_coords(:, 1), lig_prox_coords(:, 2), lig_prox_coords(:, 3), ...
    %     50, 'b', 'filled', 'DisplayName', 'Proximal Coordinates');
    % 
    % % Plot distal ligament coordinates
    % scatter3(lig_dist_coords(:, 1), lig_dist_coords(:, 2), lig_dist_coords(:, 3), ...
    %     50, 'r', 'filled', 'DisplayName', 'Distal Coordinates');
    % 
    % % Add connections between proximal and distal points
    % for i = 1:size(lig_prox_coords, 1)
    %     plot3([lig_prox_coords(i, 1), lig_dist_coords(i, 1)], ...
    %           [lig_prox_coords(i, 2), lig_dist_coords(i, 2)], ...
    %           [lig_prox_coords(i, 3), lig_dist_coords(i, 3)], ...
    %           'k-', 'LineWidth', 1.5, 'DisplayName', 'Ligament Connection');
    % end
    % 
    % % Add axis labels and legend
    % xlabel('X Coordinate');
    % ylabel('Y Coordinate');
    % zlabel('Z Coordinate');
    % title('3D Plot of Ligament Coordinates');
    % grid on;
    % legend;
    % view(3); % Set 3D view angle
    % hold off;
end