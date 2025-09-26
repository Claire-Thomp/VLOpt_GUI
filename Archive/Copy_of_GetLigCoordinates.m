function [final_prox_coord, final_dist_coord, final_prox_coord_flipped, final_dist_coord_flipped] = GetLigCoordinates(OuterAPLength, OuterMLLength, InnerAPLength, InnerMLLength, disc_height, ...
    OuterLamenaThickness, InnerLamenaThickness, numPoints, flipped, numEllipses, coupled)

ellipseParams = [
    OuterAPLength/2,OuterMLLength/2,1,1;   % Outermost
    (OuterAPLength/2)-OuterLamenaThickness,(OuterMLLength/2)-OuterLamenaThickness, -1, -1;  % Second outer
    (InnerAPLength/2)+InnerLamenaThickness,(InnerMLLength/2)+InnerLamenaThickness, -1, -1   % Second inner;
    InnerAPLength/2,InnerMLLength/2,1,1   % Innermost
];

dim = 3;  % Assuming each coordinate is 3D

% Preallocate arrays
lig_prox_coord = zeros(numPoints, dim, numEllipses);
lig_dist_coord = zeros(numPoints, dim, numEllipses);

% Loop through ellipses
for i = 1:numEllipses
    Halfa = ellipseParams(i,1);
    Halfb = ellipseParams(i,2);
    angle = ellipseParams(i,3);
    direction = ellipseParams(i,4);
    
    % Proximal points
    prox = get_prox_point(Halfa, Halfb, numPoints);
    prox(abs(prox) < 0.01) = 0;
    prox(end,:) = [];  % Remove duplicate or final point
    lig_prox_coord(:,:,i) = prox;
    
    % Distal points
    dist = get_dist_point(disc_height, prox, angle, direction, Halfa, Halfb, numPoints);
    dist(abs(dist) < 0.01) = 0;
    dist(end,:) = [];  % Remove final point
    lig_dist_coord(:,:,i) = dist;
end

if flipped == 0
    final_prox_coord = lig_prox_coord;
    final_dist_coord = lig_dist_coord; 
    final_prox_coord_flipped = zeros(numPoints, dim, numEllipses);
    final_dist_coord_flipped = zeros(numPoints, dim, numEllipses);
elseif flipped == 1 % flipped ligs tracking the same as originals 
    
    %Flip all Z values
    flipped_prox = [lig_prox_coord(:,1,:), lig_prox_coord(:,2,:), -lig_prox_coord(:,3,:)];
    flipped_dist = [lig_dist_coord(:,1,:), lig_dist_coord(:,2,:), -lig_dist_coord(:,3,:)];
    if coupled == 1
        %Preallocate final matrices
        n = size(lig_prox_coord, 1);
        final_prox_coord = zeros(2*numPoints, dim,numEllipses);
        final_dist_coord = zeros(2*numPoints, dim, numEllipses);
        final_prox_coord_flipped = zeros(numPoints, dim, numEllipses);
        final_dist_coord_flipped = zeros(numPoints, dim, numEllipses);
        
        % Interleave original and flipped points
        final_prox_coord(1:2:end, :,:) = lig_prox_coord;
        final_prox_coord(2:2:end, :,:) = flipped_prox;
        
        final_dist_coord(1:2:end, :,:) = lig_dist_coord;
        final_dist_coord(2:2:end, :,:) = flipped_dist;
    else
        final_prox_coord = lig_prox_coord;
        final_dist_coord = lig_dist_coord; 
        final_prox_coord_flipped = flipped_prox;
        final_dist_coord_flipped = flipped_dist;
    end
else
    %Straight Flipped Ligs
    flipped_prox = [lig_prox_coord(:,1,:), lig_prox_coord(:,2,:), -lig_prox_coord(:,3,:)];
    flipped_dist = [lig_prox_coord(:,1,:), lig_prox_coord(:,2,:), lig_prox_coord(:,3,:)];
   
    n = size(lig_prox_coord, 1);
    final_prox_coord = zeros(2*numPoints, dim,numEllipses);
    final_dist_coord = zeros(2*numPoints, dim, numEllipses);
    
    % Interleave original and flipped points
    final_prox_coord(1:2:end, :,:) = lig_prox_coord;
    final_prox_coord(2:2:end, :,:) = flipped_prox;
    
    final_dist_coord(1:2:end, :,:) = lig_dist_coord;
    final_dist_coord(2:2:end, :,:) = flipped_dist;
end  
    
    function [lig_prox_coord] = get_prox_point(APLength, MLLength, numPoints)
      
        %equation of the line 
        Halfa = APLength;
        Halfb = MLLength;
        %must be odd
        numPoints = numPoints+1;
        arcLength = @(theta) sqrt((Halfa * sin(theta)).^2 + (Halfb * cos(theta)).^2);
        theta = linspace(pi/2, 2.5*pi, 1000); % Dense set of angles for high accuracy
        dArc = arcLength(theta);
        cumArc = cumtrapz(theta, dArc); % Cumulative arc length
        totalArc = cumArc(end);
        evenArc = linspace(0, totalArc, numPoints);
        evenTheta = interp1(cumArc, theta, evenArc); 

        % Generate x and y coordinates
        proxx = Halfa * cos(evenTheta);
        proxy = Halfb * sin(evenTheta);
        proxz = ones(size(proxx)) * (disc_height / 2); % z-coordinates

        lig_prox_coord = [proxx', proxy', proxz'];
        
        
        % scatter(lig_prox_coord(:, 1),lig_prox_coord(:,2));
        % hold on
    end

    function [lig_dist_points] = get_dist_point(disc_height, lig_prox_coord, multiplier, middle,MLLength,APLength, numPoints)
            
          
        b= APLength;
        a = MLLength;
        
        if middle == 1
            angle = deg2rad(45); %rad
        else 
            angle = -deg2rad(45);
        end
        

        
            x0 = lig_prox_coord(1,1);
            y0 = lig_prox_coord(1,2);
            x = disc_height*tan(angle);
            [~, ~, delta_t]=move_on_ellipse(a, b, x0, y0, x, multiplier);
            
            
            numPoints = numPoints+1;
            arcLength = @(theta) sqrt((a * sin(theta)).^2 + (b * cos(theta)).^2);
            theta = linspace(pi/2+delta_t, 2.5*pi+delta_t, 1000); % Dense set of angles for high accuracy
            dArc = arcLength(theta);
            cumArc = cumtrapz(theta, dArc); % Cumulative arc length
            totalArc = cumArc(end);
            evenArc = linspace(0, totalArc, numPoints);
            evenTheta = interp1(cumArc, theta, evenArc); 
    
            % Generate x and y coordinates
            distx = a * cos(evenTheta);
            disty = b * sin(evenTheta);
            
            distz = ones(size(distx)) * -(disc_height / 2);
            lig_dist_points = [distx', disty', distz'];
      
        
        
        % scatter(lig_dist_points(:,1),lig_dist_points(:,2));
    end



    function [x_new, y_new, delta_t] = move_on_ellipse(a, b, x0, y0, x, direction)
        %Ensure the initial point is on the ellipse
   
    
        % Compute initial parameter t0
        t0 = atan2(y0 / b, x0 / a);
    
        % Compute radius of curvature at (x0, y0)
        R = (a^2 * b^2) / ((a^2 * sin(t0)^2 + b^2 * cos(t0)^2)^(3/2));
    
        % Compute angular step
        delta_t = x / R; 
    
        % Adjust the angle based on direction
        if direction == 1 
            t_new = t0 - delta_t;
        elseif direction == -1
            t_new = t0 + delta_t;
        end
    
        % Compute new coordinates
        x_new = a * cos(t_new);
        y_new = b * sin(t_new);
    end
    % % 3D Visualization of ligament coordinates
    % figure;
    % hold on;
    % grid on;
    % axis equal;
    % view(3);
    % 
    % % Plot proximal points (top)
    % scatter3(final_prox_coord(:,1), final_prox_coord(:,2), final_prox_coord(:,3), ...
    %     50, 'b', 'filled', 'DisplayName', 'Proximal Points');
    % 
    % % Plot distal points (bottom)
    % scatter3(final_dist_coord(:,1), final_dist_coord(:,2), final_dist_coord(:,3), ...
    %     50, 'r', 'filled', 'DisplayName', 'Distal Points');
    % 
    % % Connect proximal and distal points with lines
    % for i = 1:size(final_prox_coord,1)
    %     line([final_prox_coord(i,1), final_dist_coord(i,1)], ...
    %          [final_prox_coord(i,2), final_dist_coord(i,2)], ...
    %          [final_prox_coord(i,3), final_dist_coord(i,3)], ...
    %          'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    % end
    % 
    % xlabel('X (mm)');
    % ylabel('Y (mm)');
    % zlabel('Z (mm)');
    % title('3D Visualization of Ligament Coordinates');
    % hold off;

end