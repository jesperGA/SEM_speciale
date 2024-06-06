function [Object_nodes, scale] = insertObject(mesh, study)
    % insertObject Inserts an object into the mesh for duct flow simulations
    % Inputs:
    %   mesh - Struct containing the mesh details
    %   study - Struct containing the study parameters
    % Outputs:
    %   Object_nodes - Indices of the nodes within the object
    %   scale - Scaling factors for the penalty method

    LX = max(mesh.X(:, 2)) - min(mesh.X(:, 2));
    LY = max(mesh.X(:, 3)) - min(mesh.X(:, 3));

    if strcmp(study.example, 'Duct')
        if length(study.ObjectCoords) == 4 % Square object
            Object_nodes = mesh.X( ...
                mesh.X(:, 2) >= study.ObjectCoords(1) & mesh.X(:, 2) <= study.ObjectCoords(2) & ...
                mesh.X(:, 3) >= study.ObjectCoords(3) & mesh.X(:, 3) <= study.ObjectCoords(4), ...
                1);
        elseif length(study.ObjectCoords) == 3 % Circular object
            % Extract circle parameters
            centerX = study.ObjectCoords(1);
            centerY = study.ObjectCoords(2);
            radius = study.ObjectCoords(3);

            % Calculate the distance from each node to the center of the circle
            distances = sqrt((mesh.X(:, 2) - centerX).^2 + (mesh.X(:, 3) - centerY).^2);

            % Find nodes where the distance is less than or equal to the radius
            insideCircle = distances <= radius;
            Object_nodes = mesh.X(insideCircle, 1);

            % Scaling factor decreases linearly from 1 at the center to 0 at the perimeter
            scale = 2 - (2 * distances(insideCircle) / radius);
            scale(scale > 1) = 1;
        end
    end
end
