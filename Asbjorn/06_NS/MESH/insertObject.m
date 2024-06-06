function [Object_nodes, scale] = insertObject(mesh,study)

LX = max(mesh.X(:,2))-min(mesh.X(:,2));
LY = max(mesh.X(:,3))-min(mesh.X(:,3));

    if length(study.ObjectCoords)==4 % square
        Object_nodes = mesh.X( ...
            mesh.X(:,2)>=study.ObjectCoords(1) & mesh.X(:,2)<=study.ObjectCoords(2) & ...
            mesh.X(:,3)>=study.ObjectCoords(3) & mesh.X(:,3)<=study.ObjectCoords(4), ...
            1);
    elseif length(study.ObjectCoords)==3 % circle
        % Extract circle parameters
        centerX = study.ObjectCoords(1);
        centerY = study.ObjectCoords(2);
        radius = study.ObjectCoords(3);

        % Calculate the distance from each node to the center of the circle
        distances = sqrt((mesh.X(:,2) - centerX).^2 + (mesh.X(:,3) - centerY).^2);

        % Find nodes where the distance is less than or equal to the radius
        insideCircle = distances <= radius;
        Object_nodes = mesh.X(insideCircle, 1);

        if study.Cloak
            radius = 0.03757483113;

            % Calculate the distance from each node to the center of the circle
            distances = sqrt((mesh.X(:,2) - centerX).^2 + (mesh.X(:,3) - centerY).^2);
        
            % Find nodes where the distance is less than or equal to the radius
            insideCircle = distances <= radius;
            Object_nodes = mesh.X(insideCircle, 1);

            nodeCoords1 = [centerX, centerY+radius];
            nodeCoords2 = [centerX, centerY-radius];
            nodeCoords3 = [centerX+3*study.ObjectCoords(3), centerY];

            % Find nodes inside the triangle
            insideTriangle = inpolygon(mesh.X(:,2), mesh.X(:,3), [nodeCoords1(1), nodeCoords2(1), nodeCoords3(1), nodeCoords1(1)], [nodeCoords1(2), nodeCoords2(2), nodeCoords3(2), nodeCoords1(2)]);
            Object_nodes = [Object_nodes; mesh.X(insideTriangle, 1)];
        end

        % % Calculate penalty scale for each node within the circle
        % % Scaling factor decreases linearly from 1 at the center to 0 at the perimeter
        % scale = 1 - (distances(insideCircle) / radius);

        % % Scaling factor is alpha_upper from r=0 to r/2 and then decreases linearly to 0 at the perimeter
        % scale = 2 - (2*distances(insideCircle) /(radius));
        % scale(scale>1)=1;

        scale = ones(length(Object_nodes),1);

end
end