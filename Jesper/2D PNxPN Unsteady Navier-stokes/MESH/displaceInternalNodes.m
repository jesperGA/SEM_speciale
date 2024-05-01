function Mesh = displaceInternalNodes(Mesh, dx, dy)
    % Get the number of nodes horizontally and vertically
    numXNodes = length(unique(Mesh.Xv(:,2)));
    numYNodes = length(unique(Mesh.Xv(:,3)));
    
    % Find the center node coordinates
    centerX = median(unique(Mesh.Xv(:,2)));
    centerY = median(unique(Mesh.Xv(:,3)));

    % Displacement should be maximum at the center and reduce to the edges
    for i = 1:size(Mesh.Xv, 1)
        x = Mesh.Xv(i, 2);
        y = Mesh.Xv(i, 3);

        % Calculate distance ratio from the center for x and y
        xRatio = 1 - abs(x - centerX) / (max(unique(Mesh.Xv(:,2))) - centerX);
        yRatio = 1 - abs(y - centerY) / (max(unique(Mesh.Xv(:,3))) - centerY);

        % Apply scaled displacement based on distance from center
        Mesh.Xv(i, 2) = x + dx * xRatio;
        Mesh.Xv(i, 3) = y + dy * yRatio;
    end
end
