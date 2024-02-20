function [mesh] = mesh2D(LX, LY, NELX, NELY, N)

    dxe = LX / NELX;
    dye = dxe;
    NEL = NELX * NELY;
    NGLL = N + 1; % Number of GLL nodes per element
    
    % Generate the Spectral Element mesh
    [mesh.IX, x, y] = MeshBox(LX, LY, NELX, NELY, NGLL, 0);
    
    %-------------------------------------------------------------------------%
    %                             Node coordinates and numbers                %
    %-------------------------------------------------------------------------%
    node_coordinates = [x, y];

    % Adjust y-coordinates
    element_nodes = [mesh.IX(:, :, 3), mesh.IX(:, :, 4)];
    xx = node_coordinates(element_nodes, 1);
    yy = node_coordinates(element_nodes, 2);
    ratios = 1 + 1/4 * sinpi(xx(yy == 1)) ./ (1/2);
    % ratios(end/2+1:end)= ratios(end/2:-1:1);
    for i = NGLL:-1:1
        yy(yy == yy(i * N + 1)) = 1/2 + (yy(yy == yy(i * N + 1)) - 1/2) .* ratios;    
    end

    node_coordinates(element_nodes, 2) = yy;

    mesh.X = [(1:length(x))', node_coordinates(:, 1), node_coordinates(:, 2)];
    
    %-------------------------------------------------------------------------%
    %                             Material parameters                         %
    %-------------------------------------------------------------------------%
    
    mesh.Material = [1.1, 1.2, 1.3, 1.4];

end