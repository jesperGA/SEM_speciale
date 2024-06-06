function [mesh] = mesh2D(LX, LY, NELX, NELY, N, distribution)
    % Creates a 2D spectral element mesh for a given domain and distribution type.
    
    % Calculate the element dimensions
    dxe = LX / NELX;
    dye = dxe; % Assuming square elements for simplicity, adjust if necessary for non-square

    % Compute the number of GLL (Gauss-Lobatto-Legendre) nodes per element
    Npoints = N + 1;

    % Generate the spectral element mesh coordinates and connectivity
    [mesh.IX, x, y] = MeshBox(LX, LY, NELX, NELY, Npoints, distribution);

    % Combine node numbers with their coordinates
    node_numbers = (1:length(x))';
    node_coordinates = [x, y];

    mesh.X = [node_numbers, node_coordinates(:, 1), node_coordinates(:, 2)];
end
