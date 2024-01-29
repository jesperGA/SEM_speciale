function [mesh] = mesh2D(LX, LY, NELX, NELY, N)

    dxe = LX / NELX;
    dye = dxe;
    NEL = NELX * NELY;
    NGLL = N + 1; % Number of GLL nodes per element
    
    % Generate the Spectral Element mesh
    [mesh.IX, x, y] = MeshBox(LX, LY, NELX, NELY, NGLL);
    
    %-------------------------------------------------------------------------%
    %                             Node coordinates and numbers                %
    %-------------------------------------------------------------------------%
    node_coordinates = [x, y];
    node_numbers = 1:size(node_coordinates, 1);
    
    % Adjust y-coordinates
    element_nodes = [mesh.IX(:, :, 3), mesh.IX(:, :, 4)];
    xx = node_coordinates(element_nodes, 1);
    yy = node_coordinates(element_nodes, 2);
    ratios = 1 + 1/4 * sinpi(xx(yy == 1)) ./ (1/2);
    for i = NGLL:-1:1
        yy(yy == yy(i * N + 1)) = 1/2 + (yy(yy == yy(i * N + 1)) - 1/2) .* ratios;    
    end
    
    node_coordinates(element_nodes, 2) = yy;
    
    mesh.X = [(1:length(x))', node_coordinates(:, 1), node_coordinates(:, 2)];
    
    %-------------------------------------------------------------------------%
    %                             Material parameters                         %
    %-------------------------------------------------------------------------%
    
    mesh.Material = [1.1, 1.2, 1.3, 1.4];
    
    %-------------------------------------------------------------------------%
    %                             Boundary Conditions                         %
    %-------------------------------------------------------------------------%
    % South, East, North, West BCs
    idx = 0;
    F = @(X, Y) sin(X) .* exp(-Y);
    mesh.bound = zeros(4*(NELX*N+1)-4,3);
    % South BCs
    idx=NELX*N+1;
    nodes = mesh.X(mesh.X(:,3)==0,1);
    mesh.bound(1:idx,1) =nodes;
    mesh.bound(1:idx,3) =  F(mesh.X(nodes,2),mesh.X(nodes,3));
    % East BCs
    idx=idx + NELY*N;
    nodes = mesh.X(mesh.X(:,2)==max(mesh.X(:,2)),1);
    mesh.bound(idx-NELX*N:idx,1) =nodes;
    mesh.bound(idx-NELX*N:idx,3) = F(mesh.X(nodes,2),mesh.X(nodes,3));
    % North BCs
    idx=idx + NELY*N;
    nodes = [mesh.IX(:,end,3);mesh.IX(2:end,end,4)]; nodes=nodes(end:-1:1);
    mesh.bound(idx-NELX*N:idx,1) =nodes;
    mesh.bound(idx-NELX*N:idx,3) = F(mesh.X(nodes,2),mesh.X(nodes,3));
    % West BCs
    idx=idx + NELY*N;
    nodes = mesh.X(mesh.X(:,2)==0,1);nodes=nodes(end:-1:2);
    mesh.bound(idx-NELX*N:idx-1,1) =nodes;
    mesh.bound(idx-NELX*N:idx-1,3) = F(mesh.X(nodes,2),mesh.X(nodes,3));

end