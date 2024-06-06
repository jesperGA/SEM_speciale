function [mesh] = boundaryConditions(study, mesh, NELX, NELY, U1, U2, P)
    % Applies boundary conditions to the mesh based on the study configuration.

    % Retrieve boundary nodes for setting conditions
    nodes = getBoundaryNodes(mesh, study);
    if strcmp(study.example, 'Pipe') && isfield(mesh, 'newBoundary')
        % Update boundary nodes for 'Pipe' example if new boundaries are defined
        nodes = union(nodes, mesh.newBoundary);
    end

    % Initialize boundary condition arrays for velocity components
    mesh.bound_u1 = zeros(length(nodes), 3);
    mesh.bound_u2 = zeros(length(nodes), 3);

    % Assign node indices
    mesh.bound_u1(:,1) = nodes;
    mesh.bound_u2(:,1) = nodes;

    % Apply velocity boundary conditions
    mesh.bound_u1(:,3) = U1(mesh.X(nodes, 2), mesh.X(nodes, 3), 0);
    mesh.bound_u2(:,3) = U2(mesh.X(nodes, 2), mesh.X(nodes, 3), 0);

    % Initialize pressure boundary conditions array
    mesh.bound_p = [];

end

function nodes = getBoundaryNodes(mesh, study)
    % Identifies boundary nodes based on the study configuration.

    % South boundary nodes
    nodes = mesh.X(mesh.X(:,3) == min(mesh.X(:,3)), 1);

    % East boundary nodes (exclude for 'Pipe')
    if ~strcmp(study.example, 'Pipe')
        eastNodes = mesh.X(mesh.X(:,2) == max(mesh.X(:,2)), 1);
        nodes = [nodes; eastNodes(2:end)];  % Skip duplicated corner node
    end

    % North boundary nodes (reversed order)
    northNodes = mesh.X(mesh.X(:,3) == max(mesh.X(:,3)), 1);
    nodes = [nodes; northNodes(end:-1:1)];  % Add in reversed order
    
    % West boundary nodes (reversed order, skip corners)
    westNodes = mesh.X(mesh.X(:,2) == min(mesh.X(:,2)), 1);
    nodes = [nodes; westNodes(end-1:-1:2)];  % Skip both corners
end
