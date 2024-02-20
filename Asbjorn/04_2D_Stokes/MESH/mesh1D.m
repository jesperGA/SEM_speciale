function mesh = mesh1D(L, ne, N)
    % meshBraggGrating: Generates a mesh for Bragg grating simulation.
    % 
    % Inputs:
    %   L: Length of the Bragg grating.
    %   ne: Number of elements (subdivisions) along the length of the Bragg grating.
    %   N: Polynomial order.
    %
    % Outputs:
    %   mesh: Structure containing the mesh information.
    %     mesh.N: Polynomial order.
    %     mesh.X: Node coordinates matrix (nn x 2).
    %     mesh.IX: Element connectivity matrix (ne x (N + 3)).
    %     mesh.mu: Design variable vector (ne x 1).
    %     mesh.Material: Matrix (2 x 4) defining material properties for two different materials.
    %     mesh.inactive_idx: Indices of inactive elements in the mesh.
    %     mesh.active_idx: Indices of active elements in the mesh.

    % Number of nodes
    nn = ne * N + 1;

    % Containers
    mesh.N = N;
    mesh.X = zeros(nn, 2);
    mesh.IX = zeros(ne, N + 3);
    mesh.mu = zeros(ne, 1);

    % Define material properties
    material1 = [1, 1, 1, 1];
    material2 = [1, 1.5, 1, 1.3];
    mesh.Material = [material1; material2];

    % Fill Node and Element containers
    xi = lglnodes(N);
    mesh.X(1, 1:2) = [1, 0];
    dx_e = L / ne;
    for e = 1:ne
        index_nodes = (e - 1) * N + 2 : e * N + 1;
        x_old = mesh.X(index_nodes(1) - 1, 2);
        mesh.X(index_nodes, 1) = [(e - 1) * N + 2 : e * N + 1];
        mesh.X(index_nodes, 2) = x_old + (xi(2:end) + 1) / 2 * dx_e;
        edof = [(e - 1) * N + 1, e * N + 1];
        mesh.IX(e, 1) = e;
        mesh.IX(e, 2:N + 2) = edof(1):edof(2);
    end

    % Material parameter 1
    mesh.IX(:,end) = 1;

    % bc
    mesh.bound=[1 0;
                nn 0];
end
