function [mesh] = meshStaggered(study, LX, LY, NELX, NELY, N)
    % meshStaggered Generates staggered meshes for different simulations based on the study settings.
    % Inputs:
    %   study - Struct containing the study parameters
    %   LX - Length of the domain in the x-direction
    %   LY - Length of the domain in the y-direction
    %   NELX - Number of elements in the x-direction
    %   NELY - Number of elements in the y-direction
    %   N - Polynomial degree inside each element
    % Outputs:
    %   mesh - Struct containing the mesh details

    % Initialize standard mesh and a refined mesh for pressure calculation
    mesh = mesh2D(LX, LY, NELX, NELY, N, 0);
    if strcmp(study.FEM, 'LFEM') || strcmp(study.FEM, 'QFEM')
        meshP = mesh2D(LX, LY, NELX, NELY, N-1, 2);
    else
        meshP = mesh2D(LX, LY, NELX, NELY, N-2, 2);
    end

    % Adjust mesh coordinates for specific examples
    if strcmp(study.example, 'TaylorVortex') || strcmp(study.example, 'LidDriven')
        adjustment = LX / 2;
        mesh.X(:, 2:3) = mesh.X(:, 2:3) - adjustment;
        meshP.X(:, 2:3) = meshP.X(:, 2:3) - adjustment;
    end

    % Assign mesh element/node indices
    mesh.Xp = meshP.X;
    mesh.IXp = meshP.IX;

    % Material parameters and element dimensions
    switch study.example
        case 'TaylorVortex'
            mesh.Material = [1, 1]; % mu, rho
        case 'Duct'
            mesh.Material = [1e-3, 1]; % mu, rho
        case 'LidDriven'
            mesh.Material = [1e-2, 1]; % mu, rho
    end

    % Calculate element dimensions
    for e = 1:size(mesh.IX, 3)
        nen = mesh.IX(:, :, e);
        xy = mesh.X(nen, 2:3);
        mesh.L1(e) = max(xy(:, 1)) - min(xy(:, 1));
        mesh.L2(e) = max(xy(:, 2)) - min(xy(:, 2));
    end
end
