function [mesh] = meshStaggered(study, LX, LY, NELX, NELY, N)
    % Generates staggered meshes for different simulations based on the study settings.

    % Initialize standard mesh and a refined mesh for pressure calculation
    [mesh] = mesh2D(LX, LY, NELX, NELY, N, 0);
    if strcmp(study.FEM,'LFEM') || strcmp(study.FEM,'QFEM')
        [meshP] = mesh2D(LX, LY, NELX, NELY, N-1, 2);
    else 
        [meshP] = mesh2D(LX, LY, NELX, NELY, N-2, 2);
    end


    % Adjust mesh coordinates for specific examples
    if strcmp(study.example,'Roenquist_NS') || strcmp(study.example,'LidDriven')
        adjustment = LX / 2;
        mesh.X(:,2:3) = mesh.X(:,2:3) - adjustment;
        meshP.X(:,2:3) = meshP.X(:,2:3) - adjustment;
    end

    % Assign mesh element/node indices
    mesh.Xp = meshP.X;
    mesh.IXp = meshP.IX;

    % % Remove elements for specific example
    % if strcmp(study.example, 'Pipe')
    %     removeEleNr = round(size(mesh.IX, 3) / 4);
    %     [mesh.X, mesh.IX, mesh.newBoundary] = remove_element(mesh.X, mesh.IX, removeEleNr);
    %     [mesh.Xp, mesh.IXp] = remove_element(mesh.Xp, mesh.IXp, removeEleNr);
    % end

    % Material parameters and element dimensions
    if strcmp(study.example,'Roenquist_NS') || strcmp(study.example,'LidDriven')
        mesh.Material = [1, 1]; %mu, rho
    elseif strcmp(study.example,'Pipe')
        mesh.Material = [1e-3, 1]; %mu, rho
    end

    for e = 1:size(mesh.IX, 3)
        nen = mesh.IX(:,:,e);
        xy = mesh.X(nen, 2:3);
        mesh.L1(e) = max(xy(:,1)) - min(xy(:,1));
        mesh.L2(e) = max(xy(:,2)) - min(xy(:,2));
    end
end

function [X, IX, newBoundary] = remove_element(X, IX, n_remove)
    % Remove specified element from the mesh and adjust connectivity.
    remove_dofs = IX(:,:,n_remove);
    IX = IX(:, :, [1:n_remove-1, n_remove+1:end]);
    newBoundary = intersect(remove_dofs(:), IX(:));
    remove_dofs = setdiff(remove_dofs(:), IX(:));

    rows_to_keep = true(size(X, 1), 1);
    rows_to_keep(remove_dofs) = false;
    X = X(rows_to_keep, :);

    % Adjust indexing for removed elements
    mapping_from = X(:,1);
    mapping_to = (1:length(X))';
    IX = replace_numbers(IX, mapping_from, mapping_to);
    X(:,1) = (1:length(X))';

    [~, idx] = ismember(newBoundary, mapping_from);
    newBoundary = mapping_to(idx);    
end

function IX_new = replace_numbers(IX, mapping_from, mapping_to)
    % Replace old indices with new indices in the index matrix
    IX_1d = IX(:);
    [~, idx] = ismember(IX_1d, mapping_from);
    IX_new = mapping_to(idx);
    IX_new = reshape(IX_new, size(IX));
end
