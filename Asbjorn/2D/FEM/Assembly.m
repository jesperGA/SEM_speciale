function opt = Assembly(mesh, study, opt)
    % Assembly: Assembles global matrices for a finite element analysis
    % 
    % Input:
    % mesh - Struct containing mesh data
    % study - Struct containing study parameters
    % opt - Struct for storing assembled matrices and vectors
    % 
    % Output:
    % opt - Updated struct containing assembled matrices and vectors
    % opt.K, opt.M, opt.C, opt.Null, opt.P

    %% Global numbers (for safe keeping)
    opt.nel = size(mesh.IX, 3);
    opt.neqn = size(mesh.X, 1);

    %% SUPPORTS (The N-Matrix)
    [N, opt.g] = createSupportMatrix(mesh, opt);
    opt.Null = spdiags(N,0,opt.neqn,opt.neqn);
    
    %% LOAD VECTOR
    opt.f = zeros(opt.neqn, 1); 

    %% MATRICES (STIFFNESS, MASS and DAMPING)
    ldof = (study.N + 1) ^ 2;

    % Initialize arrays for triplets
    [I, J, AE, BE, ntriplets] = initializeTriplets(opt, ldof);

    % Loop over elements and integrate
    for e = 1:opt.nel
        % Get element data
        [xy, edof] = getElementData(mesh, study, e);

        % Calculate element matrices
        [A, B] = calculateElementMatrices(xy, study, edof);

        % Add to global system
        [I, J, AE, BE, ntriplets] = addToGlobalSystem(A, B, edof, I, J, AE, BE, ntriplets, ldof);
    end

    % Assemble global matrices
    opt = assembleGlobalMatrices(I, J, AE, BE, opt, ntriplets);

    % Modify for boundary conditions
    opt = applyBoundaryConditions(opt);
end

function [N, g] = createSupportMatrix(mesh, opt)
    N = ones(opt.neqn, 1);
    g = sparse(opt.neqn, 1);
    for i = 1:size(mesh.bound, 1)
        N(mesh.bound(i, 1)) = 0;
        g(mesh.bound(i, 1)) = mesh.bound(i, 3);
    end
    opt.Null = spdiags(N, 0, opt.neqn, opt.neqn);
end

function [I, J, AE, BE, ntriplets] = initializeTriplets(opt, ldof)
    I = zeros(opt.nel * ldof * ldof, 1);
    J = zeros(opt.nel * ldof * ldof, 1);
    AE = zeros(opt.nel * ldof * ldof, 1);
    BE = zeros(opt.nel * ldof * ldof, 1);
    ntriplets = 0;
end

function [xy, edof] = getElementData(mesh, study, e)
    nen = mesh.IX(:,:,e);
    xy = mesh.X(nen, 2:3);
    edof = reshape(nen, [], 1);
end

function [A, B] = calculateElementMatrices(xy, study, edof)
    [xi, w, h] = GetGLL(study.N + 1);
    x = reshape(xy(:, 1), study.N + 1, study.N + 1);
    y = reshape(xy(:, 2), study.N + 1, study.N + 1);
    [A, B] = elementMatrix2D(x, y, xi, w, h', study.N);
end

function [I, J, AE, BE, ntriplets] = addToGlobalSystem(A, B, edof, I, J, AE, BE, ntriplets, ldof)
    for krow = 1:ldof
        for kcol = 1:ldof
            ntriplets = ntriplets + 1;
            I(ntriplets) = edof(krow);
            J(ntriplets) = edof(kcol);
            AE(ntriplets) = A(krow, kcol);
            BE(ntriplets) = B(krow, kcol);
        end
    end
end

function opt = assembleGlobalMatrices(I, J, AE, BE, opt, ntriplets)
    ind = find(I > 0);
    opt.A = sparse(I(ind), J(ind), AE(ind), opt.neqn, opt.neqn);
    opt.B = sparse(I(ind), J(ind), BE(ind), opt.neqn, opt.neqn);
end

function opt = applyBoundaryConditions(opt)
    % applyBoundaryConditions: Applies boundary conditions to the system matrices
    % 
    % Input:
    % opt - Struct containing system matrices and boundary conditions
    % 
    % Output:
    % opt - Updated struct with modified system matrices and load vector

    % Identify free degrees of freedom
    free_dofs = diag(opt.Null);

    % Create a diagonal sparse matrix of boundary values
    G = spdiags(opt.g, 0, opt.neqn, opt.neqn);

    % Modify the load vector
    opt.f = opt.B * opt.f;                   % Apply mass matrix to load vector
    opt.f = opt.f - sum(opt.A * G, 2);       % Adjust for stiffness matrix
    opt.f(~free_dofs) = opt.g(~free_dofs);   % Enforce boundary values on the load vector

    % Modify the stiffness matrix for boundary conditions
    opt.A = opt.Null' * opt.A * opt.Null - (opt.Null - speye(opt.neqn, opt.neqn));
end
