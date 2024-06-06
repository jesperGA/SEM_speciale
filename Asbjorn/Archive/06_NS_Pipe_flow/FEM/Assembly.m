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

    N = study.N;
    numElems = size(mesh.IX, 3);

    %% Global numbers (for safe keeping)
    opt.nel = numElems;
    opt.neqn_u = size(mesh.X, 1);
    opt.neqn_p = size(mesh.Xp, 1);
    opt.neqn = 2 * opt.neqn_u + opt.neqn_p;

    %% SUPPORTS (The N-Matrix)
    [opt.Null_u1, g_u1] = createSupportMatrix(mesh.bound_u1, opt.neqn_u);
    [opt.Null_u2, g_u2] = createSupportMatrix(mesh.bound_u2, opt.neqn_u);
    [opt.Null_p, g_p] = createSupportMatrix(mesh.bound_p, opt.neqn_p);
    opt.g = [g_u1; g_u2; g_p];
    
    %% LOAD VECTOR
    opt.f_u1 = study.F1(mesh.X(:,2),mesh.X(:,3));
    opt.f_u2 = study.F2(mesh.X(:,2),mesh.X(:,3));

    %% MATRICES 
    % Initialize matrix dimensions and get Gauss-Lobatto-Legendre points
    ldof_u = (N + 1) ^ 2;    
    ldof_p = (N-2 + 1) ^ 2;

    [xi, w, D] = GetGLL(N + 1);
    [xi_gl,w_gl]=lgwt(N-2+1,-1,1); 
    
    % Preallocate arrays
    xy = zeros(ldof_u, 2, numElems);
    edof_u = zeros(ldof_u, numElems);
    edof_p = zeros(ldof_p, numElems);
    A = zeros(ldof_u, ldof_u, numElems);
    B = zeros(ldof_u, ldof_u, numElems);
    D1 = zeros(ldof_p, ldof_u, numElems);
    D2 = zeros(ldof_p, ldof_u, numElems);
    M = zeros(ldof_p, ldof_p, numElems);
    
    % Compute element matrices in parallel
    for e = 1:numElems
        
        % Get element data
        [xy(:,:,e), edof_u(:,e)] = getElementData(mesh.X, mesh.IX, e);
        [~,    edof_p(:,e)] = getElementData(mesh.Xp, mesh.IXp, e);

        % Calculate stifness and mass matrices
        [A(:,:,e), B(:,:,e)] = elementMatrix2D_mex(xy(:,:,e), xi, w, D, N);

        % Calculate coupling matrices
        [D1(:,:,e), D2(:,:,e)] = calculateDMatrices(xy(:,:,e), N, xi, w, xi_gl, w_gl);

        % Precondition Matrix (Pressure mesh mass matrix)
        % [M(:,:,e)] = preCondMatr(w_gl, N-2, xy(:,:,e), ldof_p);

    end
    opt.Be = B;
    
    % Assemble global matrices
    [opt.A, opt.B] = assembleGlobalMatrices(edof_u, edof_u, A, B, ldof_u, ldof_u, numElems, opt.neqn_u, opt.neqn_u);
    [opt.D1, opt.D2] = assembleGlobalMatrices(edof_u, edof_p, D1, D2, ldof_u, ldof_p, numElems, opt.neqn_u, opt.neqn_p);
    % [opt.M, ~] = assembleGlobalMatrices(edof_p, edof_p, M, M, ldof_p, ldof_p, numElems, opt.neqn_p, opt.neqn_p);


    %% Assemble global system

    if strcmp(study.timeInt,'BDF1_AB3')
        gamma0 = 1;
    elseif strcmp(study.timeInt,'BDF3_EX3')
        gamma0 = 11/6;
    end

    dt = study.t(2) - study.t(1);
    H = opt.A + gamma0 * 1 / dt * opt.B;
    
    % Zero matrices
    Zer1 = sparse(opt.neqn_u,opt.neqn_u); 
    Zer2 = sparse(opt.neqn_p,opt.neqn_p);
    Zer3 = sparse(opt.neqn_u,opt.neqn_p);

    alpha = zeros(size(H));
    % Calculate the linear indices
    linearIndices = sub2ind(size(alpha), mesh.IX(:,:,5), mesh.IX(:,:,5));
    % Use linear indices to get the elements
    alpha(linearIndices)=1;

    opt.LHS = [H+alpha      Zer1       -opt.D1.';
               Zer1         H+alpha    -opt.D2.';
              -opt.D1      -opt.D2      Zer2];

    opt.Null = [opt.Null_u1    Zer1           Zer3
                Zer1           opt.Null_u2    Zer3
                Zer3'           Zer3'         opt.Null_p];

    % Modify for boundary conditions
    [opt.LHS_BC, ~] = applyBoundaryConditions(opt.LHS, zeros(opt.neqn,1), opt.Null, opt.g, opt.neqn);    

end

function [Null, g] = createSupportMatrix(bound, neqn)

    % Initialize vectors
    N = ones(neqn, 1);
    g = sparse(neqn, 1);
    if numel(bound)>0
        N(bound(:, 1)) = 0;
        g(bound(:, 1)) = bound(:, 3);
    end
    Null = spdiags(N, 0, neqn, neqn);
end


function [xy, edof] = getElementData(X, IX, e)
    nen = IX(:,:,e);
    xy = X(nen, 2:3);
    edof = reshape(nen, [], 1);
end

function [A_global, B_global] = assembleGlobalMatrices(edof1, edof2, ae, be, ldof1, ldof2, numElems, neqn1, neqn2)
    % Initialize arrays for triplets
    % Zero the arrays for the triplets
    I = zeros(numElems*ldof1*ldof2,1);
    J = zeros(numElems*ldof1*ldof2,1);
    A = zeros(numElems*ldof1*ldof2,1);
    B = zeros(numElems*ldof1*ldof2,1);

    % ASSEMBLY
    for e = 1:numElems

        % add to global system
        offset = (e - 1) * ldof2 * ldof1;
        [krow, kcol] = meshgrid(1:ldof2, 1:ldof1);
        indices = offset + (1:numel(krow))';
        linearIndices = sub2ind([ldof2, ldof1], krow, kcol);
        I(indices) = edof2(krow(:),e);
        J(indices) = edof1(kcol(:),e);

        a = ae(:,:,e);
        b = be(:,:,e);

        A(indices) = a(linearIndices);
        B(indices) = b(linearIndices);

    end

    ind = find(I>0);
    A_global = sparse(I(ind),J(ind),A(ind),neqn2,neqn1);
    B_global = sparse(I(ind),J(ind),B(ind),neqn2,neqn1);
end

