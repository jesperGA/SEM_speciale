function opt = Assembly(mesh, study, opt)
    % Assembly: Assembles global matrices for a finite element analysis
    % 
    % Input:
    % mesh - Struct containing mesh data
    % study - Struct containing study parameters
    % 
    % Output:
    % opt - Updated struct containing assembled matrices and vectors
    % opt.A, opt.B, opt.C, opt.Null, opt.P

    N = study.N;
    numElems = size(mesh.IX, 3);

    %% Global numbers (for safe keeping)
    opt.nel = numElems;
    opt.neqn_u = size(mesh.X, 1);
    opt.neqn_p = size(mesh.Xp, 1);
    opt.neqn = 2 * opt.neqn_u + opt.neqn_p;

    %% Assemble Null matrices for boundary conditions
    [opt.Null_u1, g_u1] = createSupportMatrix(mesh.bound_u1, opt.neqn_u);
    [opt.Null_u2, g_u2] = createSupportMatrix(mesh.bound_u2, opt.neqn_u);
    [opt.Null_p, g_p] = createSupportMatrix(mesh.bound_p, opt.neqn_p);
    opt.g = [g_u1; g_u2; g_p];
    
    %% LOAD VECTOR
    opt.f_u1 = sparse(study.F1(mesh.X(:,2), mesh.X(:,3), 0));
    opt.f_u2 = sparse(study.F2(mesh.X(:,2), mesh.X(:,3), 0));

    %% MATRICES 
    % Initialize matrix dimensions and get Gauss-Lobatto-Legendre points
    ldof_u = (N + 1) ^ 2;
    if strcmp(study.FEM,'LFEM') || strcmp(study.FEM,'QFEM')
        ldof_p = (N-1 + 1) ^ 2;
    else
        ldof_p = (N-2 + 1) ^ 2;
    end

    %% Prepare Gauss points, weights, and derivatives
    [xi, w, D] = GetGLL(study.N + 1);
    if strcmp(study.FEM,'LFEM') || strcmp(study.FEM,'QFEM')
        [xi_gl, w_gl] = lgwt(N, -1, 1);
    else
        [xi_gl, w_gl] = lgwt(N-1, -1, 1);
    end

    %% Preallocate arrays for element calculations
    % Preallocate arrays
    edof_u = zeros(ldof_u, numElems);
    edof_p = zeros(ldof_p, numElems);
    a  = zeros(ldof_u, ldof_u, numElems);
    b  = zeros(ldof_u, ldof_u, numElems);
    grad1 = zeros(ldof_u, ldof_u, numElems);
    grad2 = zeros(ldof_u, ldof_u, numElems);
    d1 = zeros(ldof_p, ldof_u, numElems);
    d2 = zeros(ldof_p, ldof_u, numElems);
    
    % Compute element matrices in parallel
    parfor e = 1:numElems
        
        % Get element data
        [xy, edof_u(:,e)] = getElementData(mesh.X, mesh.IX, e);
        [~ , edof_p(:,e)] = getElementData(mesh.Xp, mesh.IXp, e);

        % Calculate stifness and mass matrices
        [a(:,:,e), b(:,:,e), grad1(:,:,e), grad2(:,:,e), J(:,:,e)] = elementMatrix2D(xy, xi, w, D, N);

        % Calculate coupling matrices
        [d1(:,:,e), d2(:,:,e)] = calculateDMatrices(xy, N, xi, w, xi_gl, w_gl);
    end
    opt.b = b;
    opt.grad1 = grad1;
    opt.grad2 = grad2;
    opt.J = J;
    
    % Assemble global matrices
    [opt.A, opt.B] = assembleGlobalMatrices(edof_u, edof_u, a, b, ldof_u, ldof_u, numElems, opt.neqn_u, opt.neqn_u);
    [opt.Grad1, opt.Grad2] = assembleGlobalMatrices(edof_u, edof_u, grad1, grad2, ldof_u, ldof_u, numElems, opt.neqn_u, opt.neqn_u);
    [opt.D1, opt.D2] = assembleGlobalMatrices(edof_u, edof_p, d1, d2, ldof_u, ldof_p, numElems, opt.neqn_u, opt.neqn_p);

    %% Assemble global system

    dt = study.t(2) - study.t(1);
    mu = mesh.Material(1);
    rho = mesh.Material(2);
    if strcmp(study.timeInt,'BDF1_AB3')
        gamma0 = 1;
        H = mu * opt.A + gamma0 * rho / dt * opt.B;
    elseif strcmp(study.timeInt,'BDF3_EX3')
        gamma0 = 11/6;
        H = mu * opt.A + gamma0 * rho / dt * opt.B;
    elseif strcmp(study.timeInt,'FullyExplicit')
        H = rho / dt * opt.B;
    end

    
    % Zero matrices
    Zer1 = sparse(opt.neqn_u,opt.neqn_u); 
    Zer2 = sparse(opt.neqn_p,opt.neqn_p);
    Zer3 = sparse(opt.neqn_u,opt.neqn_p);

    % Brinkman term
    opt.Alpha = zeros(size(H));
    if strcmp(study.example,'Pipe')
        linearIndices = sub2ind(size(opt.Alpha), mesh.Object_nodes, mesh.Object_nodes);
        opt.Alpha(linearIndices)=study.alpha * mesh.scale; % Brinkman penalty
        opt.Alpha = sparse(opt.Alpha);
    end

    % System Left Hand Side assembly
    opt.LHS = [H+opt.Alpha      Zer1       -opt.D1.';
               Zer1         H+opt.Alpha    -opt.D2.';
              -opt.D1      -opt.D2      Zer2];
   
    % Boundary conditions integration into the system
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

