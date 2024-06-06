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
    opt.neqn_u = size(mesh.X, 1);
    opt.neqn_p = size(mesh.Xp, 1);
    opt.neqn = 2 * opt.neqn_u + opt.neqn_p;


    %% SUPPORTS (The N-Matrix)
    [N_u1, g_u1] = createSupportMatrix(mesh.bound_u1, opt.neqn_u);
    [N_u2, g_u2] = createSupportMatrix(mesh.bound_u2, opt.neqn_u);
    [N_p, g_p] = createSupportMatrix(mesh.bound_p, opt.neqn_p);
    opt.g = [g_u1; g_u2; g_p];
    opt.Null_u1 = spdiags(N_u1,0,opt.neqn_u,opt.neqn_u);
    opt.Null_u2 = spdiags(N_u2,0,opt.neqn_u,opt.neqn_u);
    opt.Null_p = spdiags(N_p,0,opt.neqn_p,opt.neqn_p);
    
    %% LOAD VECTOR
    opt.f_u1 = zeros(opt.neqn_u, 1); 
    opt.f_u2 = zeros(opt.neqn_u, 1); 

    opt.f_u1 = study.F1(mesh.X(:,2),mesh.X(:,3));
    opt.f_u2 = study.F2(mesh.X(:,2),mesh.X(:,3));

    %% MATRICES (STIFFNESS and MASS)
    ldof = (study.N + 1) ^ 2;

    % Initialize arrays for triplets
    [I, J, AE, BE, ntriplets] = initializeTriplets(opt, ldof);
    opt.B = zeros(ldof,ldof,opt.nel);

    % Loop over elements and integrate
    for e = 1:opt.nel
        % Get element data
        [xy, edof] = getElementData(mesh.X, mesh.IX, e);

        % Calculate element matrices
        [A, opt.Be(:,:,e)] = calculateElementMatrices(xy, study.N);

        % Add to global system
        [I, J, AE, BE, ntriplets] = addToGlobalSystem(A, opt.Be(:,:,e), edof, I, J, AE, BE, ntriplets, ldof);


                % Precondition Matrix (Pressure mesh mass matrix)
        % [M(:,:,e)] = preCondMatr(w_gl, N-2, xy(:,:,e), ldof_p);
    end
    % [opt.M, ~] = assembleGlobalMatrices(edof_p, edof_p, M, M, ldof_p, ldof_p, numElems, opt.neqn_p, opt.neqn_p);

    % Assemble global matrices
    opt = assembleGlobalMatrices(I, J, AE, BE, opt, ntriplets);

    %% Coupling MATRICES (D1,D2)
    ldof_u = (study.N + 1) ^ 2;
    ldof_p = (study.N-2 + 1) ^ 2;

    % Initialize arrays for triplets
    % Zero the arrays for the triplets
    I = zeros(opt.nel*ldof_u*ldof_p,1);
    J = zeros(opt.nel*ldof_u*ldof_p,1);
    DE1 = zeros(opt.nel*ldof_u*ldof_p,1);
    DE2 = zeros(opt.nel*ldof_u*ldof_p,1);
    ntriplets = 0;

    [xi, w] = lglnodes(study.N); % Legendre-Gauss-Lobatto nodes and weights
    [xi_gl,w_gl]=lgwt(study.N-2+1,-1,1); 

    % Loop over elements and integrate
    for e = 1:opt.nel
        % Get element data
        [xy_u, edof_u] = getElementData(mesh.X,  mesh.IX,  e);
        [~,    edof_p] = getElementData(mesh.Xp, mesh.IXp, e);

        % Calculate coupling matrices
        [D1, D2] = calculateDMatrices(xy_u, study.N, xi, w, xi_gl, w_gl);

        % add to global system
        for krow = 1:ldof_p
            for kcol = 1:ldof_u
                ntriplets = ntriplets+1;
                I(ntriplets) = edof_p(krow);
                J(ntriplets) = edof_u(kcol);
                DE1(ntriplets) = D1(krow,kcol);
                DE2(ntriplets) = D2(krow,kcol);
            end
        end
    end
    ind = find(I>0);
    opt.D1 = sparse(I(ind),J(ind),DE1(ind),opt.neqn_p,opt.neqn_u);
    opt.D2 = sparse(I(ind),J(ind),DE2(ind),opt.neqn_p,opt.neqn_u);


    %% Precondition Matrix (Pressure mesh mass matrix)
    ldof = (study.N - 2 + 1) ^ 2;

    % Initialize arrays for triplets
    [I, J, ME, ~, ntriplets] = initializeTriplets(opt, ldof);
    
    [xy, ~] = getElementData(mesh.X, mesh.IX, 1);
    L1 = max(xy(:,1)) - min(xy(:,1));
    L2 = max(xy(:,2)) - min(xy(:,2));

    % Loop over elements and integrate
    for e = 1:opt.nel

        % Get element data
        [xy, edof] = getElementData(mesh.Xp, mesh.IXp, e);

        [~, w]=lgwt(study.N + 1 - 2,-1,1); 

        [M] = preCondMatr(w, study.N-2, L1, L2);

        % add to global system (I,J,[SE])
        for krow = 1:ldof_p
            for kcol = 1:ldof_p
                ntriplets = ntriplets+1;
                I(ntriplets) = edof(krow);
                J(ntriplets) = edof(kcol);
                ME(ntriplets) = M(krow,kcol);
            end
        end
    end

    ind = find(I>0);
    opt.M = sparse(I(ind),J(ind),ME(ind),opt.neqn_p,opt.neqn_p);

    %% Assemble global system
    rho = mesh.Material(1);
    mu  = mesh.Material(2);

    if strcmp(study.timeInt,'BDF1_AB3')
        gamma0 = 1;
    elseif strcmp(study.timeInt,'BDF3_EX3')
        gamma0 = 11/6;
    end

    dt = study.t(2) - study.t(1);
    H = mu * opt.A + gamma0 * rho / dt * opt.B;
    
    % Zero matrices
    Zer1 = sparse(opt.neqn_u,opt.neqn_u); 
    Zer2 = sparse(opt.neqn_p,opt.neqn_p);
    Zer3 = sparse(opt.neqn_u,opt.neqn_p);

    opt.LHS = [H            Zer1        -opt.D1.';
           Zer1         H           -opt.D2.';
           -opt.D1      -opt.D2      Zer2];

    opt.Null = [opt.Null_u1    Zer1           Zer3
                Zer1           opt.Null_u2    Zer3
                Zer3'           Zer3'         opt.Null_p];

    % Modify for boundary conditions
    [opt.LHS_BC, ~] = applyBoundaryConditions(opt.LHS, zeros(opt.neqn,1), opt.Null, opt.g, opt.neqn);    

end

function [N, g] = createSupportMatrix(bound, neqn)
    N = ones(neqn, 1);
    g = sparse(neqn, 1);
    for i = 1:size(bound, 1)
        N(bound(i, 1)) = 0;
        g(bound(i, 1)) = bound(i, 3);
    end
end

function [I, J, AE, BE, ntriplets] = initializeTriplets(opt, ldof)
    I = zeros(opt.nel * ldof * ldof, 1);
    J = zeros(opt.nel * ldof * ldof, 1);
    AE = zeros(opt.nel * ldof * ldof, 1);
    BE = zeros(opt.nel * ldof * ldof, 1);
    ntriplets = 0;
end

function [xy, edof] = getElementData(X, IX, e)
    nen = IX(:,:,e);
    xy = X(nen, 2:3);
    edof = reshape(nen, [], 1);
end

function [A, B] = calculateElementMatrices(xy, N)
    
    [xi, w, h] = GetGLL(N + 1);

    x = reshape(xy(:, 1), N + 1, N + 1);
    y = reshape(xy(:, 2), N + 1, N + 1);
    
    [A, B] = elementMatrix2D(x, y, xi, w, h, N);
    
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
    opt.A = sparse(I(ind), J(ind), AE(ind), opt.neqn_u, opt.neqn_u);
    opt.B = sparse(I(ind), J(ind), BE(ind), opt.neqn_u, opt.neqn_u);
end
