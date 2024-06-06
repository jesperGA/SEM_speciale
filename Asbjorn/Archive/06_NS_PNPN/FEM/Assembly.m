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
    ldof = (N + 1) ^ 2;

    % Initialize arrays for triplets
    [I, J, AE, BE, ntriplets] = initializeTriplets(opt, ldof);
    grad1E = zeros(opt.nel * ldof * ldof, 1);
    grad2E = zeros(opt.nel * ldof * ldof, 1);
    gradgrad1E = zeros(opt.nel * ldof * ldof, 1);
    gradgrad2E = zeros(opt.nel * ldof * ldof, 1);
    opt.B = zeros(ldof,ldof,opt.nel);

    [xi, w, h] = GetGLL(N + 1);
    for j = 1:N+1
        for i = 1:length(xi)
            DD(i,j) = doublederivativeLagrange(j, xi', xi(i));
        end
    end

    % Loop over elements and integrate
    for e = 1:opt.nel
        % Get element data
        [xy, edof] = getElementData(mesh.X, mesh.IX, e);

        x = reshape(xy(:, 1), N + 1, N + 1);
        y = reshape(xy(:, 2), N + 1, N + 1);
        
        [A, opt.Be(:,:,e), grad1, grad2, gradgrad1, gradgrad2] = elementMatrix2D(x, y, xi, w, h, DD, N);

        for krow = 1:ldof
            for kcol = 1:ldof
                ntriplets = ntriplets + 1;
                I(ntriplets) = edof(krow);
                J(ntriplets) = edof(kcol);
                AE(ntriplets) = A(krow, kcol);
                BE(ntriplets) = opt.Be(krow, kcol,e);
                grad1E(ntriplets) = grad1(krow, kcol);
                grad2E(ntriplets) = grad2(krow, kcol);
                gradgrad1E(ntriplets) = gradgrad1(krow, kcol);
                gradgrad2E(ntriplets) = gradgrad2(krow, kcol);
            end
        end
    end

    ind = find(I > 0);
    opt.A = sparse(I(ind), J(ind), AE(ind), opt.neqn_u, opt.neqn_u);
    opt.B = sparse(I(ind), J(ind), BE(ind), opt.neqn_u, opt.neqn_u);
    opt.grad1 = sparse(I(ind), J(ind), grad1E(ind), opt.neqn_u, opt.neqn_u);
    opt.grad2 = sparse(I(ind), J(ind), grad2E(ind), opt.neqn_u, opt.neqn_u);
    opt.gradgrad1 = sparse(I(ind), J(ind), gradgrad1E(ind), opt.neqn_u, opt.neqn_u);
    opt.gradgrad2 = sparse(I(ind), J(ind), gradgrad2E(ind), opt.neqn_u, opt.neqn_u);

    %% Coupling MATRICES (D1,D2)
    ldof_u = (N + 1) ^ 2;
    ldof_p = (N + 1) ^ 2;

    % Initialize arrays for triplets
    % Zero the arrays for the triplets
    I = zeros(opt.nel*ldof_u*ldof_p,1);
    J = zeros(opt.nel*ldof_u*ldof_p,1);
    DE1 = zeros(opt.nel*ldof_u*ldof_p,1);
    DE2 = zeros(opt.nel*ldof_u*ldof_p,1);
    ntriplets = 0;

    [xi, w] = lglnodes(N); % Legendre-Gauss-Lobatto nodes and weights

    % Loop over elements and integrate
    for e = 1:opt.nel
        % Get element data
        [xy_u, edof_u] = getElementData(mesh.X,  mesh.IX,  e);
        [~,    edof_p] = getElementData(mesh.Xp, mesh.IXp, e);

        % Calculate coupling matrices
        [D1, D2] = calculateDMatrices(xy_u, N, xi, w, xi, w);

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
    ldof = (N + 1) ^ 2;

    % Initialize arrays for triplets
    [I, J, ME, ~, ntriplets] = initializeTriplets(opt, ldof);
    
    [xy, ~] = getElementData(mesh.X, mesh.IX, 1);
    L1 = max(xy(:,1)) - min(xy(:,1));
    L2 = max(xy(:,2)) - min(xy(:,2));

    % Loop over elements and integrate
    for e = 1:opt.nel

        % Get element data
        [xy, edof] = getElementData(mesh.Xp, mesh.IXp, e);

        [~, w]=GetGLL(N + 1); 
        [~, w]=lgwt(N + 1 ,-1,1); 

        [M] = preCondMatr(w, N, L1, L2);

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

    dt = study.t(2) - study.t(1);

    H = mu * opt.A + rho / dt * opt.B;

    
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

function ddl = doublederivativeLagrange(j, xi, x)
% Calculates the second derivative of the j-th Lagrange interpolation polynomial at a point x
%
% Inputs:
%   j - Index of the polynomial
%   xi - Vector of x coordinates
%   x - The point at which to evaluate the derivative
%
% Output:
%   ddl - The value of the second derivative at point x

k = length(xi); % Number of points

% Initialize the second derivative sum
ddl = 0;
for i = 1:k
    if i ~= j
        sum_m = 0;
        term1 = 1 / (xi(j) - xi(i));
        for m = 1:k
            if m ~= j && m ~= i
            term = 1 / (xi(j) - xi(m));
                for n = 1:k
                    if n ~= j && n ~= i && n ~= m
                        term = term * (x - xi(n)) / (xi(j) - xi(n));
                    end
                end
                sum_m = sum_m + term; 
            end
        end
        ddl = ddl + term1 * sum_m;
    end
end
end