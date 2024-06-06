function opt = Assembly(mesh,study,opt)
% function opt = AssembleMindlin(mesh,study,opt)
% 
% Method that returns assemblies of K,M,C and P
% output: opt.K, opt.M, opt.C, opt.Null, opt.P

%% Global numbers (for safe keeping)
opt.nel = size(mesh.IX,1);
opt.neqn = size(mesh.X,1);

%% SUPPORTS (The N-Matrix)
Nmatr = ones(opt.neqn,1);
opt.g=zeros(opt.neqn,1);
for i=1:size(mesh.bound,1)
   Nmatr(mesh.bound(i,1)) = 0;
   opt.g(mesh.bound(i,1)) = mesh.bound(i,2);
end
opt.Null = spdiags(Nmatr,0,opt.neqn,opt.neqn);

%% LOAD VECTOR

opt.f = zeros(opt.neqn,1); 
opt.f(:)=study.f;

%% MATRICES (STIFFMESS, MASS and DAMPING)
ldof=mesh.N+1;

% Zero the arrays for the triplets
I = zeros(opt.nel*ldof*ldof,1);
J = zeros(opt.nel*ldof*ldof,1);
AE = zeros(opt.nel*ldof*ldof,1);
BE = zeros(opt.nel*ldof*ldof,1);

ntriplets = 0;

% Loop over element and integrate
for e=1:opt.nel

    % Get coordinates
    x = mesh.X(mesh.IX(e,2):mesh.IX(e,end-1),2);
    
    % Get element dofs
    edof=[mesh.IX(e,2):mesh.IX(e,study.N+2)];

    [xi,w]=lglnodes(study.N);
    % w(:)=sum(w)/(length(x));
    % xi(:) = linspace(xi(1),xi(end),length(xi));
    
    [A,B] = elementMatrixBar(x,xi,w);
    % [A,B] = elementMatrixBarTrapz(x,xi,w);
    % [A,B] = elementMatrixBarFEM(x);


    
    % add to global system (I,J,[KE,ME,CE])
    for krow = 1:ldof
        for kcol = 1:ldof
            ntriplets = ntriplets+1;
            I(ntriplets) = edof(krow);
            J(ntriplets) = edof(kcol);
            AE(ntriplets) = A(krow,kcol);
            BE(ntriplets) = B(krow,kcol);
        end
    end
    
end
ind = find(I>0);
opt.A = sparse(I(ind),J(ind),AE(ind),opt.neqn,opt.neqn);
opt.B = sparse(I(ind),J(ind),BE(ind),opt.neqn,opt.neqn);

opt.A_noBCs = opt.A;
opt.B_noBCs = opt.B;

% Modidy for BCs
opt.A = opt.Null'*opt.A*opt.Null - (opt.Null-speye(opt.neqn,opt.neqn));
opt.B = opt.Null'*opt.B*opt.Null;
opt.f = opt.Null*opt.f;
end