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

%% MATRICES (STIFFMESS, MASS and DAMPING)
ldof=(study.N+1)^2;

% Zero the arrays for the triplets
I = zeros(opt.nel*ldof*ldof,1);
J = zeros(opt.nel*ldof*ldof,1);
KE = zeros(opt.nel*ldof*ldof,1);
ME = zeros(opt.nel*ldof*ldof,1);
CE = zeros(opt.nel*ldof*ldof,1);
ntriplets = 0;

% Loop over element and integrate
for e=1:opt.nel

    % element nodes
    nen = mesh.IX(:,:,e);
    
    % Get coordinates
    xy = mesh.X(nen,2:3);
    
    % Get element dofs
    for i=1:(study.N+1)^2
        edof(2*i-1) = 2*nen(i)-1;
        edof(2*i-0) = 2*nen(i)-0;
    end

    % Get material parameters
    % matID = mesh.IX(e,6);
    matID = 1;
    thk = mesh.Material(matID,1);        E   = mesh.Material(matID,2);
    nu  = mesh.Material(matID,3);        rho = mesh.Material(matID,4);  

    [xi,w,h] = GetGLL(study.N+1);
    % w(:)=sum(w)/(length(x));
    % xi(:) = linspace(xi(1),xi(end),length(xi));
    
    % x=xy(:,1);
    % y=xy(:,2);
    x = reshape(xy(:,1),study.N+1,study.N+1);
    y = reshape(xy(:,2),study.N+1,study.N+1);
    [A,B] = elementMatrix2D(x,y,xi,w,h',study.N);
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

% Modidy for BCs
opt.A = opt.Null'*opt.A*opt.Null - (opt.Null-speye(opt.neqn,opt.neqn));
opt.B = opt.Null'*opt.B*opt.Null;
opt.f = opt.Null*opt.f;
end