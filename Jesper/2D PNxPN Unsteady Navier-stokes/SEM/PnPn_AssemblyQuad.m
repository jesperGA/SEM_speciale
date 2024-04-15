function [opt,study] = PnPn_AssemblyQuad(mesh,opt,study)

n_GLL = study.n_GLL; %Specify number of GLL points
N = n_GLL-1;
ldof =(n_GLL^2); %One dof per gll in every element

opt.nel = size(mesh.IXp,3);
opt.neqnV = size(mesh.Xv,1);
opt.neqnP = opt.neqnV;

%% Boundary conditions
if isfield(mesh,'bound')
    Null1 = ones(opt.neqnV,1);
    Null2 = ones(opt.neqnV,1);
    opt.g1=zeros(opt.neqnV,1);
    opt.g2=zeros(opt.neqnV,1);
    for i=1:size(mesh.bound,1)
        if mesh.bound(i,2) == 1
            Null1(1*mesh.bound(i,1)) = 0;
            opt.g1(1*mesh.bound(i,1)) = mesh.bound(i,3);
        elseif mesh.bound(i,2) == 2
            Null2(1*mesh.bound(i,1)) = 0;
            opt.g2(1*mesh.bound(i,1)) = mesh.bound(i,3);
        end
        if strcmp(study.BC_type,'static') == 1
            opt.g(1*mesh.bound(i,1)) = mesh.bound(i,3);
        elseif strcmp(study.BC_type,'dynamic') == 1
            opt.g_sys = @(x,y,t) [-cos(x).*sin(y).*exp(-2.*t);
                sin(x).*cos(y).*exp(-2.*t)];
        end
    end
    opt.Null2 = spdiags(Null2,0,opt.neqnV,opt.neqnV);
    opt.Null1 = spdiags(Null1,0,opt.neqnV,opt.neqnV);
end

if strcmp(study.study_type,'unsteady') == 1


    if strcmp(study.p_type,'liddriven') == 1
    study.U0 = [repmat(study.U10,opt.neqnV,1);
        repmat(study.U20,opt.neqnV,1)];
    elseif strcmp(study.p_type,'roenquist2') == 1
    study.U0 = opt.g_sys(mesh.Xv(:,2),mesh.Xv(:,3),study.t(1));
    end
end

%%

I = zeros(opt.nel*ldof*ldof,1);
J = zeros(opt.nel*ldof*ldof,1);
KE = zeros(opt.nel*ldof*ldof,1);
ME = zeros(opt.nel*ldof*ldof,1);
DE1 = zeros(opt.nel*ldof*ldof,1);
DE2 = zeros(opt.nel*ldof*ldof,1);
DD1 = zeros(opt.nel*ldof*ldof,1);
DD2 = zeros(opt.nel*ldof*ldof,1);
ntriplets = 0;

%Define GLL weights and points. Collected from study struct.
w = study.w;
xi = study.xi;

for e = 1:opt.nel
    
    fprintf('Assembling element %d of %d\r', e, opt.nel);

    nenV = mesh.IXv(:,:,e);
    nenV = nenV(:);

    edofV = zeros(ldof,1);
    indicesV = 1:length(nenV);
    edofV(indicesV) = nenV;
    %Find x and y for each node.
    xy = mesh.Xv(nenV,2:3);
    x = reshape(xy(:,1),N+1,N+1);
    y = reshape(xy(:,2),N+1,N+1);


    [me,ke,de,dde] = PnPn_twoD_element_matrices(x,y,n_GLL,w,xi);
    opt.ME(:,:,e) = me;
    opt.DE1v(:,:,e) = de{1};
    opt.DE2v(:,:,e) = de{2};


    for krow = 1:ldof
        for kcol = 1:ldof
            ntriplets = ntriplets+1;
            I(ntriplets) = edofV(krow);
            J(ntriplets) = edofV(kcol);
            KE(ntriplets) = ke(krow,kcol);
            DE1(ntriplets) = de{1}(krow,kcol);
            DE2(ntriplets) = de{2}(krow,kcol);
            DD1(ntriplets) = dde{1}(krow,kcol);
            DD2(ntriplets) = dde{2}(krow,kcol);
            % mass matrix
            ME(ntriplets) = me(krow,kcol);

        end
    end

clc
end
ind = find(I>0);
opt.K = sparse(I(ind),J(ind),KE(ind),opt.neqnV,opt.neqnV);
opt.M = sparse(I(ind),J(ind),ME(ind),opt.neqnV,opt.neqnV);
opt.DE1 = sparse(I(ind),J(ind),DE1(ind),opt.neqnV,opt.neqnV);
opt.DE2 = sparse(I(ind),J(ind),DE2(ind),opt.neqnV,opt.neqnV);
opt.DD1 = sparse(I(ind),J(ind),DD1(ind),opt.neqnV,opt.neqnV);
opt.DD2 = sparse(I(ind),J(ind),DD2(ind),opt.neqnV,opt.neqnV);
% opt.P = sparse(opt.neqnV,1);
% Initial conditions for transient problem.

end