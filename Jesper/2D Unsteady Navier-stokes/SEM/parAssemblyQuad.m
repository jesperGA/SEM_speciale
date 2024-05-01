function [opt,study] = parAssemblyQuad(mesh,opt,study)

n_GLL = study.n_GLL; %Specify number of GLL points
n_GL = study.n_GL;
Nv = n_GLL-1;
Np = n_GL-1;
ldofv =(n_GLL^2); %One dof per gll in every element
ldofp = (n_GL^2);

opt.nel = size(mesh.IXp,3);
opt.neqnV = size(mesh.Xv,1);
opt.neqnP = size(mesh.Xp,1);

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

I = zeros(opt.nel*ldofv*ldofv,1);
J = zeros(opt.nel*ldofv*ldofv,1);
KE = zeros(opt.nel*ldofv*ldofv,1);
ME = zeros(opt.nel*ldofv*ldofv,1);

Ip = zeros(opt.nel*ldofv*ldofp,1);
Jp = zeros(opt.nel*ldofv*ldofp,1);
DE1 = zeros(opt.nel*ldofv*ldofp,1);
DE2 = zeros(opt.nel*ldofv*ldofp,1);

Mh = zeros(opt.nel*ldofp*ldofp,1);
Ip2 = Mh;
Jp2 = Mh;


% ntriplets = 0;
% ntripletsP = 0;
% ntripletsP2 = 0;


%Define GLL weights and points. Collected from study struct.
w = study.w;
xi = study.xi;
zeta = study.zeta;
wp = study.wp;

IXv = mesh.IXv;
IXp = mesh.IXp;
Xv = mesh.Xv;
Xp = mesh.Xp;

elementKE = cell(opt.nel, 1);
elementME = cell(opt.nel, 1);
elementDE1 = cell(opt.nel, 1);
elementDE2 = cell(opt.nel, 1);
elementMh = cell(opt.nel, 1);
fprintf('Computing matrices')
parfor e = 1:opt.nel

    nenV = IXv(:,:,e);
    nenV = nenV(:);

    nenP = IXp(:,:,e);
    nenP = nenP(:);

    xy = Xv(nenV,2:3);
    x = reshape(xy(:,1),Nv+1,Nv+1);
    y = reshape(xy(:,2),Nv+1,Nv+1);

    xyp = Xp(nenP,2:3);
    xp = reshape(xyp(:,1),Np+1,Np+1);
    yp = reshape(xyp(:,2),Np+1,Np+1);

    [me,ke,deP,mhat,deV] = twoD_element_matrices(x,y,xp,yp,[],[],n_GLL,w,wp,xi,zeta);
    elementME{e} = me;
    elementKE{e} = ke;
    elementDE1{e} = deP{1};
    elementDE2{e} = deP{2};
    elementMh{e} = mhat;

    ME_loc(:,:,e) = me;
    DE1v(:,:,e) = deV{1};
    DE2v(:,:,e) = deV{2};

end
fprintf('- DONE\n')
opt.ME = ME_loc;
opt.DE1v = DE1v;
opt.DE2v = DE2v;
for e = 1:opt.nel
    fprintf('Assembling element %d of %d\r', e, opt.nel);

    me = elementME{e};
    ke = elementKE{e};
    deP{1} = elementDE1{e};
    deP{2} = elementDE2{e};
    mhat = elementMh{e};

    nenV = IXv(:,:,e);
    nenV = nenV(:);

    nenP = IXp(:,:,e);
    nenP = nenP(:);

    edofP = zeros(ldofp,1);
    indicesP = 1:length(nenP);
    edofP(indicesP) = nenP;

    edofV = zeros(ldofv,1);
    indicesV = 1:length(nenV);
    edofV(indicesV) = nenV;

    % Offset to start filling the arrays from the correct position
    offset = (e - 1) * ldofv^2;

    % Create expanded grid indices for the current element
    [krow, kcol] = meshgrid(1:ldofv, 1:ldofv);
    krow = krow(:);
    kcol = kcol(:);
    linearIndices = sub2ind([ldofv, ldofv], krow, kcol);

    % Compute linear indices for storage
    indices = offset + (1:ldofv^2)';

    % Assign to I, J, KE, ME
    I(indices) = edofV(krow);
    J(indices) = edofV(kcol);
    KE(indices) = ke(linearIndices);
    ME(indices) = me(linearIndices);

    % DE1 and DE2
    offsetP = (e - 1) * ldofp * ldofv;
    [prow, kcol] = meshgrid(1:ldofp, 1:ldofv);
    indicesP = offsetP + (1:numel(prow))';
    linearIndices = sub2ind([ldofp, ldofv], prow, kcol);
    Ip(indicesP) = edofP(prow(:));
    Jp(indicesP) = edofV(kcol(:));
    DE1(indicesP) = deP{1}(linearIndices);
    DE2(indicesP) = deP{2}(linearIndices);

    % Mh
    offsetP2 = (e - 1) * ldofp^2;
    [prow, pcol] = meshgrid(1:ldofp, 1:ldofp);
    indicesP2 = offsetP2 + (1:numel(prow))';
    linearIndices = sub2ind([ldofp, ldofp], prow, pcol);
    Ip2(indicesP2) = edofP(prow(:));
    Jp2(indicesP2) = edofP(pcol(:));
    Mh(indicesP2) = mhat(linearIndices);
    clc
end
ind = find(I>0);
indp = find(Ip>0);
indp2 = find(Ip2>0);
opt.K = sparse(I(ind),J(ind),KE(ind),opt.neqnV,opt.neqnV);
opt.M = sparse(I(ind),J(ind),ME(ind),opt.neqnV,opt.neqnV);
% opt.P = sparse(opt.neqnV,1);
opt.DE1 = sparse(Ip(indp),Jp(indp),DE1(indp),opt.neqnP,opt.neqnV);
opt.DE2 = sparse(Ip(indp),Jp(indp),DE2(indp),opt.neqnP,opt.neqnV);
opt.Mh = sparse(Ip2(indp2),Jp2(indp2),Mh(indp2),opt.neqnP,opt.neqnP);
% Initial conditions for transient problem.

end