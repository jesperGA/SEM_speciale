function [opt,study] = AssemblyQuad(mesh,opt,study)

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
        opt.g(1*mesh.bound(i,1)) = mesh.bound(i,3);
    end
    opt.Null2 = spdiags(Null2,0,opt.neqnV,opt.neqnV);
    opt.Null1 = spdiags(Null1,0,opt.neqnV,opt.neqnV);
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


ntriplets = 0;
ntripletsP = 0;
ntripletsP2 = 0;


%Define GLL weights and points. Collected from study struct.
w = study.w;
xi = study.xi;
zeta = study.zeta;
wp = study.wp;

for e = 1:opt.nel

    nenV = mesh.IXv(:,:,e);
    nenV = nenV(:);

    nenP = mesh.IXp(:,:,e);
    nenP = nenP(:);

    % generate edof
    % for i = 1:length(nen)
    %     edof(2*i-1) = 2*nen(i)-1;
    %     edof(2*i-0) = 2*nen(i)-0;
    % end
    %Calculate dofs for each elements.
    edofP = zeros(ldofp,1);
    indicesP = 1:length(nenP);
    edofP(indicesP) = nenP;

    edofV = zeros(ldofv,1);
    indicesV = 1:length(nenV);
    edofV(indicesV) = nenV;
    %Find x and y for each node.
    xy = mesh.Xv(nenV,2:3);
    x = reshape(xy(:,1),Nv+1,Nv+1);
    y = reshape(xy(:,2),Nv+1,Nv+1);

    % x = fliplr(x);y=fliplr(y);
    % matID = mesh.IX(e,end);

    xyp = mesh.Xp(nenP,2:3);
    xp = reshape(xyp(:,1),Np+1,Np+1);
    yp = reshape(xyp(:,2),Np+1,Np+1);

    [me,ke,de,mhat] = twoD_element_matrices(x,y,xp,yp,[],[],n_GLL,w,wp,xi,zeta);

    for krow = 1:ldofv
        for kcol = 1:ldofv
            ntriplets = ntriplets+1;
            I(ntriplets) = edofV(krow);
            J(ntriplets) = edofV(kcol);
            KE(ntriplets) = ke(krow,kcol);
            % mass matrix
            ME(ntriplets) = me(krow,kcol);

        end
    end

    for prow=1:ldofp
        for kcol = 1:ldofv
            ntripletsP =ntripletsP+1;
            Ip(ntripletsP) = edofP(prow);
            Jp(ntripletsP) = edofV(kcol);
            DE1(ntripletsP) = de{1}(prow,kcol);
            DE2(ntripletsP) = de{2}(prow,kcol);
        end
    end


    for prow = 1:ldofp
        for pcol = 1:ldofp
            ntripletsP2 =ntripletsP2+1;
            Ip2(ntripletsP2) = edofP(prow);
            Jp2(ntripletsP2) = edofP(pcol);
            Mh(ntripletsP2) = mhat(prow,pcol);
        end
    end
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