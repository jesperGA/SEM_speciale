function [opt,study] = AssemblyQuad(mesh,opt,study)

n_GLL = study.n_GLL; %Specify number of GLL points
N = n_GLL-1;
ldof =(n_GLL^2); %One dof per gll in every element

opt.nel = size(mesh.IX,3);
opt.neqn = size(mesh.X,1);

%% Boundary conditions
if isfield(mesh,'bound')
    Null = ones(opt.neqn,1);
    opt.g=zeros(opt.neqn,1);
    for i=1:size(mesh.bound,1)
        Null(1*mesh.bound(i,1)-(1-mesh.bound(i,2))) = 0;
        opt.g(1*mesh.bound(i,1)-(1-mesh.bound(i,2))) = mesh.bound(i,3);
    end
    opt.Null = spdiags(Null,0,opt.neqn,opt.neqn);
end
%%

I = zeros(opt.nel*ldof*ldof,1);
J = zeros(opt.nel*ldof*ldof,1);
KE = zeros(opt.nel*ldof*ldof,1);
ME = zeros(opt.nel*ldof*ldof,1);

ntriplets = 0;


%Define GLL weights and points. Collected from study struct.
w = study.w;
xi = study.xi;

for e = 1:opt.nel

    nen = mesh.IX(:,:,e);
    nen = nen(:);

    % generate edof
    % for i = 1:length(nen)
    %     edof(2*i-1) = 2*nen(i)-1;
    %     edof(2*i-0) = 2*nen(i)-0;
    % end
    %Calculate dofs for each elements.
    edof = zeros(ldof,1);
    indices = 1:length(nen);
    edof(indices) = nen;

    %Find x and y for each node.
    xy = mesh.X(nen,2:3);
    x = reshape(xy(:,1),N+1,N+1);
    y = reshape(xy(:,2),N+1,N+1);
    % matID = mesh.IX(e,end);


    [me,ke] = twoD_helmholtz(x,y,[],[],n_GLL,w,xi);

    for krow = 1:ldof
        for kcol = 1:ldof
            ntriplets = ntriplets+1;
            I(ntriplets) = edof(krow);
            J(ntriplets) = edof(kcol);
            KE(ntriplets) = ke(krow,kcol);
            % mass matrix
            ME(ntriplets) = me(krow,kcol);

        end

    end

end
ind = find(I>0);
opt.K = sparse(I(ind),J(ind),KE(ind),opt.neqn,opt.neqn);
opt.M = sparse(I(ind),J(ind),ME(ind),opt.neqn,opt.neqn);
opt.P = sparse(opt.neqn,1);

% Initial conditions for transient problem.

end