function [opt,study] = Assembly1bar(mesh,opt,study)

n_GLL = study.n_GLL; %Specify number of GLL points
ldof =n_GLL*1;

opt.nel = size(mesh.IX,1);
opt.neqn = size(mesh.X,1);

I = zeros(opt.nel*ldof*ldof,1);
J = zeros(opt.nel*ldof*ldof,1);
KE = zeros(opt.nel*ldof*ldof,1);
ME = zeros(opt.nel*ldof*ldof,1);

ntriplets = 0;


%Define GLL weights and points. Collected from study struct.
w = study.w;
xi = study.xi;

for e = 1:opt.nel

    nen = mesh.IX(e,2:end-1);
    edof = nen; %Only true for 1D.
    x = mesh.X(nen,2);

    matID = mesh.IX(e,end);

    [me,ke] = oneDbar_SE(x,[],[],n_GLL,w,xi);
    % [me,ke] = trapz_SE(x,[],[],n_GLL,w,xi);


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

% Initial conditions for transient problem.

end