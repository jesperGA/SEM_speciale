function [Cv] = adv_mat_assembly(mesh,opt,study,U)

n_GLL = study.n_GLL;
ldof = (n_GLL)^2;

U1 = U(1:opt.neqnV);U2 = U(opt.neqnV+1:end);

ntriplets = 0;

I = zeros(opt.nel*ldof*ldof,1);
J = zeros(opt.nel*ldof*ldof,1);
CE = zeros(opt.nel*ldof*ldof,1);

for e = 1:opt.nel
    nen = mesh.IXv(:,:,e);
    nen = nen(:);

    edof = zeros(ldof,1);
    indices = 1:length(nen);
    edof(indices) = nen;

    v1 = sparse(1:ldof,1:ldof,U1(edof),ldof,ldof);
    v2 = sparse(1:ldof,1:ldof,U2(edof),ldof,ldof);
    
    me = opt.ME(:,:,e);
    de1 = opt.DE1v(:,:,e);
    de2 = opt.DE2v(:,:,e);

    ce = me*(v1*de1+v2*de2);
    
    for krow = 1:ldof
        for kcol = 1:ldof
            ntriplets = ntriplets+1;
            I(ntriplets) = edof(krow);
            J(ntriplets) = edof(kcol);
            CE(ntriplets) = ce(krow,kcol);

        end
    end

end
ind = find(I>0);
C_temp = sparse(I(ind),J(ind),CE(ind),opt.neqnV,opt.neqnV);
Cv =[C_temp*U1;
    C_temp*U2
    sparse(opt.neqnP,1)];


end