function [Cv] = adv_mat_assembly(IXv,ME,DE1v,DE2v,n_GLL,neqn,U)
neqnV = neqn(1);
neqnP = neqn(2);

nel = size(ME,3);

ldof = (n_GLL)^2;

U1 = U(1:neqnV);U2 = U(neqnV+1:end);

ntriplets = 0;

I = zeros(nel*ldof*ldof,1);
J = zeros(nel*ldof*ldof,1);
CE = zeros(nel*ldof*ldof,1);

for e = 1:nel
    nen = IXv(:,:,e);
    nen = nen(:);

    edof = zeros(ldof,1);
    indices = 1:length(nen);
    edof(indices) = nen;

    v1 = sparse(1:ldof,1:ldof,U1(edof),ldof,ldof);
    v2 = sparse(1:ldof,1:ldof,U2(edof),ldof,ldof);
    
    me = ME(:,:,e);
    de1 = DE1v(:,:,e);
    de2 = DE2v(:,:,e);

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
C_temp = sparse(I(ind),J(ind),CE(ind),neqnV,neqnV);
Cv =[C_temp*U1;
    C_temp*U2
    sparse(neqnP,1)];


end