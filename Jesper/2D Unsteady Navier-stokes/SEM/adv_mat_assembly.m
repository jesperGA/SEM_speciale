function [Cv] = adv_mat_assembly(IXv, ME_DE1, ME_DE2,n_GLL,neqn,U,S)
neqnV = neqn(1);
neqnP = neqn(2);

nel = size(ME_DE1,3);

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

    ce = v1 * ME_DE1(:,:,e) + v2 * ME_DE2(:,:,e);

    % Create all combinations of krow and kcol
    [krow, kcol] = meshgrid(1:ldof, 1:ldof);
    krow = krow(:);
    kcol = kcol(:);

    % Vectorize the assignment to I, J, CE
    I(ntriplets+1:ntriplets+ldof^2) = edof(krow);
    J(ntriplets+1:ntriplets+ldof^2) = edof(kcol);
    CE(ntriplets+1:ntriplets+ldof^2) = ce(sub2ind([ldof, ldof], krow, kcol));

    % Update ntriplets
    ntriplets = ntriplets + ldof^2;

end
ind = find(I>0);
C_temp = sparse(I(ind),J(ind),CE(ind),neqnV,neqnV);
Cv =[C_temp*U1;
    C_temp*U2
    sparse(neqnP,1)];


end