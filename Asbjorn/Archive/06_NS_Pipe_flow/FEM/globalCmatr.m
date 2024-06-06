function C_global = globalCmatr(ldof, nel, neqn_u, Be,  u1, u2, D_hat, L1, L2, IX)
    %% C Matrix (Convection)

    I =  zeros(nel * ldof * ldof, 1);
    J =  zeros(nel * ldof * ldof, 1);
    CE = zeros(nel * ldof * ldof, 1);
    ntriplets = 0;

    % Loop over elements and integrate
    for e = 1:nel

        % Get element data
        edof = reshape(IX(:,:,e), [], 1);

        C = calculateCMatr(L1(e), L2(e), ldof, Be(:,:,e), u1(edof), u2(edof), D_hat);

        % add to global system (I,J,[SE])
        for krow = 1:ldof
            for kcol = 1:ldof
                ntriplets = ntriplets+1;
                I(ntriplets) = edof(krow);
                J(ntriplets) = edof(kcol);
                CE(ntriplets) = C(krow,kcol);
            end
        end
    end

    ind = find(I>0);
    C_global = sparse(I(ind),J(ind),CE(ind),neqn_u,neqn_u);
end