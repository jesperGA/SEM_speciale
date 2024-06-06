function Cv = globalCmatr_ADJ(I, J, BD1, BD2, numElems, u1, u2, ldof, neqn)
    % Compute the global convection matrix for given velocity fields u1 and u2.
    % Inputs are precomputed indices I, J and matrix coefficients BD1, BD2.
    % The output is a sparse matrix C representing the convection effect in the PDE.
    neqn_u = neqn(1);neqnP = neqn(2);
    convectionEntries = zeros(size(I));
    index = 0;
    for e = 1:numElems
        elementIndices = (index + 1):(index + ldof^2);
        elementDofs = I(elementIndices);

        % Compute the contribution of each element to the convection matrix
        convectionEntries(elementIndices) = (u1(elementDofs) .* BD1(elementIndices) + u2(elementDofs) .* BD2(elementIndices));
        index = index + ldof^2;
    end

    % Assemble the sparse global convection matrix from computed entries
    C_temp = sparse(I, J, convectionEntries, neqn_u, neqn_u);

    Cv =[C_temp*u1;
    C_temp*u2
    sparse(neqnP,1)];
end