function C = calculateCMatr(L1, L2, ldof, B, u1, u2, D_hat)

    I = eye(length(D_hat));
    D1 = kron(I, D_hat) * 2/L1;
    D2 = kron(D_hat, I) * 2/L2;

    V1 = sparse(1:ldof,1:ldof,u1,ldof,ldof);
    V2 = sparse(1:ldof,1:ldof,u2,ldof,ldof);

    C = B * (V1* D1 + V2 * D2);

end