function [A,B] = elementMatrixBar(x,xi,w);

    N=length(x)-1;

    % Initialize element matrices
    A = zeros(N+1, N+1);
    B = zeros(N+1, N+1);
    J = (x(end)-x(1))/2;

    for i = 1 : N+1
        for j = 1 : N+1
            d=i==j; % Kronecker delta
            B(i,j) = B(i,j) + w(j) * J * d;
            for k = 1 : N+1
                dli_xi_k = DerivativeLagrangePoly(N, xi, k, i);
                dlj_xi_k = DerivativeLagrangePoly(N, xi, k, j);
                A(i,j) = A(i,j) + exp(x(1)+(xi(k)+1)*J) * w(k) * dlj_xi_k * dli_xi_k * 1 / J;
            end
        end
    end
end


function dl = DerivativeLagrangePoly(N, xi, i, k)
    dl=0;
    for j=1:N+1
        d_ij = d_ij_function(i-1,j-1,N,xi);
        lk = CardinalPolynomial(xi,k,xi(j));
        dl = dl + d_ij*lk;
    end
    
end

function d_ij = d_ij_function(i,j,N,xi)
    d_ij = 0;
    if i==0 && j==0
        d_ij = -1/4*N*(N+1);
    elseif i>=0 && i<=N && j>=0 && j<=N && i~=j
        L_N_xi_i = LegendrePoly(N,xi(i+1));
        L_N_xi_j = LegendrePoly(N,xi(j+1));
        d_ij = L_N_xi_i/L_N_xi_j*(1/(xi(i+1)-xi(j+1)));
    elseif i>=1 && i==j && j<=N-1
        d_ij = 0;
    elseif i==N && j==N
        d_ij = 1/4*N*(N+1);
    end

end

function L_N = LegendrePoly(N,xi)
    L=zeros(N+1,1);
    L(1) = 1;
    L(2) = xi;
    for i=3:N+1
        n=i-1;
        L(i) = 1/n*((2*n-1)*xi*L(i-1)-(n-1)*L(i-2));
    end
    L_N=L(end);
end