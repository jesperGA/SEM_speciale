function [D1, D2] = calculateDMatrices(xy, N, xi, w, zeta, omega)
    
    L1 = max(xy(:,1)) - min(xy(:,1));
    L2 = max(xy(:,2)) - min(xy(:,2));

    for j = 1:N+1
        for i = 1:length(zeta)
            pi_j_zeta_i(i,j) = Lagrange(j, xi, zeta(i));
            D_j_zeta_i(i,j) = derivativeLagrange(j, xi, zeta(i));
        end
    end
    I_tilde = omega .* pi_j_zeta_i;
    D_tilde = omega .* D_j_zeta_i;
    
    D1 = L2/2 * kron(I_tilde, D_tilde);
    D2 = L1/2 * kron(D_tilde, I_tilde);
end

function dl = derivativeLagrange(j, xi, x)
% Calculates the j'th Lagrange interpolation polynomial's derivative at a point x
%
% Inputs:
%   j - Derivative of j'th polynomium
%   xi - vector of xi coordinates
%   x - the point at which to evaluate the derivative
%
% Output:
%   dl - the value of the derivative at point z

dl = 0;
k = length(xi);
for i = 1:k
    if i ~= j
        term = 1 / (xi(j) - xi(i));
        for m = 1:k
            if m ~= j && m ~= i
                term = term * (x - xi(m)) / (xi(j) - xi(m));
            end
        end
        dl = dl + term;
    end
end
end

function l = Lagrange(j, xi, x)
% Computes the j-th Lagrange basis polynomial at a point x.
%
% Inputs:
%   j - j'th Lagrange polynomium
%   xi - vector of xi coordinates
%   x - the point at which to evaluate the basis polynomial
%
% Output:
%   l - the value of the i-th basis polynomial at point z

k = length(xi);
l = 1;
for m = 1:k
    if m ~= j
        l = l * (x - xi(m)) / (xi(j) - xi(m));
    end
end
end