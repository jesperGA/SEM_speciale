function [A,B] = elementMatrix2D(x,y,xis,w,D,N);

    ldof=(N+1)^2;
    A = zeros(ldof,ldof);
    B = zeros(ldof,ldof);
    % fe = zeros(ldof,1);
    % 
    % 
    % xr = zeros(N+1,N+1);
    % xs = zeros(N+1,N+1);
    % yr = zeros(N+1,N+1);
    % ys = zeros(N+1,N+1);
    % J = zeros(2,2);
    % 
    % for i = 1:N+1
    %     xi = xis(i);
    %     for j = 1:N+1
    %        eta = xis(j);
    %       [J,N] = shapeFunc(x,y,eta,xi);
    %        Bb(1:3,:) = Bm(1:3,:);
    %        dJ = det(J);
    %        kb = kb+Bb'*Dm*Bb*kbw(i)*kbw(j)*dJ;
    %        fe = fe + N(1,:)'*q*kbw(i)*kbw(j)*dJ;
    %     end
    % end
    % 
    % for p = 1:N+1
    %     for q = 1:N+1
    %         for m = 1:N+1
    %             xr(p,q) = xr(p,q) + h(p,m)*x(m,q);
    %             xs(p,q) = xs(p,q) + h(q,m)*x(p,m);
    %             yr(p,q) = yr(p,q) + h(p,m)*y(m,q);
    %             ys(p,q) = ys(p,q) + h(q,m)*y(p,m);
    %         end
    %         J(p,q) = xr(p,q)*ys(p,q) - xs(p,q)*yr(p,q);
    %     end
    % end
    % 
    % % Initialize element matrices
    % A = zeros(N+1, N+1);
    % B = zeros(N+1, N+1);
    % J = (x(end)-x(1))/2;
    % k=1

    % for i = 1 : N+1
    %     for j = 1 : N+1
    %             dlj_xi_i(i,j) = DerivativeLagrangePoly(N, xis, i, j);
    %     end
    % end
    % [xis,w,h]=GetGLL(N+1)
    % [x,w]=lglnodes(N)
    % dlj_xi_i
    % 
    % for i = 1 : N+1
    %     for j = 1 : N+1
    %         d=i==j; % Kronecker delta
    %         B(i,j) = B(i,j) + w(j) * J * d;
    %         for k = 1 : N+1
    %             dli_xi_k = DerivativeLagrangePoly(N, xis, k, i)
    %             dlj_xi_k = DerivativeLagrangePoly(N, xis, k, j);
    %             A(i,j) = A(i,j) + exp(x(1)+(xis(k)+1)*J) * w(k) * dlj_xi_k * dli_xi_k * 1 / J;
    %         end
    %     end
    % end
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

function [J,N] = shapeFunc(x,y,eta,xi)
J = zeros(2,2);
N = zeros(2,length(x)*2);

% The nodal shape functions
NN = zeros(4,1);
NN(1,1) = ((1 - xi) * (1 - eta)) / 4;
NN(2,1) = ((1 + xi) * (1 - eta)) / 4;
NN(3,1) = ((1 + xi) * (1 + eta)) / 4;
NN(4,1) = ((1 - xi) * (1 + eta)) / 4;

% The shape functions for the displ. and two rotations
for i = 1:4
    N(1,i*3-2) = NN(i,1);
    N(2,i*3-1) = NN(i,1);
    N(3,i*3) = NN(i,1);
end

% Jacobian
J(1,1) = 1/4*[(eta-1) (1-eta) (1+eta) (-1-eta)]*x.';
J(1,2) = 1/4*[(eta-1) (1-eta) (1+eta) (-1-eta)]*y.';
J(2,1) = 1/4*[(xi-1) (-1-xi) (1+xi) (1-xi)]*x.';
J(2,2) = 1/4*[(xi-1) (-1-xi) (1+xi) (1-xi)]*y.';

% inverse of J
gamma = inv(J);

% Derivative of shapefunction in local space
dN = 1/4*[(eta-1) (1-eta) (1+eta) (-1-eta)
          (xi-1) (-1-xi) (1+xi) (1-xi)];

end