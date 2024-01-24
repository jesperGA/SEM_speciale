function [A,B] = elementMatrix2D(x,y,xis,w,D,N);

    x = flipud(fliplr(x));
    y = flipud(fliplr(y));
    ldof=(N+1)^2;
    A = zeros(N+1,N+1,N+1,N+1);
    B = zeros(N+1,N+1,N+1,N+1);
    grad = zeros(N+1,N+1,N+1,N+1,2);


    % fe = zeros(ldof,1);
    % 
    % 
    xr = zeros(N+1,N+1);
    xs = zeros(N+1,N+1);
    yr = zeros(N+1,N+1);
    ys = zeros(N+1,N+1);
    J = zeros(2,2);


    for p = 1:N+1
        for q = 1:N+1
            for m = 1:N+1
                xr(p,q) = xr(p,q) + D(p,m)*x(m,q);
                xs(p,q) = xs(p,q) + D(q,m)*x(p,m);
                yr(p,q) = yr(p,q) + D(p,m)*y(m,q);
                ys(p,q) = ys(p,q) + D(q,m)*y(p,m);
            end
            % J(p,q) = xr(p,q)*ys(p,q) - xs(p,q)*yr(p,q);
        end
    end
    J = xr.*ys - xs.*yr;

    for p = 1:N+1
        for q = 1:N+1
            for m = 1:N+1
                for n = 1:N+1
                    d_qn=q==n; % Kronecker delta
                    d_pm=q==n; % Kronecker delta
                    grad(p,q,m,n,1) = ys(p,q) * D(p,m) * d_qn - yr(p,q) * d_pm * D(q,n);
                    grad(p,q,m,n,2) = xr(p,q) * d_pm * D(q,n) - xs(p,q) * D(p,m) * d_qn;
                end
            end
        end
    end

    for i = 1:N+1
        for j = 1:N+1
            for m = 1:N+1
                for n = 1:N+1
                    for p = 1:N+1
                        for q = 1:N+1
                            A(i,j,m,n) = A(i,j,m,n) + w(p) * w(q) * 1/J(p,q) * dot(grad(p,q,i,j,:),grad(p,q,m,n,:));
                        end
                    end
                    d_im=i==m; % Kronecker delta
                    d_jn=j==n; % Kronecker delta
                    B(i,j,m,n) = B(i,j,m,n) + w(i) * w(j) * J(i,j) * d_im * d_jn;
                end
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