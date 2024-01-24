function [me,ke] = twoD_helmholtz(x,y,E,dens,n_gll,w,xi)

ldof = numel(x);

ke = zeros([n_gll*n_gll,n_gll*n_gll]);
me = ke;



N=n_gll-1;
% Calc all langragian derivatives: MAYBE CALCULATE EVEN FURTHER BACK? iS
% THE SAME FOR ALL ELEMENTS.
for i = 1:N+1
    for j = 1:N+1
        dl(i,j) = dlagrange(xi,N,j,i);
    end
end

[dJ,xr,xs,yr,ys] = Jac2D(x,y,dl,N);
gradl = zeros([n_gll,n_gll,n_gll,n_gll,2]);
kroen = eye(n_gll,n_gll);

for p=1:n_gll
    for q = 1:n_gll
        for m=1:n_gll
            for n = 1:n_gll
                gradl(p,q,m,n,:) = [ys(p,q)*dl(p,m)*kroen(q,n)-yr(p,q)*kroen(p,m)*dl(q,n)
                                    xr(p,q)*kroen(p,m)*dl(q,n)-xs(p,q)*dl(p,m)*kroen(q,n)];
            end
        end
    end
end



%Calculate contributions to the jacobian.


for i = 1:n_gll
    for j = 1:n_gll
        for m = 1:n_gll
            for n = 1:n_gll
                row = (i-1)*n_gll + j;
                col = (m-1)*n_gll + n;
                for p = 1:n_gll
                    for q=1:n_gll
                        ke(row,col) =ke(row,col) + w(p)*w(q)*(1/dJ(p,q))*dot(gradl(p,q,i,j,:),gradl(p,q,m,n,:));
                    end
                end

                me(row,col) = w(i)*w(j)*abs(dJ(i,j))*kroen(i,m)*kroen(j,n);
            end
        end
    end
end




end

function L_out = legendre_pol(N,xi)
%Calculated acording to Igel PPT page 44
L = zeros([N+1,1]);
L(1) = 1;
L(2) = xi;

for i = 3:N+1
    n=i-1;
    L(i) = 1/n*((2*n-1)*xi*L(i-1)-(n-1)*L(i-2));
end

L_out = L(end);

end

function dij = dij_calc(n,xi,i,j)
% xi is a vector with 2 scalar values. xi = [xi_i, xi_j]
if i == n && j==n
    dij = 1/4*n*(n+1);
elseif i == 0 && j==0
    dij = -1/4*n*(n+1);
elseif i==j
    dij = 0;
else
    li = legendre_pol(n,xi(1));
    lj = legendre_pol(n,xi(2));
    dij = (li/lj)*(1/(xi(1)-xi(2)));
end

end

function dl = dlagrange(xi_vec,N,k,i)
dl = 0;
for j=1:N+1
    dij = dij_calc(N,xi_vec([i,j]),i-1,j-1);
    l = cardinalP(xi_vec,k,xi_vec(j));
    dl = dl + dij*l;
end

end

function [l]=cardinalP(x,i,t)

l=ones(size(t)); % Accepting both a row vector and a column vector
m=length(x);
% for j=[1:i-1 i+1:m]
for  j = 1:m
    if j==i
    else
        l=l.*(t-x(j))./(x(i)-x(j));
    end
end
end

function [dJ,xr,xs,yr,ys] = Jac2D(x,y,dl,N)
xr = zeros(N+1,N+1);xs = xr;yr = xr; ys = xr; % Initialize vectors
for p = 1:N+1
    for q=1:N+1
        for m=1:N+1
            xr(p,q) = xr(p,q) + dl(p,m)*x(m,q);
            xs(p,q) = xs(p,q) + dl(q,m)*x(p,m);
            yr(p,q) = yr(p,q) + dl(p,m)*y(m,q);
            ys(p,q) = ys(p,q) + dl(q,m)*y(p,m);
        end
    end
end

dJ = xr.*ys-xs.*yr;
end
