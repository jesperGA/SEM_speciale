function [me,ke] = twoD_helmholtz(x,y,E,dens,me_Gauss,w,xi)

ldof = numel(x);

ke = zeros(me_Gauss);
me = zeros([me_Gauss,1]);
L = x(end)-x(1);
D = 1111;



N=me_Gauss-1;
% Calc all langragian derivatives:
for i = 1:N+1
    for j = 1:N+1
        dl(i,j) = dlagrange(xi,N,i,j);
    end
end

%Calculate contributions to the jacobian.
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
% 
% J = [sum(sum(xr)) sum(sum(yr))
%     sum(sum(xs)) sum(sum(ys))];
% det(J)

dJ = xr.*ys-xs.*yr;
% det(J)

% for i = 1:N+1
%     for j = 1:N+1
%         temp = 0;
%         for k = 1:N+1
%             %Calculate lagrange derivatives
%             dl1 = dl(j,k); dl2 = dl(i,k);
%             %Compute sum
%             temp = temp +w(k)*iJ*dl1*dl2;
%         end
%         % ke(i,j) =  exp(x(i))*temp;
%         ke(i,j) =  temp;
%     end
% end
% 
for m = 1:ldof
    for n = 1:ldof

    me(i) = J * w(i);
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

function [J] = Jac2D(xi,eta,x,y,dl)

end
