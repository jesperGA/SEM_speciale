function [me,ke] = oneDbar_SE(x,E,dens,me_Gauss,w,xi)

ke = zeros(me_Gauss);
me = zeros([me_Gauss,1]);
L = x(end)-x(1);

J = L/2; iJ = J^-1;

N=me_Gauss-1;

for i = 1:N+1
    for j = 1:N+1
        temp = 0;
        for k = 1:N+1
            %Calculate lagrange derivatives
            dl1 = dlagrange(xi,N,j,k); dl2 = dlagrange(xi,N,i,k);
            %Compute sum
            %With continous coeffecient:
            temp = temp +exp(x(1)+(xi(k)+1)*J)*w(k)*iJ*dl1*dl2;

            %With discrete coeffecient, 
        end
        % ke(i,j) =  exp(x(i))*temp;
        ke(i,j) =  temp;
    end
end

for i = 1:N+1
    me(i) = J * w(i);
end
me = speye(N+1).*me;


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
