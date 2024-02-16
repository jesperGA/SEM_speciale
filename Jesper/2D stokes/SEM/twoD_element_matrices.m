function [me,ke,de,mhat] = twoD_element_matrices(x,y,xp,yp,E,dens,n_gll,w,wp,xi,zeta)

% ldof = numel(x);

ke = zeros([n_gll*n_gll,n_gll*n_gll]);
me = ke; %de{1} = ke;de{2} = ke;


L1 = x(end,1)-x(1,1);
L2 = y(1,end)-y(1,1);



N=n_gll-1;
N_gl = numel(wp);
mhat = zeros(N_gl,N_gl);
% Calc all langragian derivatives: MAYBE CALCULATE EVEN FURTHER BACK? iS
% THE SAME FOR ALL ELEMENTS.
for i = 1:N+1
    for j = 1:N+1
        dl(i,j) = dlagrange(xi,N,j,i);
    end
end
% dl = dl;


[dJ,xr,xs,yr,ys] = Jac2D(x,y,dl,N);
gradl = zeros([n_gll,n_gll,n_gll,n_gll,2]);
kroen = eye(n_gll,n_gll);
%%
%Calculate the gradient operator of the test function. 4th order tensor
%with 2 spatial components
for p=1:n_gll
    for q = 1:n_gll
        for m=1:n_gll
            for n = 1:n_gll
                gradl(p,q,m,n,1) = ys(p,q)*dl(p,m)*kroen(q,n)-yr(p,q)*kroen(p,m)*dl(q,n);
                gradl(p,q,m,n,2) = xr(p,q)*kroen(p,m)*dl(q,n)-xs(p,q)*dl(p,m)*kroen(q,n);
            end
        end
    end
end
%Compute ke and me.
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

%%
%CALCULATE dlP, the derivative of the lagrange polynomial from v in the p
%grid.
% LagrangeFormInterpolation(x_data,y_data,x_fit);
for i =1:N+1
    for j=1:N_gl
        I_tilde(j,i) = wp(j).*cardinalP(xi,i,zeta(j));
    end
    dlP(:,i) =wp.*LagrangeFormInterpolation(xi,dl(:,i).',zeta).';
end

%Only for rectlinear elements.
for i = 1:N_gl
    for j = 1:N_gl
        for m = 1:N_gl
            for n = 1:N_gl
                row = (i-1)*N_gl + j;
                col = (m-1)*N_gl + n;
                mhat(row,col) = wp(i)*wp(j)*abs(dJ(i,j))*kroen(i,m)*kroen(j,n);
            end
        end
    end
end
%CALC D only for rectilinear elements. For deformed element use dJ for P
%mesh.
de{1} = kron(L2/2.*I_tilde,dlP);de{2} = kron(L1/2.*dlP,I_tilde);
% for i = 1:N_gl
%     for j = 1:N+1
%         for m = 1:N_gl
%             for n = 1:N+1
%                 % Compute row and column indices for mapping into C
%                 row = (i-1)*(N_gl) + m;
%                 col = (j-1)*(N+1) + n;
%
%                 % Compute the Kronecker product and insert into C
%                 de{1}(row, col) = dJ(i,m)*I_tilde(i, j) * dlP(m, n);
%                 de{2}(row, col) = dJ(m,i)*dlP(i,j)*I_tilde(m,n);
%             end
%         end
%     end
% end
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
%The derivative of the ith lagrange polynomial to the kth xi.
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
%Calculated according to Rønquist (2.xx)
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
function fit = LagrangeFormInterpolation(knots,ydata,t)
%LagrangeFormInterpolation: CalculateS the values of the interpolating
%                           polynomial in Lagrange form
%
%   fit = LagrangeFormInterpolation(knots,ydata,t)
%
% Input:
%   knots=[x0 x1 ... xn]   is a row of n+1 knot-values
%   ydata=[y0 y1 ... yn]   is a row of the corresponding n+1 y-values
%   t=[t1 ... tm]          is a row of all the m values that the inter-
%                          polating polynomial is to be evaluated in
% Output:
%   fit=[P(t1) ... P(tm)]  a row with the m function values of the
%                          interpolating polynomial

m=length(knots);
if (m~=length(ydata) )
    disp('Der skal være lige mange knuder og ydata');
    fit=NaN; return;
end

cardinals=zeros(m,length(t));
for i=1:m
    cardinals(i,:)=cardinalP(knots,i,t)';
end
Y=repmat(ydata',1,length(t));
fit=sum(cardinals.*Y);

% --------- Subfunction ----------------------------------
end
