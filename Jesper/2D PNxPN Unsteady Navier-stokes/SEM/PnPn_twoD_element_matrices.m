function [me,ke,de,dde] = PnPn_twoD_element_matrices(x,y,n_gll,w,xi)

% ldof = numel(x);

ke = zeros([n_gll*n_gll,n_gll*n_gll]);
me = ke; %de{1} = ke;de{2} = ke;
de{1} = ke;de{2} = ke; %Initilization of divergence/gradient matrix for v-grid.

L1 = x(end,1)-x(1,1);
L2 = y(1,end)-y(1,1);

N=n_gll-1;
% N_gl = numel(wp);
% Calc all langragian derivatives: MAYBE CALCULATE EVEN FURTHER BACK? iS
% THE SAME FOR ALL ELEMENTS.
for i = 1:N+1
    for j = 1:N+1
        dl(i,j) = dlagrange(xi,N,j,i);
        ddl (i,j) = ddlagrange(xi,N,i,j);
        % Ip(i,j) = w(i)*w(j)*cardinalP(xi,i,xi(j));
    end
end
% dl = dl;

[dJ,xr,xs,yr,ys] = Jac2D(x,y,dl,N);
gradl = zeros([n_gll,n_gll,n_gll,n_gll,2]);
grad2l  = gradl;
kroen = eye(n_gll,n_gll);
%%
%Calculate the gradient operator of the test function. 4th order tensor
%with 2 spatial components. Machen Sie mich bitte effizienter.
%Gradient operator also used for derivative operator matrix
for p=1:n_gll
    for q = 1:n_gll
        for m=1:n_gll
            for n = 1:n_gll
                row = (p-1)*n_gll + q;
                col = (m-1)*n_gll + n;

                gradl(p,q,m,n,1) = ys(p,q)*dl(p,m)*kroen(q,n)-yr(p,q)*kroen(p,m)*dl(q,n);
                gradl(p,q,m,n,2) = xr(p,q)*kroen(p,m)*dl(q,n)-xs(p,q)*dl(p,m)*kroen(q,n);
                de{1}(row,col) = (1/dJ(p,q))*gradl(p,q,m,n,1);
                de{2}(row,col) = (1/dJ(p,q))*gradl(p,q,m,n,2);

                grad2l(p,q,m,n,1) = ys(p,q)*ddl(p,m)*kroen(q,n)-yr(p,q)*kroen(p,m)*ddl(q,n);
                grad2l(p,q,m,n,2) = xr(p,q)*kroen(p,m)*ddl(q,n)-xs(p,q)*ddl(p,m)*kroen(q,n);

                dde{1}(row,col) = (1/dJ(p,q))*grad2l(p,q,m,n,1);
                dde{2}(row,col) = (1/dJ(p,q))*grad2l(p,q,m,n,2);

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

%Calculate dJ*w*cardinal for all GLL. in doubt of the addition dJ.
for i = 1:N+1
    for j = 1:N+1
        Ip(i,j) = cardinalP(xi,i,xi(j));
    end
end
%Derivative matrix on the V-grid
deV{1} = 2/L1.*kron(Ip,dl);deV{2} = 2/L2.*kron(dl,Ip);
% de{2} = test1;de{1} = test2;

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

function l1 = ddlagrange(xi,N,k,i)
%The second derivative of the kth lagrange polynomial to the ith xi.
l1 = 0;
l2 = 0;
l3 = 1;
for ii = 1:N+1
    if ii==k
    else
        for jj = 1:N+1
            if jj == ii || jj== k
            else
                for kk = 1:N+1
                    if kk==ii || kk== jj || kk == k
                    else
                        l3 = l3*(xi(i)-xi(kk))/(xi(k)-xi(kk));
                    end
                end
                l2 = l2 + 1/(xi(k)-xi(jj))*l3;
                l3 = 1;
            end
        end
        l1 = l1+1/(xi(k)-xi(ii))*l2;
        l2 = 0;
    end
end

end

function fit = LagrangeFormInterpolation(knots,ydata,t)
%Source: Introduction to Numerical algorithms Course @ The techical
%University of Denmark
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
