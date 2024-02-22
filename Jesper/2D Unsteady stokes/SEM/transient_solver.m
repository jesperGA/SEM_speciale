function [opt] = transient_solver(opt,study,mesh)

pref_dof = mesh.pref_dof;
dt = study.dt;
ndt = study.nt;

xx = mesh.Xv(:,2);yy = mesh.Xv(:,3);
xxp = mesh.Xp(:,2);yyp = mesh.Xp(:,3);

if strcmp(study.p_type,'roenquist') == 1
    opt.P1 = 2+pi*cos(pi*xx).*sin(pi*yy);
    opt.P2 = pi*sin(pi*xx).*cos(pi*yy);
elseif strcmp(study.p_type,'bercover') == 1
    % p1 = @(x,y) 128*(x.^2.*(x-1).*12.*(2.*y-1)+2.*(y-1).*(2.*y-1).*y.*(12.*x.^2-12.*x+2));
    % opt.P1 = p1(xx,yy);
    % opt.P2 = p1(yy,xx);
    % opt.P1 = (1536 .* xx .^ 2 .* (xx - 1) .^ 2 .* (2 .* yy - 1)) + (256 .*...
    %     (yy - 1) .* (2 .* yy - 1) .* yy .* (12 .* xx .^ 2 - 12 .* xx + 2));
    % opt.P2 = (1536 .* yy .^ 2 .* (yy - 1) .^ 2 .* (2 .* xx - 1)) + (256 .*...
    %     (xx - 1) .* (2 .* xx - 1) .* xx .* (12 .* yy .^ 2 - 12 .* yy + 2));

    F1 = @(X1,X2) 128 .* (X1.^2 .* (X1 - 1) .^2 .*12 .* (2 .* X2 - 1) + ...
        2 .* (X2 - 1) .* (2 .* X2 - 1) .* X2 .* (12 * X1 .^ 2 - 12 .* X1 + 2))+X2-1/2;
    F2 = @(X1,X2)   F1(X2,X1);

    opt.P1 = F1(xx,yy);
    opt.P2 = F2(xx,yy);
else
    error('Load case not defined')
end

if strcmp(study.int_type,'BDFk') == 1
    gamma = [1, 1 , 1];
    b0 = gamma(2)/gamma(1);
    a_k = gamma(3:end)./gamma(1);
    %Only for k=1 order in BDF-k [b0 = gamma(2)/gamma(1), a_k =gamma(3:end)./gamma(1)];


    % DEFINE BETA VARIABLES FOR BDF.
    beta = [flip(a_k)/b0, 1/b0];

    H = beta(end)/dt*opt.M+opt.M;

    %SET UP GLOBAL MATRICES FOR BCs.
    null_sys = [opt.Null1, zeros(opt.neqnV),zeros(opt.neqnV,opt.neqnP);
        zeros(opt.neqnV),opt.Null2,zeros(opt.neqnV,opt.neqnP);
        zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnV),speye(opt.neqnP,opt.neqnP)];
    g_sys =  [opt.g1;opt.g2;zeros(opt.neqnP,1)];


    %Assembly LHS(sys_mat) and part of RHS (B_mat) of the direct equation
    sys_mat = [H,zeros(opt.neqnV,opt.neqnV),-opt.DE1.';
        zeros(opt.neqnV,opt.neqnV), H,-opt.DE2.';
        -opt.DE1,-opt.DE2,zeros(opt.neqnP,opt.neqnP)];
    B_mat = [opt.M,zeros(opt.neqnV,opt.neqnV),zeros(opt.neqnV,opt.neqnP);
        zeros(opt.neqnV,opt.neqnV),opt.M,zeros(opt.neqnV,opt.neqnP);
        zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnP)];


    free = diag(null_sys);
    sys_org = sys_mat;

    sys_mat = null_sys'*sys_mat*null_sys-(null_sys-speye(size(null_sys)));
    opt.sys_mat = sys_mat;

    BCs = sys_org*g_sys;

    U(:,1) = study.U0;

    for i = 2:ndt

        Un = [U(1:2*opt.neqnV);zeros(opt.neqnP,1)];


        P = B_mat*([opt.P1;opt.P2;zeros(opt.neqnP,1)]-beta(end-1)/dt*Un);
        P = P-BCs;
        P(~free) = 0;
        P(find(g_sys)) = g_sys(find(g_sys));
        if strcmp(study.solve_type,'direct') == 1

            U = sys_mat \ (P);

        elseif strcmp(study.solve_type,'uzawa') == 1



        end

        opt.Pr(:,i) = U(end-opt.neqnP+1:end);
        opt.Pr(:,i) = opt.Pr(:,i)-opt.Pr(pref_dof);
        opt.U(:,i) = U(1:2*opt.neqnV);
    end
else
    disp('Time integration scheme not implemented')
end
end