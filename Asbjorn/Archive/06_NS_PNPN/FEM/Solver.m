function opt = Solver(mesh,study,opt)

% function opt = Solver(mesh,study,p)
%
% Input
% study.analysis =  'static'
%                   ! More will come from next week !
%
% study.X with X = neig (int), omega (float(:))
%
% mesh: mesh data (X,IX,bound,....)
% opt: Matrices (K,C,M,P,Null)
%
% Output
% opt.U , opt.w , (opt.pdof)

H1 = opt.LHS_BC(1:opt.neqn_u,1:opt.neqn_u);
H2 = opt.LHS_BC(opt.neqn_u+1:2*opt.neqn_u,opt.neqn_u+1:2*opt.neqn_u);
D1 = -opt.LHS_BC(2*opt.neqn_u+1:end,1:opt.neqn_u);
D2 = -opt.LHS_BC(2*opt.neqn_u+1:end,opt.neqn_u+1:2*opt.neqn_u);

dt = study.t(2) - study.t(1);

rho = mesh.Material(1);
mu = mesh.Material(2);

opt.u1 = zeros(opt.neqn_u,length(study.t));
opt.u2 = zeros(opt.neqn_u,length(study.t));
opt.p  = zeros(opt.neqn_p,length(study.t));
g_u1 = sparse(opt.neqn_u,1);
g_u2 = sparse(opt.neqn_u,1);
g_p = sparse(opt.neqn_p,1);

gamma0 = 11/6;
alpha0 = 3;
alpha1 = -3/2;
alpha2 = 1/3;
beta0 = 3;
beta1 = -3;
beta2 = 1;

ldof = (study.N + 1) ^ 2;

[~,~,D_hat] = GetGLL(study.N+1);

u1_oldold = study.U1(mesh.X(:,2),mesh.X(:,3),0-2*dt);
u1_old = study.U1(mesh.X(:,2),mesh.X(:,3),0-dt);
u1 = study.U1(mesh.X(:,2),mesh.X(:,3),0);
u2_oldold = study.U2(mesh.X(:,2),mesh.X(:,3),0-2*dt);
u2_old = study.U2(mesh.X(:,2),mesh.X(:,3),0-dt);
u2 = study.U2(mesh.X(:,2),mesh.X(:,3),0);
p = sparse(opt.neqn_p,1);

Coldold = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  u1_oldold, u2_oldold, D_hat, mesh.X, mesh.IX);
Cold    = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  u1_old,    u2_old,    D_hat, mesh.X, mesh.IX);
C       = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  u1,        u2,        D_hat, mesh.X, mesh.IX);

% [L,Up,P] = lu(opt.LHS_BC);
% warning('off','all')

free_dofs = diag(opt.Null);

n1 = zeros(length(mesh.X),1);
n2 = zeros(length(mesh.X),1);

% South
South = mesh.X(mesh.X(:,3)==min(mesh.X(:,3)),1);
n2(South) = -1;
% East
East = mesh.X(mesh.X(:,2)==max(mesh.X(:,2)),1);
n1(East) = 1;
% North
North = mesh.X(mesh.X(:,3)==max(mesh.X(:,3)),1);
n2(North) = 1;
% West
West = mesh.X(mesh.X(:,2)==min(mesh.X(:,2)),1);
n1(West) = -1;

for n = 1:length(study.t)

    v_hat1 =  alpha0 * u1 + alpha1 * u1_old + alpha2 * u1_oldold + dt * (beta0 * C * u1 + beta1 * Cold * u1_old + beta2 * Coldold * u1_oldold + opt.f_u1);
    v_hat2 =  alpha0 * u2 + alpha1 * u2_old + alpha2 * u2_oldold + dt * (beta0 * C * u2 + beta1 * Cold * u2_old + beta2 * Coldold * u2_oldold + opt.f_u2);

    [curl1, curl2] =               calculateCurl(u1,        u2,        opt.grad1, opt.grad2, opt.gradgrad1, opt.gradgrad2);
    [curl1_old, curl2_old] =       calculateCurl(u1_old,    u2_old,    opt.grad1, opt.grad2, opt.gradgrad1, opt.gradgrad2);
    [curl1_oldold, curl2_oldold] = calculateCurl(u1_oldold, u2_oldold, opt.grad1, opt.grad2, opt.gradgrad1, opt.gradgrad2);

    nodes = mesh.bound_u1(:, 1);
    g_u1(nodes) = study.U1(mesh.X(nodes,2),mesh.X(nodes,3),study.t(n));
    g_u2(nodes) = study.U2(mesh.X(nodes,2),mesh.X(nodes,3),study.t(n));

    BC1 = 1/dt * (g_u1 - v_hat1) - (beta0 * curl1 + beta1 * curl1_old + beta2 * curl1_oldold) - opt.f_u1;
    BC2 = 1/dt * (g_u2 - v_hat2) - (beta0 * curl2 + beta1 * curl2_old + beta2 * curl2_oldold) - opt.f_u2;

    opt.p(:,n) = opt.A \ ( - opt.B * (n1 .* BC1 + n2 .* BC2) + (opt.D1 * v_hat1 + opt.D2 * v_hat2) ./ dt);
    opt.p(:,n) = opt.p(:,n) - mean(opt.p(:,n) - study.P(mesh.Xp(:,2),mesh.Xp(:,3),study.t(n)));

    v_hathat1 = v_hat1 - dt * D1 * opt.p(:,n);
    v_hathat2 = v_hat2 - dt * D2 * opt.p(:,n);

    opt.u1(:,n) = (1 - dt * gamma0 * H1) \ (dt * v_hathat1);
    opt.u2(:,n) = (1 - dt * gamma0 * H2) \ (dt * v_hathat2);


    % if length(mesh.bound_p) >= 1
    %     nodes = mesh.bound_p(:, 1);
    %     g_p(nodes) = study.P(mesh.Xp(nodes,2),mesh.Xp(nodes,3),study.t(n));
    % end

    % g = [g_u1; g_u2; g_p];
    %
    % RHS = RHS - opt.LHS * g;       % Adjust for stiffness matrix
    % RHS(~free_dofs) = 0;
    % RHS(find(g)) = g(find(g));   % Enforce boundary values on the load vector

    if(mod(n,100))==0
        clc
        fprintf('Time integration %i %c done\n',round(n/length(study.t)*100),'%');
    end

    if any(isnan([opt.u1(:,n); opt.u2(:,n); opt.p(:,n)]))
        warning('The vector contains NaN values.');
        return
    end

    Coldold = Cold;
    Cold = C;
    C = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  full(opt.u1(:,n)), full(opt.u2(:,n)), D_hat, mesh.X, mesh.IX);
    u1_oldold = u1_old;
    u1_old = u1;
    u1 = opt.u1(:,n);
    u2_oldold = u2_old;
    u2_old = u2;
    u2 = opt.u2(:,n);

end

end
