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

U_big = zeros(opt.neqn,length(study.t));
g_u1 = sparse(opt.neqn_u,1);
g_u2 = sparse(opt.neqn_u,1);
g_p = sparse(opt.neqn_p,1);

if strcmp(study.timeInt,'BDF1_AB3')
    alpha0 = 1;
    alpha1 = 0;
    alpha2 = 0;
    beta0 = 23/12;
    beta1 = -4/3;
    beta2 = 5/12;
elseif strcmp(study.timeInt,'BDF3_EX3')
    alpha0 = 3;
    alpha1 = -3/2;
    alpha2 = 1/3;
    beta0 = 3;
    beta1 = -3;
    beta2 = 1;
end


ldof = (study.N + 1) ^ 2;
    
[~,~,D_hat] = GetGLL(study.N+1);

u1_oldold = study.U1(mesh.X(:,2),mesh.X(:,3),0-2*dt);
u1_old = study.U1(mesh.X(:,2),mesh.X(:,3),0-dt);
u1 = study.U1(mesh.X(:,2),mesh.X(:,3),0);
u2_oldold = study.U2(mesh.X(:,2),mesh.X(:,3),0-2*dt);
u2_old = study.U2(mesh.X(:,2),mesh.X(:,3),0-dt);
u2 = study.U2(mesh.X(:,2),mesh.X(:,3),0);
p = sparse(opt.neqn_p,1);

Coldold = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  u1_oldold, u2_oldold, D_hat, mesh.L1, mesh.L2, mesh.IX);
Cold    = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  u1_old,    u2_old,    D_hat, mesh.L1, mesh.L2, mesh.IX);
C       = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  u1,        u2,        D_hat, mesh.L1, mesh.L2, mesh.IX);

% [L,Up,P] = lu(opt.LHS_BC);
Dsysmat = decomposition(opt.LHS_BC);

% warning('off','all')

free_dofs = diag(opt.Null);

prev_mes = 0;
fprintf('SOLVER \n ')

boundaryNodes = mesh.bound_u1(:, 1);

for n = 1:length(study.t)

    g1 = opt.B * (opt.f_u1 + rho / dt * (alpha0 * u1 + alpha1 * u1_old + alpha2 * u1_oldold)) - study.Re * (beta0 * C * u1 + beta1 * Cold * u1_old + beta2 * Coldold * u1_oldold);
    g2 = opt.B * (opt.f_u2 + rho / dt * (alpha0 * u2 + alpha1 * u2_old + alpha2 * u2_oldold)) - study.Re * (beta0 * C * u2 + beta1 * Cold * u2_old + beta2 * Coldold * u2_oldold);

    RHS = [g1; g2; sparse(opt.neqn_p,1)];

    g_u1(boundaryNodes) = study.U1(mesh.X(boundaryNodes,2),mesh.X(boundaryNodes,3),study.t(n));
    g_u1(1,1)=0;
    g_u1(311,1)=0;
    g_u2(boundaryNodes) = study.U2(mesh.X(boundaryNodes,2),mesh.X(boundaryNodes,3),study.t(n));
    g = [g_u1; g_u2; g_p];

    RHS = RHS - opt.LHS * g;       % Adjust for stiffness matrix
    RHS(~free_dofs) = 0;
    RHS(find(g)) = g(find(g));   % Enforce boundary values on the load vector

    % temp = L \ (P*RHS);
    % U = (Up \ temp);
    U = Dsysmat\RHS;

    if mod(n,100) == 0
        status_procent = (n/length(study.t))*100;
        mess = sprintf('%3.2f percent',status_procent);
        back_space = repmat('\b',1,prev_mes);
        fprintf(back_space);
        fprintf(mess);
        prev_mes = length(mess);
    end

    if any(isnan(U))
        error('The vector contains NaN values.');
    end

    u1_oldold = u1_old;
    u1_old = u1;
    u1 = U(1:opt.neqn_u);
    u2_oldold = u2_old;
    u2_old = u2;
    u2 = U(opt.neqn_u+1:2*opt.neqn_u);
    Coldold = Cold;
    Cold = C;
    C = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  full(u1), full(u2), D_hat, mesh.L1, mesh.L2, mesh.IX);
    p = U(2*opt.neqn_u+1:end);

    U_big(:,n) = [u1; u2; p];
end

opt.u1 = U_big(1:opt.neqn_u,:);
opt.u2 = U_big(opt.neqn_u+1:2*opt.neqn_u,:);
opt.p = U_big(2*opt.neqn_u+1:end,:);
opt.p = opt.p-mean(opt.p);


end