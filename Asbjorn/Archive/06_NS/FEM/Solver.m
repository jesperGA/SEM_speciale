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

S = D1 * (H1 \ D1') + D2 * (H2 \ D2');
S = 1/2*(S+S');

dt = study.t(2) - study.t(1);

rho = mesh.Material(1);
mu = mesh.Material(2);

opt.u1 = zeros(opt.neqn_u,length(study.t));
opt.u2 = zeros(opt.neqn_u,length(study.t));
opt.p  = zeros(opt.neqn_p,length(study.t));
g_u1 = sparse(opt.neqn_u,1);
g_u2 = sparse(opt.neqn_u,1);
g_p = sparse(opt.neqn_p,1);

MAXIT = 1e2;
tol = 1e-10;

E = opt.D1 * (opt.B \ opt.D1') + opt.D2 * (opt.B \ opt.D2'); %OBS w/o BCs enforced
preCond = (mu * inv(opt.M) + rho/dt * inv(E)); % Lav evt. noget chol() el. lign.

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

% u1_oldold=u1_oldold.*0;u1_old=u1_old.*0; %%%% SLETT

u1 = study.U1(mesh.X(:,2),mesh.X(:,3),0);
u2_oldold = study.U2(mesh.X(:,2),mesh.X(:,3),0-2*dt);
u2_old = study.U2(mesh.X(:,2),mesh.X(:,3),0-dt);
u2 = study.U2(mesh.X(:,2),mesh.X(:,3),0);
p = sparse(opt.neqn_p,1);

Coldold = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  u1_oldold, u2_oldold, D_hat, mesh.X, mesh.IX);
Cold    = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  u1_old,    u2_old,    D_hat, mesh.X, mesh.IX);
C       = globalCmatr_mex(ldof, opt.nel, opt.neqn_u, opt.Be,  u1,        u2,        D_hat, mesh.X, mesh.IX);

[L,Up,P] = lu(opt.LHS_BC);
warning('off','all')

free_dofs = diag(opt.Null);

for n = 1:length(study.t)

    g1 = opt.B * (opt.f_u1 + rho / dt * (alpha0 * u1 + alpha1 * u1_old + alpha2 * u1_oldold)) - study.Re * (beta0 * C * u1 + beta1 * Cold * u1_old + beta2 * Coldold * u1_oldold);
    g2 = opt.B * (opt.f_u2 + rho / dt * (alpha0 * u2 + alpha1 * u2_old + alpha2 * u2_oldold)) - study.Re * (beta0 * C * u2 + beta1 * Cold * u2_old + beta2 * Coldold * u2_oldold);

    RHS = [g1; g2; sparse(opt.neqn_p,1)];

    nodes = mesh.bound_u1(:, 1);
    g_u1(nodes) = study.U1(mesh.X(nodes,2),mesh.X(nodes,3),study.t(n));
    g_u2(nodes) = study.U2(mesh.X(nodes,2),mesh.X(nodes,3),study.t(n));

    if length(mesh.bound_p) >= 1
        nodes = mesh.bound_p(:, 1);
        g_p(nodes) = study.P(mesh.Xp(nodes,2),mesh.Xp(nodes,3),study.t(n));
    end

    g = [g_u1; g_u2; g_p];

    RHS = RHS - opt.LHS * g;       % Adjust for stiffness matrix
    RHS(~free_dofs) = 0;
    RHS(find(g)) = g(find(g));   % Enforce boundary values on the load vector

    if strcmp(study.solver,'Direct')
            
        % spparms('spumoni',2)
        temp = L \ (P*RHS);
        U = (Up \ temp);
        % U = opt.LHS_BC\RHS;
        % fprintf('Solution found in %f s\n',time);

        opt.u1(:,n) = U(1:opt.neqn_u);
        opt.u2(:,n) = U(opt.neqn_u+1:2*opt.neqn_u);
        opt.p(:,n)= U(2*opt.neqn_u+1:end);
        opt.p(:,n) = opt.p(:,n) - mean(opt.p(:,n) - study.P(mesh.Xp(:,2),mesh.Xp(:,3),study.t(n)));
    else
        RHS1 = RHS(1:opt.neqn_u);
        RHS2 = RHS(opt.neqn_u+1:2*opt.neqn_u);
        RHS3 = RHS(2*opt.neqn_u+1:end);
        RHS_Uzawa = - D1 * (H1 \ RHS1) - D2 * (H2 \ RHS2) - RHS3;
        if strcmp(study.solver,'pcg')
            if n==2;
                preCond = inv(preCond);
            end
            opt.p(:,n) = pcg(S, RHS_Uzawa, tol, MAXIT, preCond, [] , p);
    
        elseif strcmp(study.solver,'Uzawa')
    
            [opt.p(:,n),conv_message] = pcgUzawa(S, RHS_Uzawa, preCond,  p, H1, H2, D1, D2, tol, MAXIT);
    
        elseif strcmp(study.solver,'UzawamodJGA')
    
            opt.p(:,n) = pcg_mod(S, RHS_Uzawa, preCond, opt.p(:,n-1), H1, H2, D1, D2, tol, MAXIT);
    
        end

        opt.p(:,n) = opt.p(:,n) - mean(opt.p(:,n) - study.P(mesh.Xp(:,2),mesh.Xp(:,3),study.t(n)));

        opt.u1(:,n) = H1 \ (D1' * opt.p(:,n) + RHS1);
        opt.u2(:,n) = H2 \ (D2' * opt.p(:,n) + RHS2);
    end

    if(mod(n,100))==0
        clc
        if strcmp(study.solver,'Uzawa')
            disp(conv_message)
        end
        fprintf('Time integration %i %c done\n',round(n/length(study.t)*100),'%');
    end


    if any(isnan([opt.u1(:,n); opt.u2(:,n); opt.p(:,n)]))
        error('The vector contains NaN values.');
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
    p = opt.p(:,n);

end
        RHS1 = RHS(1:opt.neqn_u);
        RHS2 = RHS(opt.neqn_u+1:2*opt.neqn_u);
        RHS3 = RHS(2*opt.neqn_u+1:end);
% sum(abs(H1*opt.u1(:,n)-D1'*opt.p(:,n)-RHS1))
% sum(abs(H2*opt.u2(:,n)-D2'*opt.p(:,n)-RHS2))
sum(abs(-D1*opt.u1(:,n)-D2*opt.u2(:,n)-RHS3))

end